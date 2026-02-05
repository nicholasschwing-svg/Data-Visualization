function RawMultiBandViewer(initial)
% RawMultiBandViewer
% Multi-band RAW viewer + CERBERUS + MX20 HSI context viewer.
%
% Requires helper functions in the same folder:
%   fridge_init_from_raw.m
%   fridge_read_frame.m
%   fridge_derive_times_from_hdrs.m
%   viewer_format_pixel_value.m
%
% FRIDGE: LWIR/MWIR/SWIR/MONO/VIS-COLOR RAW stacks (ENVI BSQ, .raw + .hdr)
%   - Frames treated as time samples.
%   - Per-frame timestamps taken from ENVI "band names" when present.
%
% CERBERUS / MX20: calibrated ENVI-style cubes (.hsic + .hdr) shown as
% quicklook context images.

    if nargin < 1 || isempty(initial)
        initial = struct();
    end

    % Normalize map keys so callers can provide either char or string
    % modality names without triggering containers.Map indexing errors.
    labelCtor = pickLabelCtor();
    makeLabel = @(parent, varargin) labelCtor(parent, varargin{:});

    persistent lastWindowState;

    existingFig = findall(0, 'Type', 'figure', 'Tag', 'RawMultiBandViewer');
    if ~isempty(existingFig)
        upd = getappdata(existingFig(1), 'RMBV_Update');
        if isa(upd, 'function_handle')
            try
                upd(initial);
                figure(existingFig(1));
                return;
            catch
                % Fall back to rebuilding if the existing viewer cannot be
                % updated with the new selection.
            end
        end
        try
            lastWindowState = existingFig(1).WindowState;
        catch
            % Ignore missing property on older MATLAB versions.
        end
        delete(existingFig(1));
    end

    %======================== UI LAYOUT ===================================
    f = uifigure('Name','AARO Multi-Band Viewer', ...
        'Position',[80 80 1280 900], ...
        'Tag','RawMultiBandViewer');
    if ~isempty(lastWindowState) && isprop(f,'WindowState')
        try
            f.WindowState = lastWindowState;
        catch
            % Ignore window-state reapply failures.
        end
    end

    f.CloseRequestFcn = @(src,~)onCloseRequest(src);

    % 3x3 page grid (header, image grid, controls)
    page = uigridlayout(f,[3,3]);
    page.RowHeight   = {'fit', '1x', 'fit'};
    page.ColumnWidth = {'1x','1x','1x'};

    % Header row with title + return button
    headerRow = uigridlayout(page,[1,3]);
    headerRow.Layout.Row    = 1;
    headerRow.Layout.Column = [1 3];
    headerRow.ColumnWidth   = {'1x','fit','fit'};
    headerRow.Padding       = [8 8 8 8];

    header = makeLabel(headerRow, ...
        'Text','Multiband FRIDGE + HSI viewer (driven by timeline selection).', ...
        'FontWeight','bold','HorizontalAlignment','left');
    header.Layout.Row    = 1;
    header.Layout.Column = 1;

    makeLabel(headerRow,'Text','');  % spacer

    btnReturn = uibutton(headerRow, 'Text','Return to Timeline', ...
        'Enable','off', ...
        'ButtonPushedFcn',@(~,~)returnToTimeline());
    btnReturn.Layout.Row    = 1;
    btnReturn.Layout.Column = 3;
    btnReturn.FontWeight    = 'bold';

    % Axes names and grid positions
    modalities = {'LWIR','MWIR','SWIR','MONO','VIS-COLOR'};
    modalities = normalizeModalities(modalities);
    axMap          = containers.Map('KeyType','char','ValueType','any');
    frameLabelMap  = containers.Map('KeyType','char','ValueType','any');
    fileLabelMap   = containers.Map('KeyType','char','ValueType','any');
    panelMap       = containers.Map('KeyType','char','ValueType','any');
    hsiControlMap  = containers.Map('KeyType','char','ValueType','any');

    function ctor = pickLabelCtor()
        % Prefer uilabel when available, otherwise use the class constructor
        % with explicit Parent assignment, with a final fallback to uicontrol
        % text for legacy runtimes.
        if exist('uilabel','builtin') || exist('uilabel','file')
            ctor = @uilabel;
            return;
        end
        if exist('matlab.ui.control.Label','class')
            ctor = @(parent, varargin) matlab.ui.control.Label('Parent', parent, varargin{:});
            return;
        end
        ctor = @fallbackLabel;
    end

    function lbl = fallbackLabel(parent, varargin)
        lbl = uicontrol('Parent', parent, 'Style','text');
        for ii = 1:2:numel(varargin)
            name = varargin{ii};
            if ii+1 <= numel(varargin)
                val = varargin{ii+1};
            else
                break;
            end
            try
                set(lbl, name, val);
            catch
                % Best-effort application of properties; ignore ones that
                % uicontrol text labels do not support.
            end
        end
    end

    function kOut = keyify(kIn)
        % Normalize any provided key (char, string, cellstr) to a single
        % character vector so containers.Map indexing never sees a
        % multi-element input that would trigger "Only one level of
        % indexing is supported by a containers.Map." errors.
        if isstring(kIn)
            if numel(kIn) >= 1
                kOut = char(kIn(1));
            else
                kOut = '';
            end
        elseif iscellstr(kIn) || (iscell(kIn) && numel(kIn)==1 && ischar(kIn{1}))
            kOut = kIn{1};
        else
            kOut = kIn;
        end
    end

    function modsOut = normalizeModalities(modsIn)
        if isstring(modsIn)
            modsOut = cellstr(modsIn);
        else
            modsOut = modsIn;
        end
        for jj = 1:numel(modsOut)
            modsOut{jj} = keyify(modsOut{jj});
        end
    end

    function sz = parseResolution(str, defaultSz)
        tok = regexp(str, '^\s*(\d+)\s*[xX]\s*(\d+)\s*$', 'tokens', 'once');
        if isempty(tok)
            sz = defaultSz;
        else
            sz = [str2double(tok{1}), str2double(tok{2})];
            if any(isnan(sz) | sz <= 0)
                sz = defaultSz;
            end
        end
    end

    % Safe map accessors to avoid "key not present" runtime popups when
    % callers provide unexpected modality spellings.
    function tf = hasKey(mapObj, k)
        k = keyify(k);
        try
            tf = isKey(mapObj, k);
        catch ME
            if contains(ME.message, 'Only one level of indexing')
                tf = false;
            else
                rethrow(ME);
            end
        end
    end

    function v = getOr(mapObj, k, defaultVal)
        k = keyify(k);
        if nargin < 3
            defaultVal = [];
        end
        try
            if isKey(mapObj, k)
                v = mapObj(k);
            else
                v = defaultVal;
            end
        catch ME
            if contains(ME.message, 'Only one level of indexing')
                v = defaultVal;
            else
                rethrow(ME);
            end
        end
    end

    function mapOut = normalizeMapKeys(mapIn)
        % Rebuild maps with char keys so callers can supply string keys
        % without triggering "key not present" errors later.
        if isempty(mapIn) || ~isa(mapIn, 'containers.Map')
            mapOut = mapIn;
            return;
        end

        mapOut = containers.Map('KeyType','char','ValueType','any');
        keysIn = mapIn.keys;
        for mm = 1:numel(keysIn)
            kIn = keysIn{mm};
            kOut = keyify(kIn);
            mapOut(kOut) = mapIn(kIn);
        end
    end

    function mapOut = ensureMapKeys(mapIn, keyList, defaultVal)
        % Ensure a map contains the expected modality keys so downstream
        % lookups never trigger missing-key errors even if upstream helpers
        % used alternate spellings or omitted entries.
        mapOut = containers.Map('KeyType','char','ValueType','any');
        for ii = 1:numel(keyList)
            k = keyify(keyList{ii});
            if isa(mapIn,'containers.Map') && isKey(mapIn, k)
                mapOut(k) = mapIn(k);
            else
                mapOut(k) = defaultVal;
            end
        end
    end

    % Image grid that is rebuilt based on available panes
    imgGrid = uigridlayout(page,[1,1]);
    imgGrid.Layout.Row    = 2;
    imgGrid.Layout.Column = [1 3];
    imgGrid.RowHeight     = {'1x'};
    imgGrid.ColumnWidth   = {'1x'};

    % Off-screen parent used to hold inactive panels so the grid only sees
    % panes that are actually visible for the current selection.
    hiddenBin = uipanel(f, 'Visible','off', 'Position',[0 0 1 1]);

    %----------------------------------------------------------------------
    % FRIDGE panes: panel + inner grid (label row + axes row)
    %----------------------------------------------------------------------
    for i = 1:numel(modalities)
        pnl = uipanel(hiddenBin);

        pGrid = uigridlayout(pnl,[4,1]);
        pGrid.RowHeight   = {'fit','fit','1x','fit'};
        pGrid.ColumnWidth = {'1x'};

        modName = keyify(modalities{i});
        if strcmp(modName,'VIS-COLOR')
            dispName = 'FRIDGE VIS';
        else
            dispName = ['FRIDGE ' modName];
        end

        lblTop = makeLabel(pGrid, ...
            'Text', dispName, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');
        lblTop.Layout.Row    = 1;
        lblTop.Layout.Column = 1;

        lblFrame = makeLabel(pGrid, ...
            'Text', 'Frame: -', ...
            'HorizontalAlignment','center', ...
            'FontWeight','bold');
        lblFrame.Layout.Row    = 2;
        lblFrame.Layout.Column = 1;

        ax = uiaxes(pGrid);
        ax.Layout.Row    = 3;
        ax.Layout.Column = 1;
        axis(ax,'off');
        title(ax, '');

        lblFilePane = makeLabel(pGrid, ...
            'Text','', ...
            'Interpreter','none', ...
            'HorizontalAlignment','center');
        lblFilePane.Layout.Row    = 4;
        lblFilePane.Layout.Column = 1;

        axMap(modName) = ax;
        frameLabelMap(modName) = lblFrame;
        fileLabelMap(modName)  = lblFilePane;
        panelMap(modName)      = pnl;
    end

    % HSI tab group (CERB LWIR/VNIR + MX20)
    hsiTabStash = uitabgroup(hiddenBin, 'Visible','off');
    hsiPanel = uipanel(hiddenBin);
    hsiGrid = uigridlayout(hsiPanel,[2,1]);
    hsiGrid.RowHeight   = {'fit','1x'};
    hsiGrid.ColumnWidth = {'1x'};

    lblHSIContext = makeLabel(hsiGrid, ...
        'Text','HSI Context', ...
        'FontWeight','bold', ...
        'HorizontalAlignment','center');
    lblHSIContext.Layout.Row    = 1;
    lblHSIContext.Layout.Column = 1;

    cerbTabs = uitabgroup(hsiGrid);
    cerbTabs.Layout.Row    = 2;
    cerbTabs.Layout.Column = 1;

    tabLWIR = uitab(cerbTabs,'Title','CERB LWIR');
    tabVNIR = uitab(cerbTabs,'Title','CERB VNIR');
    tabMX20 = uitab(cerbTabs,'Title','MX20 SW');
    fastTabs = containers.Map('KeyType','char','ValueType','any');
    fastAxes = containers.Map('KeyType','char','ValueType','any');
    tabHSIPlaceholder = uitab(cerbTabs,'Title','HSI Unavailable');

    placeholderGrid = uigridlayout(tabHSIPlaceholder,[1 1]);
    placeholderGrid.RowHeight   = {'1x'};
    placeholderGrid.ColumnWidth = {'1x'};

    placeholderLabel = makeLabel(placeholderGrid, ...
        'Text','No HSI context available for this selection.', ...
        'HorizontalAlignment','center', ...
        'FontAngle','italic');
    placeholderLabel.Layout.Row    = 1;
    placeholderLabel.Layout.Column = 1;
    
    % --- CERB LWIR tab: control row + axes ---
    tabLWIRGrid = uigridlayout(tabLWIR,[2 1]);
    tabLWIRGrid.RowHeight   = {'fit','1x'};
    tabLWIRGrid.ColumnWidth = {'1x'};

    createHsiControlRow(tabLWIRGrid, 'CERB_LWIR', 'CERBERUS LWIR');

    cerbAxLWIR = uiaxes(tabLWIRGrid);
    cerbAxLWIR.Layout.Row    = 2;
    cerbAxLWIR.Layout.Column = 1;
    axis(cerbAxLWIR,'off');
    title(cerbAxLWIR,'CERB LWIR');

    % --- CERB VNIR tab ---
    tabVNIRGrid = uigridlayout(tabVNIR,[2 1]);
    tabVNIRGrid.RowHeight   = {'fit','1x'};
    tabVNIRGrid.ColumnWidth = {'1x'};

    createHsiControlRow(tabVNIRGrid, 'CERB_VNIR', 'CERBERUS VNIR');

    cerbAxVNIR = uiaxes(tabVNIRGrid);
    cerbAxVNIR.Layout.Row    = 2;
    cerbAxVNIR.Layout.Column = 1;
    axis(cerbAxVNIR,'off');
    title(cerbAxVNIR,'CERB VNIR');

    % --- MX20 tab ---
    tabMX20Grid = uigridlayout(tabMX20,[2 1]);
    tabMX20Grid.RowHeight   = {'fit','1x'};
    tabMX20Grid.ColumnWidth = {'1x'};

    createHsiControlRow(tabMX20Grid, 'MX20', 'MX20 SW');

    mxAx = uiaxes(tabMX20Grid);
    mxAx.Layout.Row    = 2;
    mxAx.Layout.Column = 1;
    axis(mxAx,'off');
    title(mxAx,'MX20 SW');

    function ax = ensureFastAxis(modality)
        key = upper(modality);
        if ~isKey(fastAxes, key) || isempty(fastAxes(key)) || ~isgraphics(fastAxes(key))
            [tab, axNew] = createFastTab(key);
            fastTabs(key) = tab;
            fastAxes(key) = axNew;
        end
        ax = fastAxes(key);
    end

    function [tab, ax] = createFastTab(modality)
        key = upper(modality);
        tab = uitab(cerbTabs, 'Title', sprintf('FAST %s', key));
        grid = uigridlayout(tab,[2 1]);
        grid.RowHeight   = {'fit','1x'};
        grid.ColumnWidth = {'1x'};

        sensorKey = sprintf('FAST_%s', key);
        createHsiControlRow(grid, sensorKey, sprintf('FAST %s', key));

        ax = uiaxes(grid);
        ax.Layout.Row    = 2;
        ax.Layout.Column = 1;
        axis(ax,'off');
        title(ax, sprintf('FAST %s', key));
    end

    function createHsiControlRow(parent, sensorKey, displayName)
        ctrlRow = uigridlayout(parent, [1 3]);
        ctrlRow.Layout.Row    = 1;
        ctrlRow.Layout.Column = 1;
        ctrlRow.ColumnWidth   = {'fit','1x','fit'};
        ctrlRow.RowHeight     = {'fit'};

        btnPrevScan = uibutton(ctrlRow, 'Text','Prev Scan', 'Enable','off', ...
            'ButtonPushedFcn', @(~,~)stepScan(sensorKey, -1));
        btnPrevScan.Layout.Row    = 1;
        btnPrevScan.Layout.Column = 1;

        lblScan = makeLabel(ctrlRow, 'Text', sprintf('%s | Scan: -', displayName), ...
            'HorizontalAlignment','center');
        lblScan.Layout.Row    = 1;
        lblScan.Layout.Column = 2;

        btnNextScan = uibutton(ctrlRow, 'Text','Next Scan', 'Enable','off', ...
            'ButtonPushedFcn', @(~,~)stepScan(sensorKey, 1));
        btnNextScan.Layout.Row    = 1;
        btnNextScan.Layout.Column = 3;

        hsiControlMap(sensorKey) = struct('btnPrev', btnPrevScan, 'btnNext', btnNextScan, ...
            'label', lblScan, 'displayName', displayName);
    end

    % Create common FAST tabs up front
    ensureFastAxis('LWIR');
    ensureFastAxis('VNIR');

    panelMap('HSI') = hsiPanel;

    %----------------------------------------------------------------------
    % Controls row (bottom) – spread controls across three columns so the
    % slider stays wide and the timestamp is always visible.
    %----------------------------------------------------------------------
    ctrlWrapper = uigridlayout(page,[1,3]);
    ctrlWrapper.Layout.Row    = 3;
    ctrlWrapper.Layout.Column = [1 3];
    ctrlWrapper.RowHeight     = {'fit'};
    ctrlWrapper.ColumnWidth   = {'1.2x','2x','1x'};

    % Left column: file/memory + pixel readout
    infoCol = uigridlayout(ctrlWrapper,[2,1]);
    infoCol.RowHeight   = {'fit','fit'};
    infoCol.ColumnWidth = {'1x'};

    statusRow = uigridlayout(infoCol,[1,3]);
    statusRow.ColumnWidth = {'1x','fit','fit'};
    lblStatus = makeLabel(statusRow,'Text','Status: (no capture loaded)','HorizontalAlignment','left');
    lblFrames = makeLabel(statusRow,'Text','Frames: -');
    lblMem    = makeLabel(statusRow,'Text','');

    pixelRow = uigridlayout(infoCol,[1,2]);
    pixelRow.ColumnWidth = {'1x','1x'};
    lblPixel = makeLabel(pixelRow, ...
        'Text','Pixel: -', ...
        'HorizontalAlignment','left');
    lblValue = makeLabel(pixelRow, ...
        'Text','Value: -', ...
        'HorizontalAlignment','left');

    % Center column: navigation + slider kept wide
    navCol = uigridlayout(ctrlWrapper,[2,1]);
    navCol.RowHeight   = {'fit','fit'};
    navCol.ColumnWidth = {'1x'};

    navTop = uigridlayout(navCol,[1,2]);
    navTop.ColumnWidth = {'fit','fit'};
    btnSnapshot = uibutton(navTop,'Text','Save Snapshot','Enable','off', ...
        'ButtonPushedFcn',@(~,~)saveSnapshot());
    btnExport = uibutton(navTop,'Text','Export Video...','Enable','off', ...
        'ButtonPushedFcn',@(~,~)exportVideo());

    navBottom = uigridlayout(navCol,[1,7]);
    navBottom.ColumnWidth = {'fit','fit','fit','1x','fit','fit','fit'};
    btnPrev = uibutton(navBottom,'Text','Previous','Enable','off', ...
        'ButtonPushedFcn',@(~,~)step(-1));
    btnNext = uibutton(navBottom,'Text','Next','Enable','off', ...
        'ButtonPushedFcn',@(~,~)step(+1));
    makeLabel(navBottom,'Text','Time slider:','HorizontalAlignment','right');
    frameSlider = uislider(navBottom, ...
        'Limits',[1 2], ...
        'Value',1, ...
        'MajorTicks',[], ...
        'MinorTicks',[], ...
        'Enable','off', ...
        'ValueChangingFcn',@frameSliderChanging, ...
        'ValueChangedFcn',@frameSliderChanged);
    % Spacer to keep right-side items from crowding the slider
    makeLabel(navBottom,'Text','');
    btnJumpHsi = uibutton(navBottom,'Text','Jump to HSI','Enable','off', ...
        'Tooltip','Align FRIDGE to the current HSI timestamp', ...
        'ButtonPushedFcn',@(~,~)jumpToHsi());
    % Timestamp label sits in the right column

    % Right column: dedicated timestamp display
    timeCol = uigridlayout(ctrlWrapper,[4,1]);
    timeCol.RowHeight   = {'fit','fit','fit','fit'};
    timeCol.ColumnWidth = {'1x'};
    makeLabel(timeCol,'Text','Current time','FontWeight','bold', ...
        'HorizontalAlignment','left');
    lblTime = makeLabel(timeCol,'Text','Time: -','HorizontalAlignment','left');
    lblSelectionRange = makeLabel(timeCol,'Text','Selection: -', ...
        'HorizontalAlignment','left');
    lblDataSpan = makeLabel(timeCol,'Text','Data span: -', ...
        'HorizontalAlignment','left');

    %======================== STATE =======================================
    S = struct();
    S.files        = containers.Map(modalities, repmat({''},1,numel(modalities)));
    S.hdrs         = containers.Map(modalities, repmat({[]},1,numel(modalities)));
    S.exists       = containers.Map(modalities, num2cell(false(1,numel(modalities))));
    S.maxFrames    = containers.Map(modalities, num2cell(nan(1,numel(modalities))));
    S.frame        = 1;
    S.nFrames      = 0;            % timeline steps (time-ordered)
    S.frameCount   = 0;            % max frames across modalities
    S.dir          = '';
    S.chosen       = '';
    S.behavior     = 'hold';   % hold last frame when beyond range
    S.fridgeTimes  = [];       % legacy datetime vector
    S.fridgeTimesMap = containers.Map(modalities, repmat({datetime.empty(0,1)},1,numel(modalities)));
    S.timelineTimes  = datetime.empty(0,1);   % union of all modality times
    S.sliderMode   = 'frame';  % 'frame' (fallback) or 'time'
    S.sliderOrigin = NaT;      % reference time for slider (time mode)
    S.tStart       = NaT;
    S.tEnd         = NaT;
    S.tNow         = NaT;
    S.sliderStartTime = NaT;
    S.sliderEndTime   = NaT;
    S.dataStartTime   = NaT;
    S.dataEndTime     = NaT;
    S.sliderRangeSec = NaN;
    S.sliderStepSec  = NaN;
    S.lastFrameByMod = containers.Map(modalities, num2cell(nan(1,numel(modalities))));
    S.lastStatusByMod = containers.Map(modalities, repmat({''},1,numel(modalities)));
    S.fridgeStartTime = NaT;
    S.fridgeEndTime   = NaT;
    S.hsiEvents    = struct('sensor', {}, 'time', {}, 'path', {}, 'modality', {});
    S.hsiGroupsMap = containers.Map('KeyType','char','ValueType','any');
    S.currentHsi   = struct('sensor','', 'modality','', 'time', NaT, 'effectiveTime', NaT);
    S.hsiPreciseCache = containers.Map('KeyType','char','ValueType','any');
    S.timelineFig  = [];

    S.cerb = struct('LWIR',[],'VNIR',[]);
    S.mx20 = struct('hdr',[],'ctx',[]);
    S.fast = struct();
    S.enableHSI = true;
    S.currentHsiMap = containers.Map('KeyType','char','ValueType','any');

    targetStartTime = [];
    enableHSI = true;
    sliderInternalUpdate = false;  % prevent recursive slider callbacks

    %======================== LAYOUT HELPERS ==============================
    function tf = hasCerb(whichMod)
        tf = isfield(S, 'cerb') && isfield(S.cerb, whichMod) && ...
             ~isempty(S.cerb.(whichMod));
    end

    function tf = hasMx20()
        tf = isfield(S, 'mx20') && isfield(S.mx20, 'hdr') && ...
             ~isempty(S.mx20.hdr);
    end

    function tf = hasFast(whichMod)
        tf = isfield(S, 'fast') && isfield(S.fast, whichMod) && ...
             ~isempty(S.fast.(whichMod));
    end

    function tf = hasAnyHsi()
        if ~S.enableHSI
            tf = false;
            return;
        end
        tf = hasCerb('LWIR') || hasCerb('VNIR') || hasMx20();
        if ~tf && isfield(S,'fast')
            mods = fieldnames(S.fast);
            for ii = 1:numel(mods)
                tf = tf || hasFast(mods{ii});
            end
        end
    end

    function [tStart, tEnd] = resolveTimeDomain(initStruct)
        tStart = NaT;
        tEnd   = NaT;

        if isfield(initStruct,'tStart') && isfield(initStruct,'tEnd')
            if isdatetime(initStruct.tStart) && isdatetime(initStruct.tEnd)
                tStart = initStruct.tStart;
                tEnd   = initStruct.tEnd;
            end
        end

        if isempty(tStart) || isnat(tStart) || isempty(tEnd) || isnat(tEnd) || tEnd < tStart
            if isfield(initStruct,'hsiEvents') && ~isempty(initStruct.hsiEvents)
                evtTimes = [initStruct.hsiEvents.time];
                evtTimes = evtTimes(~isnat(evtTimes));
                if ~isempty(evtTimes)
                    tStart = min(evtTimes);
                    tEnd   = max(evtTimes);
                end
            end
        end

        if isempty(tStart) || isnat(tStart) || isempty(tEnd) || isnat(tEnd) || tEnd < tStart
            if isfield(initStruct,'fridgeStartTime') && isfield(initStruct,'fridgeEndTime')
                if isdatetime(initStruct.fridgeStartTime) && isdatetime(initStruct.fridgeEndTime)
                    tStart = initStruct.fridgeStartTime;
                    tEnd   = initStruct.fridgeEndTime;
                end
            end
        end
    end

    function updateHsiTabVisibility()
        function moveTab(tabHandle, show)
            if show
                targetParent = cerbTabs;
            else
                targetParent = hsiTabStash;
            end
            if tabHandle.Parent ~= targetParent
                tabHandle.Parent = targetParent;
            end
        end

        moveTab(tabLWIR, hasCerb('LWIR'));
        moveTab(tabVNIR, hasCerb('VNIR'));
        moveTab(tabMX20, hasMx20());
        fastKeys = fastTabs.keys;
        hasFastAny = false;
        for ii = 1:numel(fastKeys)
            moveTab(fastTabs(fastKeys{ii}), hasFast(fastKeys{ii}));
            hasFastAny = hasFastAny || hasFast(fastKeys{ii});
        end
        moveTab(tabHSIPlaceholder, ~(hasCerb('LWIR') || hasCerb('VNIR') || hasMx20() || hasFastAny));

        available = cerbTabs.Children;
        preferred = {tabLWIR, tabVNIR, tabMX20};
        for ii = 1:numel(fastKeys)
            if fastTabs(fastKeys{ii}).Parent == cerbTabs
                preferred{end+1} = fastTabs(fastKeys{ii}); %#ok<AGROW>
            end
        end
        preferred{end+1} = tabHSIPlaceholder;
        selected = [];
        for ii = 1:numel(preferred)
            if preferred{ii}.Parent == cerbTabs
                selected = preferred{ii};
                break;
            end
        end
        if ~isempty(selected)
            cerbTabs.SelectedTab = selected;
        elseif ~isempty(available)
            cerbTabs.SelectedTab = available(1);
        end
    end

    function activePanels = getActivePanels()
        activePanels = {};
        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if getOr(S.exists, m, false)
                activePanels{end+1} = m; %#ok<AGROW>
            end
        end

        if hasAnyHsi()
            activePanels{end+1} = 'HSI';
        end
    end

    function refreshMontageLayout()
        updateHsiTabVisibility();

        keysAll = panelMap.keys;
        for kk = 1:numel(keysAll)
            pnl = panelMap(keysAll{kk});
            pnl.Parent  = hiddenBin;
            pnl.Visible = 'off';
        end

        activePanels = getActivePanels();
        n = numel(activePanels);
        if n == 0
            imgGrid.RowHeight   = {'1x'};
            imgGrid.ColumnWidth = {'1x'};
            return;
        end

        maxCols = 3;
        cols = min(maxCols, n);
        rows = ceil(n / cols);
        imgGrid.RowHeight   = repmat({'1x'}, 1, rows);
        imgGrid.ColumnWidth = repmat({'1x'}, 1, cols);

        for idx = 1:n
            key = keyify(activePanels{idx});
            pnl = panelMap(key);
            pnl.Parent = imgGrid;
            pnl.Layout.Row    = ceil(idx / cols);
            pnl.Layout.Column = mod(idx-1, cols) + 1;
            pnl.Visible = 'on';
        end
    end

    %======================== CORE CALLBACKS ===============================
    function loadFromRawFile(fullRawPath)
        if isstring(fullRawPath) && isscalar(fullRawPath)
            fullRawPath = char(fullRawPath);
        end
        existingEvents = S.hsiEvents;
        savedTimeDomain = struct('tStart', S.tStart, 'tEnd', S.tEnd, ...
            'tNow', S.tNow, 'fridgeStartTime', S.fridgeStartTime, ...
            'fridgeEndTime', S.fridgeEndTime, ...
            'sliderStartTime', S.sliderStartTime, 'sliderEndTime', S.sliderEndTime, ...
            'dataStartTime', S.dataStartTime, 'dataEndTime', S.dataEndTime);
        resetUI();
        S.hsiEvents = existingEvents;
        S.tStart = savedTimeDomain.tStart;
        S.tEnd   = savedTimeDomain.tEnd;
        S.tNow   = savedTimeDomain.tNow;
        S.fridgeStartTime = savedTimeDomain.fridgeStartTime;
        S.fridgeEndTime   = savedTimeDomain.fridgeEndTime;
        S.sliderStartTime = savedTimeDomain.sliderStartTime;
        S.sliderEndTime   = savedTimeDomain.sliderEndTime;
        S.dataStartTime   = savedTimeDomain.dataStartTime;
        S.dataEndTime     = savedTimeDomain.dataEndTime;
        rebuildHsiGroups();

        if ~isfile(fullRawPath)
            uialert(f, sprintf('RAW file not found:\n%s', fullRawPath), ...
                'File Error');
            return;
        end

        [path,file,ext] = fileparts(fullRawPath);
        if ~strcmpi(ext,'.raw')
            uialert(f, 'Please provide a .raw file.', 'Not RAW');
            return;
        end

        S.dir = [path filesep];

        % Parse prefix + main modality from filename
        [prefix, pickedMod] = parsePrefixAndModality([file ext]);  % external
        if isempty(pickedMod)
            uialert(f, ...
                'Filename must contain one of LWIR/MWIR/SWIR/MONO/VIS-COLOR before .raw', ...
                'Unrecognized');
            return;
        end
        S.chosen = pickedMod;

        % -------- FRIDGE initialization moved to helper ------------------
        [filesTmp, hdrsTmp, existsTmp, maxFramesTmp, ...
         nFramesTmp, fridgeTimesTmp, fridgeTimesMapTmp] = ...
            fridge_init_from_raw(path, prefix, modalities);

        % Rebuild the returned maps so all keys are char vectors and ensure
        % every expected modality is represented. This prevents downstream
        % containers.Map indexing errors that would otherwise surface as
        % UI popups instead of images.
        filesTmp       = normalizeMapKeys(filesTmp);
        hdrsTmp        = normalizeMapKeys(hdrsTmp);
        existsTmp      = normalizeMapKeys(existsTmp);
        maxFramesTmp   = normalizeMapKeys(maxFramesTmp);
        fridgeTimesMapTmp = normalizeMapKeys(fridgeTimesMapTmp);

        S.files         = ensureMapKeys(filesTmp, modalities, '');
        S.hdrs          = ensureMapKeys(hdrsTmp, modalities, []);
        S.exists        = ensureMapKeys(existsTmp, modalities, false);
        S.maxFrames     = ensureMapKeys(maxFramesTmp, modalities, NaN);
        S.fridgeTimesMap= ensureMapKeys(fridgeTimesMapTmp, modalities, datetime.empty(0,1));
        S.nFrames       = nFramesTmp;
        S.fridgeTimes   = fridgeTimesTmp;
        S.frameCount = S.nFrames;

        refreshMontageLayout();

        % If the timeline already passed per-frame timestamps, prefer them so
        % alignment works even when headers omit band_names.
        if isfield(initial, 'fridgeTimes') && ~isempty(initial.fridgeTimes)
            S.fridgeTimes = initial.fridgeTimes(:);
            for ii = 1:numel(modalities)
                key = keyify(modalities{ii});
                S.fridgeTimesMap(key) = S.fridgeTimes;
            end
        end
        if isempty(S.fridgeTimes) && ~isempty(S.fridgeTimesMap)
            % Legacy callers may have only populated the map
            mods = S.fridgeTimesMap.keys;
            for ii = 1:numel(mods)
                mKey = keyify(mods{ii});
                if ~isempty(S.fridgeTimesMap(mKey))
                    S.fridgeTimes = S.fridgeTimesMap(mKey);
                    break;
                end
            end
        end
        % -----------------------------------------------------------------

        lblStatus.Text = 'Status: FRIDGE filenames shown in panels';
        lblFrames.Text = sprintf('Frames: %d', S.frameCount);

        if exist('memory','file') == 2 || exist('memory','builtin') == 5
            [mem, ~] = memory();
            lblMem.Text = sprintf('Avail mem: %.0f MB', ...
                mem.MemAvailableAllArrays/1024/1024);
        else
            lblMem.Text = '';
        end

        btnPrev.Enable     = 'on';
        btnNext.Enable     = 'on';
        btnSnapshot.Enable = 'on';
        btnExport.Enable   = 'on';
        updateHsiJumpAvailability();

        rebuildTimeline();
        recomputeSliderDataWindow();
        configureTimeSlider();
        if ~isempty(targetStartTime) && (isnat(S.tNow) || isempty(S.tNow))
            updateAllPanesAtTime(targetStartTime, true);
        elseif ~isnat(S.tNow)
            updateAllPanesAtTime(S.tNow, true);
        end
    end

    function step(delta)
        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime)
            return;
        end
        if isnan(S.sliderStepSec) || S.sliderStepSec <= 0
            stepSec = 1;
        else
            stepSec = S.sliderStepSec;
        end
        tNow = S.tNow;
        if isempty(tNow) || isnat(tNow)
            tNow = S.sliderStartTime;
        end
        tTarget = tNow + seconds(delta * stepSec);
        updateAllPanesAtTime(tTarget, true);
    end

    function frameSliderChanging(~, evt)
        if sliderInternalUpdate
            return;
        end
        applySliderValue(evt.Value, false);
    end

    function frameSliderChanged(src, ~)
        if sliderInternalUpdate
            return;
        end
        applySliderValue(src.Value, true);
    end

    function tf = hasFridgeTimes()
        tf = ~isempty(S.timelineTimes) && isdatetime(S.timelineTimes) && ...
             all(~isnat(S.timelineTimes));
    end

    function t = timeForFrame(idx)
        t = [];
        if ~hasFridgeTimes()
            return;
        end
        idx = min(max(1, idx), numel(S.timelineTimes));
        t = S.timelineTimes(idx);
    end

    function idx = frameForTime(tTarget)
        if hasFridgeTimes()
            [~, idx] = min(abs(S.timelineTimes - tTarget));
        else
            idx = NaN;
        end
        if isempty(idx) || isnan(idx)
            return;
        end
        idx = min(max(1, idx), S.nFrames);
    end

    function rebuildTimeline()
        % Build a combined FRIDGE time axis for export/time lookup.
        S.timelineTimes = datetime.empty(0,1);
        S.frame = 1;

        allTimes = datetime.empty(0,1);
        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if ~isKey(S.fridgeTimesMap, m)
                continue;
            end
            tVec = S.fridgeTimesMap(m);
            if ~isdatetime(tVec)
                continue;
            end
            tVec = tVec(~isnat(tVec));
            if isempty(tVec)
                continue;
            end
            allTimes = [allTimes; tVec(:)]; %#ok<AGROW>
        end

        if isempty(allTimes) && ...
                ~isnat(S.fridgeStartTime) && ...
                ~isnat(S.fridgeEndTime)
            % Synthesize evenly spaced times when only capture bounds exist.
            durSec = seconds(S.fridgeEndTime - S.fridgeStartTime);
            steps  = max(2, max(1, S.frameCount));
            offsets = linspace(0, durSec, steps);
            allTimes = S.fridgeStartTime + seconds(offsets(:));
        end

        if ~isempty(allTimes)
            allTimes        = unique(allTimes);
            allTimes        = sort(allTimes);
            S.timelineTimes = allTimes;
            S.nFrames       = numel(S.timelineTimes);
        else
            S.nFrames       = max(1, S.frameCount);
        end
    end

    function jumpToTime(tTarget)
        if isempty(tTarget)
            return;
        end
        updateAllPanesAtTime(tTarget, true);
    end

    function applySliderValue(val, isFinal)
        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime)
            return;
        end
        tTarget = S.sliderStartTime + seconds(val);
        updateAllPanesAtTime(tTarget, isFinal);
    end

    function tOut = clampTime(tIn)
        tOut = tIn;
        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime) || isempty(tIn) || isnat(tIn)
            return;
        end
        if tOut < S.sliderStartTime
            tOut = S.sliderStartTime;
        elseif tOut > S.sliderEndTime
            tOut = S.sliderEndTime;
        end
    end

    function recomputeSliderDataWindow()
        S.sliderStartTime = S.tStart;
        S.sliderEndTime   = S.tEnd;
        S.dataStartTime   = NaT;
        S.dataEndTime     = NaT;

        if isnat(S.tStart) || isnat(S.tEnd)
            return;
        end

        tCollected = datetime.empty(0,1);

        % FRIDGE timestamps
        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if ~getOr(S.exists, m, false) || ~isKey(S.fridgeTimesMap, m)
                continue;
            end
            tVec = S.fridgeTimesMap(m);
            if ~isdatetime(tVec) || isempty(tVec)
                continue;
            end
            tVec = tVec(~isnat(tVec));
            tVec = tVec(tVec >= S.tStart & tVec <= S.tEnd);
            if ~isempty(tVec)
                tCollected = [tCollected; tVec(:)]; %#ok<AGROW>
            end
        end

        % HSI group timestamps
        if S.enableHSI && ~isempty(S.hsiGroupsMap)
            hsiKeys = S.hsiGroupsMap.keys;
            for kk = 1:numel(hsiKeys)
                data = S.hsiGroupsMap(hsiKeys{kk});
                if ~isfield(data,'timesUnique') || isempty(data.timesUnique)
                    continue;
                end
                tVec = data.timesUnique(:);
                tVec = tVec(~isnat(tVec));
                tVec = tVec(tVec >= S.tStart & tVec <= S.tEnd);
                if ~isempty(tVec)
                    tCollected = [tCollected; tVec]; %#ok<AGROW>
                end
            end
        end

        if isempty(tCollected)
            return;
        end

        S.dataStartTime = min(tCollected);
        S.dataEndTime   = max(tCollected);
        S.sliderStartTime = max(S.tStart, S.dataStartTime);
        S.sliderEndTime   = min(S.tEnd,   S.dataEndTime);

        if S.sliderEndTime < S.sliderStartTime
            S.sliderStartTime = S.tStart;
            S.sliderEndTime   = S.tEnd;
        end
    end

    function configureTimeSlider()
        if isnat(S.tStart) || isnat(S.tEnd)
            frameSlider.Enable = 'off';
            frameSlider.Limits = [0 1];
            frameSlider.Value  = 0;
            lblSelectionRange.Text = 'Selection: -';
            lblDataSpan.Text = 'Data span: -';
            return;
        end

        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime)
            S.sliderStartTime = S.tStart;
            S.sliderEndTime   = S.tEnd;
        end

        rangeSec = seconds(S.sliderEndTime - S.sliderStartTime);
        if rangeSec <= 0 || isnan(rangeSec)
            rangeSec = 1;
        end
        S.sliderRangeSec = rangeSec;

        smallStepSec = min(max(0.5, rangeSec * 0.01), 1.0);
        largeStepSec = min(max(5.0, rangeSec * 0.1), 10.0);
        smallStep = min(max(smallStepSec / rangeSec, 1e-4), 0.1);
        largeStep = min(max(largeStepSec / rangeSec, 1e-3), 0.5);
        S.sliderStepSec = smallStepSec;

        frameSlider.Limits = [0 rangeSec];
        % Older MATLAB releases may not expose SliderStep on uislider.
        if isprop(frameSlider, 'SliderStep')
            frameSlider.SliderStep = [smallStep largeStep];
        end
        frameSlider.Value = 0;
        frameSlider.Enable = 'on';

        lblSelectionRange.Text = sprintf('Selection: %s', formatClockRange(S.tStart, S.tEnd));
        if ~isnat(S.dataStartTime) && ~isnat(S.dataEndTime)
            lblDataSpan.Text = sprintf('Data span: %s', formatClockRange(S.dataStartTime, S.dataEndTime));
        else
            lblDataSpan.Text = 'Data span: (no data in selection)';
        end
    end

    function updateSliderValueFromTime(tNow)
        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime) || isempty(tNow) || isnat(tNow)
            return;
        end
        sliderInternalUpdate = true;
        frameSlider.Value = seconds(tNow - S.sliderStartTime);
        sliderInternalUpdate = false;
    end

    function updateAllPanesAtTime(tNow, isFinal)
        if nargin < 2
            isFinal = true;
        end
        if isempty(tNow) || isnat(tNow) || isnat(S.sliderStartTime) || isnat(S.sliderEndTime)
            return;
        end
        tNow = clampTime(tNow);
        S.tNow = tNow;
        updateSliderValueFromTime(tNow);

        if hasFridgeTimes()
            idx = frameForTime(tNow);
            if ~isempty(idx) && ~isnan(idx)
                S.frame = idx;
            end
        end

        drawAll(tNow);
        updateTimeDisplay();
        syncHsiToTime(tNow);
    end

    function [idxSel, status] = pickNearestIndex(times, tNow, maxFrames)
        idxSel = NaN;
        status = 'missing';
        if isempty(times) || isempty(tNow) || any(isnat(tNow))
            return;
        end

        times = times(:);
        [timesSorted, ord] = sort(times);
        if tNow <= timesSorted(1)
            idxSel = ord(1);
            status = 'held';
        elseif tNow >= timesSorted(end)
            idxSel = ord(end);
            status = 'held';
        else
            diffs = abs(timesSorted - tNow);
            minDiff = min(diffs);
            idxLocal = find(diffs == minDiff, 1, 'first');
            idxSel = ord(idxLocal);
            status = 'ok';
        end

        if nargin >= 3 && ~isempty(maxFrames) && ~isnan(maxFrames)
            idxSel = min(idxSel, maxFrames);
        end
        idxSel = min(max(1, idxSel), numel(times));
    end

    function idxSel = pickNearestHSIIndex(hsiTimes, tNow)
        idxSel = NaN;
        if isempty(hsiTimes) || isempty(tNow) || any(isnat(tNow))
            return;
        end

        hsiTimes = hsiTimes(:);
        [hsiTimesSorted, ord] = sort(hsiTimes);

        if tNow <= hsiTimesSorted(1)
            idxSel = ord(1);
            return;
        elseif tNow >= hsiTimesSorted(end)
            idxSel = ord(end);
            return;
        end

        diffs = abs(hsiTimesSorted - tNow);
        minDiff = min(diffs);
        idxLocal = find(diffs == minDiff, 1, 'first');
        idxSel = ord(idxLocal);
    end

    function [scanId, label] = parseHsiScanId(pathStr, fallbackIdx)
        if nargin < 2 || isempty(fallbackIdx)
            fallbackIdx = 1;
        end
        scanId = NaN;
        label = '';
        [~, fname, ~] = fileparts(pathStr);

        patterns = { 'Scan_(\d+)', 'scan(\d+)', '_(\d{3,})$' };
        for ii = 1:numel(patterns)
            tok = regexp(fname, patterns{ii}, 'tokens', 'once', 'ignorecase');
            if ~isempty(tok)
                scanId = str2double(tok{1});
                if isnan(scanId)
                    scanId = fallbackIdx;
                end
                label = sprintf('Scan_%05d', scanId);
                return;
            end
        end

        scanId = fallbackIdx;
        label = sprintf('Scan %d', fallbackIdx);
    end

    function [sensorKey, modality] = sensorKeyForEvent(evt)
        sensorKey = '';
        modality = '';
        if ~isfield(evt,'sensor') || isempty(evt.sensor)
            return;
        end

        baseSensor = upper(evt.sensor);
        switch baseSensor
            case 'CERB'
                if isfield(evt,'modality') && ~isempty(evt.modality)
                    modality = upper(evt.modality);
                else
                    modality = 'LWIR';
                end
                sensorKey = ['CERB_' modality];
            case 'MX20'
                modality = 'SWIR';
                sensorKey = 'MX20';
            case 'FAST'
                if isfield(evt,'modality') && ~isempty(evt.modality)
                    modality = upper(evt.modality);
                else
                    modality = 'LWIR';
                end
                sensorKey = ['FAST_' modality];
        end
    end

    function [sensorBase, modality] = parseSensorKey(sensorKey)
        sensorBase = '';
        modality   = '';
        if startsWith(sensorKey, 'CERB_')
            sensorBase = 'CERB';
            modality   = extractAfter(sensorKey, 'CERB_');
        elseif strcmp(sensorKey, 'MX20')
            sensorBase = 'MX20';
            modality   = 'SWIR';
        elseif startsWith(sensorKey, 'FAST_')
            sensorBase = 'FAST';
            modality   = extractAfter(sensorKey, 'FAST_');
        end
    end

    function dispName = displayNameForSensor(sensorKey)
        dispName = sensorKey;
        if isKey(hsiControlMap, sensorKey)
            ctrl = hsiControlMap(sensorKey);
            if isfield(ctrl,'displayName') && ~isempty(ctrl.displayName)
                dispName = ctrl.displayName;
                return;
            end
        end
        if startsWith(sensorKey, 'CERB_')
            dispName = sprintf('CERBERUS %s', extractAfter(sensorKey, 'CERB_'));
        elseif startsWith(sensorKey, 'FAST_')
            dispName = sprintf('FAST %s', extractAfter(sensorKey, 'FAST_'));
        elseif strcmp(sensorKey, 'MX20')
            dispName = 'MX20 SW';
        end
    end

    function rebuildHsiGroups()
        S.hsiGroupsMap = containers.Map('KeyType','char','ValueType','any');
        S.currentHsiMap = containers.Map('KeyType','char','ValueType','any');

        if ~S.enableHSI || isempty(S.hsiEvents)
            disableAllScanControls();
            return;
        end

        tmpMap = containers.Map('KeyType','char','ValueType','any');
        for ii = 1:numel(S.hsiEvents)
            evt = S.hsiEvents(ii);
            [sensorKey, modality] = sensorKeyForEvent(evt);
            if isempty(sensorKey)
                continue;
            end
            evt.modality = modality;
            entry = getOr(tmpMap, sensorKey, struct('events', struct([])));
            entry.events = [entry.events; evt]; %#ok<AGROW>
            tmpMap(sensorKey) = entry;
        end

        tmpKeys = tmpMap.keys;
        for kk = 1:numel(tmpKeys)
            key = tmpKeys{kk};
            events = tmpMap(key).events;
            if isempty(events)
                continue;
            end
            times = [events.time]';
            if isempty(times)
                continue;
            end
            timesUnique = unique(times);
            groups = struct('time', {}, 'items', {}, 'defaultIdx', {});
            for tt = 1:numel(timesUnique)
                tVal = timesUnique(tt);
                idxMask = find(times == tVal);
                items = struct('path', {}, 'scanId', {}, 'label', {}, 'modality', {}, 'sensor', {});
                for jj = 1:numel(idxMask)
                    evt = events(idxMask(jj));
                    [scanId, label] = parseHsiScanId(evt.path, jj);
                    items(end+1) = struct('path', evt.path, 'scanId', scanId, ...
                        'label', label, 'modality', evt.modality, 'sensor', evt.sensor); %#ok<AGROW>
                end
                scanIds = [items.scanId];
                [~, ord] = sort(scanIds);
                items = items(ord);
                groups(end+1) = struct('time', tVal, 'items', items, 'defaultIdx', 1); %#ok<AGROW>
            end
            S.hsiGroupsMap(key) = struct('timesUnique', timesUnique(:), 'groups', groups, ...
                'displayName', displayNameForSensor(key));
        end

        disableAllScanControls();
    end

    function disableAllScanControls()
        if isempty(hsiControlMap)
            return;
        end
        keys = hsiControlMap.keys;
        for ii = 1:numel(keys)
            ctrl = hsiControlMap(keys{ii});
            if isfield(ctrl,'btnPrev') && isgraphics(ctrl.btnPrev)
                ctrl.btnPrev.Enable = 'off';
            end
            if isfield(ctrl,'btnNext') && isgraphics(ctrl.btnNext)
                ctrl.btnNext.Enable = 'off';
            end
            if isfield(ctrl,'label') && isgraphics(ctrl.label)
                ctrl.label.Text = sprintf('%s | Scan: -', displayNameForSensor(keys{ii}));
            end
        end
    end

    function setScanControlsState(sensorKey, groupIdx, itemIdx)
        if ~isKey(hsiControlMap, sensorKey)
            return;
        end
        ctrl = hsiControlMap(sensorKey);
        lblText = sprintf('%s | Scan: -', displayNameForSensor(sensorKey));
        prevState = 'off';
        nextState = 'off';

        data = getOr(S.hsiGroupsMap, sensorKey, []);
        if ~S.enableHSI || isempty(data) || ~isfield(data,'groups') || isempty(data.groups) || ...
                isnan(groupIdx)
            ctrl.btnPrev.Enable = prevState;
            ctrl.btnNext.Enable = nextState;
            ctrl.label.Text = lblText;
            return;
        end

        groups = data.groups;
        if groupIdx < 1 || groupIdx > numel(groups)
            ctrl.btnPrev.Enable = prevState;
            ctrl.btnNext.Enable = nextState;
            ctrl.label.Text = lblText;
            return;
        end

        group = groups(groupIdx);
        nItems = numel(group.items);
        if nItems < 1
            ctrl.btnPrev.Enable = prevState;
            ctrl.btnNext.Enable = nextState;
            ctrl.label.Text = lblText;
            return;
        end

        itemIdx = max(1, min(nItems, itemIdx));
        item = group.items(itemIdx);
        prevState = 'on';
        nextState = 'on';
        if nItems <= 1
            prevState = 'off';
            nextState = 'off';
        else
            if itemIdx <= 1
                prevState = 'off';
            end
            if itemIdx >= nItems
                nextState = 'off';
            end
        end

        lblText = sprintf('%s | %s | %s (%d/%d)', data.displayName, ...
            formatHsiTimestamp(group.time), item.label, itemIdx, nItems);

        ctrl.btnPrev.Enable = prevState;
        ctrl.btnNext.Enable = nextState;
        ctrl.label.Text = lblText;
    end

    function stepScan(sensorKey, delta)
        if nargin < 2
            delta = 0;
        end
        data = getOr(S.hsiGroupsMap, sensorKey, []);
        if isempty(data) || ~isfield(data,'groups') || isempty(data.groups)
            setScanControlsState(sensorKey, NaN, NaN);
            return;
        end

        state = getOr(S.currentHsiMap, sensorKey, struct('groupIdx', NaN, 'itemIdx', NaN));

        % Prefer the currently displayed group/item so the first click uses
        % the visible scan as its starting point instead of re-aligning to
        % the nearest timestamp and swallowing the delta.
        groupIdx = NaN;
        itemIdx = NaN;
        if isfield(state,'groupIdx') && isfield(state,'itemIdx')
            if ~isnan(state.groupIdx) && state.groupIdx >= 1 && state.groupIdx <= numel(data.groups)
                groupIdx = state.groupIdx;
                itemIdx = state.itemIdx;
            end
        end

        if isnan(groupIdx)
            tRef = timeForFrameSafe(S.frame);
            nearestIdx = pickNearestHSIIndex(data.timesUnique, tRef);
            if isnan(nearestIdx)
                nearestIdx = 1;
            end
            groupIdx = nearestIdx;
            itemIdx = data.groups(groupIdx).defaultIdx;
        end

        nItems = numel(data.groups(groupIdx).items);
        if isnan(itemIdx) || itemIdx < 1 || itemIdx > nItems
            itemIdx = data.groups(groupIdx).defaultIdx;
            nItems = numel(data.groups(groupIdx).items);
        end

        targetIdx = max(1, min(nItems, itemIdx + delta));
        updateHSIPane(sensorKey, groupIdx, targetIdx);
    end

    function updateHSIPane(sensorKey, groupIdx, itemIdx)
        data = getOr(S.hsiGroupsMap, sensorKey, []);
        if isempty(data) || ~isfield(data,'groups') || isempty(data.groups)
            setScanControlsState(sensorKey, NaN, NaN);
            return;
        end

        if groupIdx < 1 || groupIdx > numel(data.groups)
            setScanControlsState(sensorKey, NaN, NaN);
            return;
        end
        group = data.groups(groupIdx);
        if isempty(group.items)
            setScanControlsState(sensorKey, NaN, NaN);
            return;
        end

        if isnan(itemIdx) || itemIdx < 1
            itemIdx = group.defaultIdx;
        end
        itemIdx = max(1, min(numel(group.items), itemIdx));
        item = group.items(itemIdx);

        prev = getOr(S.currentHsiMap, sensorKey, struct());
        doLayout = ~isfield(prev,'path') || isempty(prev.path);
        if isfield(prev,'groupIdx') && isfield(prev,'itemIdx') && isfield(prev,'path') && ...
                prev.groupIdx == groupIdx && prev.itemIdx == itemIdx && strcmp(prev.path, item.path)
            setScanControlsState(sensorKey, groupIdx, itemIdx);
            return;
        end

        [sensorBase, modality] = parseSensorKey(sensorKey);
        switch sensorBase
            case 'CERB'
                loadCerbFromPath(modality, item.path, doLayout);
            case 'MX20'
                loadMX20FromHdr(item.path, doLayout);
            case 'FAST'
                ensureFastAxis(modality);
                loadFastFromHdr(item.path, modality, doLayout);
        end

        effTime = readHsiUnixTime(item.path);
        if isempty(effTime) || isnat(effTime)
            effTime = group.time;
        end

        newState = struct('sensor', sensorBase, 'modality', modality, 'groupIdx', groupIdx, ...
            'itemIdx', itemIdx, 'path', item.path, 'time', group.time, 'effectiveTime', effTime);
        S.currentHsiMap(sensorKey) = newState;
        S.currentHsi = newState;
        setScanControlsState(sensorKey, groupIdx, itemIdx);
    end

    function out = timeForFrameSafe(idx)
        out = NaT;
        if hasFridgeTimes() && idx <= numel(S.timelineTimes)
            out = timeForFrame(idx);
        elseif ~isnat(S.tNow)
            out = S.tNow;
        end
    end

    function s = formatHsiTimestamp(t)
        if isdatetime(t) && ~isnat(t)
            try
                tfmt = t;
                tfmt.Format = 'HH:mm:ss.SSS';
                s = char(tfmt);
                return;
            catch
            end
        end
        s = '(time unavailable)';
    end

    function syncHsiToTime(tTarget)
        if ~S.enableHSI || isempty(S.hsiGroupsMap)
            disableAllScanControls();
            return;
        end

        if nargin < 1 || isempty(tTarget)
            allTimes = datetime.empty(0,1);
            keys = S.hsiGroupsMap.keys;
            for ii = 1:numel(keys)
                data = S.hsiGroupsMap(keys{ii});
                if isfield(data,'timesUnique')
                    allTimes = [allTimes; data.timesUnique]; %#ok<AGROW>
                end
            end
            if isempty(allTimes)
                disableAllScanControls();
                return;
            end
            tTarget = min(allTimes);
        end

        keys = S.hsiGroupsMap.keys;
        for ii = 1:numel(keys)
            sensorKey = keys{ii};
            data = S.hsiGroupsMap(sensorKey);
            if isempty(data) || ~isfield(data,'timesUnique') || isempty(data.timesUnique)
                setScanControlsState(sensorKey, NaN, NaN);
                continue;
            end
            groupIdx = pickNearestHSIIndex(data.timesUnique, tTarget);
            if isnan(groupIdx)
                setScanControlsState(sensorKey, NaN, NaN);
                continue;
            end
            prev = getOr(S.currentHsiMap, sensorKey, struct('groupIdx', NaN, 'itemIdx', NaN));
            if ~isfield(prev,'groupIdx') || isnan(prev.groupIdx) || prev.groupIdx ~= groupIdx
                itemIdx = data.groups(groupIdx).defaultIdx;
            else
                itemIdx = prev.itemIdx;
                if isnan(itemIdx) || itemIdx < 1 || itemIdx > numel(data.groups(groupIdx).items)
                    itemIdx = data.groups(groupIdx).defaultIdx;
                end
            end
            updateHSIPane(sensorKey, groupIdx, itemIdx);
        end
    end

    function tHdr = readHsiUnixTime(hsicOrHdrPath)
        tHdr = [];
        if nargin < 1 || isempty(hsicOrHdrPath)
            return;
        end

        [p,n,ext] = fileparts(hsicOrHdrPath);
        if strcmpi(ext,'.hsic')
            hdrPath = fullfile(p, [n '.hdr']);
        else
            hdrPath = hsicOrHdrPath;
        end

        if ~isfile(hdrPath)
            return;
        end

        try
            txt = fileread(hdrPath);
            tok = regexp(txt, 'unixtime\s*=\s*([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once', 'ignorecase');
            if isempty(tok)
                tok = regexp(txt, 'unix time\s*=\s*([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once', 'ignorecase');
            end
            if ~isempty(tok)
                uVal = str2double(tok{1});
                if ~isnan(uVal)
                    tHdr = datetime(uVal, 'ConvertFrom', 'posixtime');
                end
            end
        catch
            % Ignore header read failures and fall back to coarse time
        end
    end

    function tEff = effectiveHsiTime(evt)
        tEff = evt.time;
        if ~isfield(evt,'path') || isempty(evt.path)
            return;
        end

        key = evt.path;
        if isKey(S.hsiPreciseCache, key)
            cached = S.hsiPreciseCache(key);
            if isdatetime(cached)
                tEff = cached;
                return;
            end
        end

        tHdr = readHsiUnixTime(evt.path);
        if ~isempty(tHdr)
            tEff = tHdr;
        end

        S.hsiPreciseCache(key) = tEff;
    end

    function updateHsiJumpAvailability()
        if ~S.enableHSI || isempty(S.hsiEvents)
            btnJumpHsi.Enable = 'off';
        else
            btnJumpHsi.Enable = 'on';
        end
    end

    function jumpToHsi()
        if isempty(S.hsiEvents)
            uialert(f, 'No HSI events are loaded to match.', 'No HSI Events');
            return;
        end

        tTarget = S.currentHsi.effectiveTime;
        if isempty(tTarget) || isnat(tTarget)
            hsiTimes = arrayfun(@(e) effectiveHsiTime(e), S.hsiEvents);
            hsiTimes = hsiTimes(~isnat(hsiTimes));
            if isempty(hsiTimes)
                uialert(f, 'HSI events are missing valid timestamps.', 'No HSI Time');
                return;
            end
            tTarget = min(hsiTimes);
        end

        jumpToTime(tTarget);
    end

    function saveSnapshot()
        defaultName = fullfile(S.dir, sprintf('montage_frame_%d.png', S.frame));
        [file, path] = uiputfile({'*.png'}, 'Save Snapshot', defaultName);
        if isequal(file,0)
            return;
        end

        opts = struct('time', timeForFrame(S.frame), 'targetSize', [1080 1920], ...
                      'outputPath', fullfile(path, file));
        try
            RMBV_Export.exportSnapshot(S, opts);
            uialert(f,'Snapshot saved.','Done');
        catch ME
            uialert(f, sprintf('Snapshot failed:\n\n%s', ME.message), 'Snapshot Error');
        end
    end

    function exportVideo()
        if S.nFrames < 1
            uialert(f, 'Load FRIDGE data before exporting.', 'No Data');
            return;
        end

        tsTag = datestr(now, 'yyyymmdd_HHMMSS');
        defaultPath = fullfile(S.dir, sprintf('montage_export_%s.mp4', tsTag));
        [file, path] = uiputfile({'*.mp4','MPEG-4 Video';'*.avi','Motion JPEG AVI'}, ...
                                 'Export Montage Video', defaultPath);
        if isequal(file,0)
            return;
        end

        timeChoices = {'Master FRIDGE timestamps (recommended)', 'Fixed FPS time base'};
        [timeChoiceIdx, okChoice] = listdlg('ListString', timeChoices, 'SelectionMode','single', ...
            'InitialValue', 1, 'PromptString','Choose export time base');
        if isempty(timeChoiceIdx) || ~okChoice
            return;
        end
        if timeChoiceIdx == 2
            timeBase = 'fixed';
        else
            timeBase = 'master';
        end

        if strcmp(timeBase, 'fixed')
            defaults = { '1920x1080', '15', '' };
            prompts  = { 'Resolution (HxW)', 'Frames per second (time spacing)', ...
                         'Time step (seconds, optional overrides FPS)' };
            dlgAns = inputdlg(prompts, 'Video Export Options (Fixed Time Base)', [1 60], defaults);
            if isempty(dlgAns)
                return;
            end
            resStr   = dlgAns{1};
            fpsStr   = dlgAns{2};
            timeStepStr = dlgAns{3};
            stepVal = 1;
        else
            defaults = { '1920x1080', '15', '1' };
            prompts  = { 'Resolution (HxW)', 'Frames per second (playback rate)', ...
                         'Every N frames from master modality' };
            dlgAns = inputdlg(prompts, 'Video Export Options (Master Time Base)', [1 60], defaults);
            if isempty(dlgAns)
                return;
            end
            resStr   = dlgAns{1};
            fpsStr   = dlgAns{2};
            stepStr  = dlgAns{3};
            timeStepStr = '';
            stepVal    = str2double(stepStr); if isnan(stepVal) || stepVal < 1, stepVal = 1; end
        end

        targetSize = parseResolution(resStr, [1080 1920]);
        fpsVal     = str2double(fpsStr); if isnan(fpsVal) || fpsVal <= 0, fpsVal = 15; end
        timeStepVal= str2double(timeStepStr); if isnan(timeStepVal) || timeStepVal <= 0, timeStepVal = []; end

        previewOpts = struct('frameStep', stepVal, 'timeStep', timeStepVal, ...
                             'fps', fpsVal, 'timeBase', timeBase);
        [nFrames, times, meta] = RMBV_Export.estimateFrameCount(S, previewOpts);
        if nFrames < 1
            uialert(f, 'No frames were selected for export. Adjust the step or time range.', 'Export Options');
            return;
        end

        estSeconds = nFrames / max(1, fpsVal);
        timeSummary = '';
        if isfield(meta,'masterModality') && ~isempty(meta.masterModality)
            timeSummary = sprintf('Master modality: %s\n', meta.masterModality);
        end
        if isfield(meta,'startTime') && isdatetime(meta.startTime) && ~isnat(meta.startTime)
            startStr = formatDatetimeSafe(meta.startTime);
            endStr   = formatDatetimeSafe(meta.endTime);
            timeSummary = sprintf('%sStart: %s\nEnd: %s\n', timeSummary, startStr, endStr);
        end
        if isfield(meta,'cadenceSeconds') && ~isnan(meta.cadenceSeconds)
            timeSummary = sprintf('%sCadence: %0.3f sec/frame (~%0.1f fps)\n', ...
                timeSummary, meta.cadenceSeconds, 1/max(meta.cadenceSeconds, eps));
        end
        msg = sprintf(['This export will write %d frames (~%0.1f seconds at %0.1f fps).\n' ...
            '%sContinue?'], nFrames, estSeconds, fpsVal, timeSummary);
        choice = questdlg(msg, 'Confirm Export Size', 'Yes', 'Cancel', 'Yes');
        if ~strcmp(choice, 'Yes')
            return;
        end

        opts = struct('targetSize', targetSize, 'fps', fpsVal, 'frameStep', stepVal, ...
                      'timeStep', timeStepVal, 'outputPath', fullfile(path,file), ...
                      'parentFigure', f, 'times', times, 'timeBase', timeBase);
        try
            RMBV_Export.exportVideo(S, opts);
            uialert(f, 'Export complete.', 'Done');
        catch ME
            uialert(f, sprintf('Export failed:\n\n%s', ME.message), 'Export Error');
        end
    end

    %======================== CERBERUS / MX20 ==============================
    function loadCerbFromPath(whichMod, fullpath, doLayout)
        if nargin < 3 || isempty(doLayout)
            doLayout = true;
        end
        if ~isfile(fullpath)
            uialert(f, sprintf('CERBERUS file not found:\n%s', fullpath), ...
                'CERBERUS Error');
            return;
        end
        [p,n,ext] = fileparts(fullpath);
        % Accept either the calibrated cube (.hsic) or a header that sits
        % next to it and map headers to their cube automatically.
        if strcmpi(ext,'.hdr')
            fullpath = fullfile(p, [n '.hsic']);
            ext = '.hsic';
        end
        if ~strcmpi(ext,'.hsic')
            uialert(f, ...
                'This tool only accepts calibrated CERBERUS cubes (.hsic).', ...
                'CERBERUS File Type');
            return;
        end
        if ~isfile(fullpath)
            uialert(f, sprintf('Expected CERBERUS cube not found:\n%s', fullpath), ...
                'CERBERUS Error');
            return;
        end
        try
            ctx = loadCerberusContext(fullpath);  % your existing helper
            ctx = rot90(ctx, -1);
        catch ME
            uialert(f, sprintf('Failed to read CERBERUS cube:\n\n%s', ME.message), ...
                'CERBERUS Error');
            return;
        end

        switch upper(whichMod)
            case 'LWIR'
                hImg = imshow(ctx, [], 'Parent', cerbAxLWIR);
                hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, 'CERB LWIR', ctx);
                hImg.PickableParts = 'all';
                hImg.HitTest       = 'on';
                [~,fnOnly,~] = fileparts(fullpath);
                title(cerbAxLWIR, ['CERB LWIR — ', fnOnly], 'Interpreter','none');
                S.cerb.LWIR = struct('path',fullpath,'ctx',ctx);
            case 'VNIR'
                hImg = imshow(ctx, [], 'Parent', cerbAxVNIR);
                hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, 'CERB VNIR', ctx);
                hImg.PickableParts = 'all';
                hImg.HitTest       = 'on';
                [~,fnOnly,~] = fileparts(fullpath);
                title(cerbAxVNIR, ['CERB VNIR — ', fnOnly], 'Interpreter','none');
                S.cerb.VNIR = struct('path',fullpath,'ctx',ctx);
        end

        if doLayout
            refreshMontageLayout();
        end
    end

    function loadMX20FromHdr(hdrOrHsicPath, doLayout)
        if nargin < 2 || isempty(doLayout)
            doLayout = true;
        end
        if ~isfile(hdrOrHsicPath)
            uialert(f, sprintf('MX20 header/file not found:\n%s', hdrOrHsicPath), ...
                'MX20 Error');
            return;
        end
        [p, n, ext] = fileparts(hdrOrHsicPath);
        if strcmpi(ext, '.hdr')
            hsicPath = fullfile(p, [n '.hsic']);
            if ~isfile(hsicPath)
                uialert(f, sprintf(['Expected MX20 .hsic next to header, but ' ...
                                    'could not find:\n%s'], hsicPath), ...
                        'MX20 Error');
                return;
            end
        else
            hsicPath = hdrOrHsicPath;
        end
        try
            ctx = loadCerberusContext(hsicPath);
            ctx = rot90(ctx, -1);
        catch ME
            uialert(f, sprintf('Failed to read MX20 cube:\n\n%s', ME.message), ...
                'MX20 Error');
            return;
        end
        hImg = imshow(ctx, [], 'Parent', mxAx);
        hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, 'MX20 SW', ctx);
        hImg.PickableParts = 'all';
        hImg.HitTest       = 'on';
        [~,fnOnly,~] = fileparts(hsicPath);
        title(mxAx, ['MX20 SW — ', fnOnly], 'Interpreter','none');
        S.mx20 = struct('hdr',hdrOrHsicPath,'ctx',ctx);

        if doLayout
            refreshMontageLayout();
        end
    end

    function loadFastFromHdr(hdrOrHsicPath, modality, doLayout)
        if nargin < 2 || isempty(modality)
            modality = 'LWIR';
        end
        if nargin < 3 || isempty(doLayout)
            doLayout = true;
        end
        key = upper(modality);
        if ~isfile(hdrOrHsicPath)
            uialert(f, sprintf('FAST header/file not found:\n%s', hdrOrHsicPath), ...
                'FAST Error');
            return;
        end
        [p, n, ext] = fileparts(hdrOrHsicPath);
        if strcmpi(ext, '.hdr')
            hsicPath = fullfile(p, [n '.hsic']);
            if ~isfile(hsicPath)
                uialert(f, sprintf(['Expected FAST .hsic next to header, but could ' ...
                                    'not find:\n%s'], hsicPath), ...
                        'FAST Error');
                return;
            end
        else
            hsicPath = hdrOrHsicPath;
        end
        try
            ctx = loadCerberusContext(hsicPath);
            ctx = rot90(ctx, -1);
        catch ME
            uialert(f, sprintf('Failed to read FAST cube:\n\n%s', ME.message), ...
                'FAST Error');
            return;
        end

        ax = ensureFastAxis(key);
        hImg = imshow(ctx, [], 'Parent', ax);
        hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, ['FAST ' key], ctx);
        hImg.PickableParts = 'all';
        hImg.HitTest       = 'on';
        [~,fnOnly,~] = fileparts(hsicPath);
        title(ax, sprintf('FAST %s — %s', key, fnOnly), 'Interpreter','none');
        S.fast.(key) = struct('hdr',hdrOrHsicPath,'ctx',ctx);

        if doLayout
            refreshMontageLayout();
        end
    end

    %======================== DRAWING / IO (FRIDGE) ========================
    function drawAll(tCurrent)
        if nargin < 1 || isempty(tCurrent) || isnat(tCurrent)
            tCurrent = S.tNow;
        end
        for i = 1:numel(modalities)
            m  = keyify(modalities{i});
            ax = axMap(m);
            frameLbl = frameLabelMap(m);
            fileLbl  = fileLabelMap(m);
            prevIdx = getOr(S.lastFrameByMod, m, NaN);
            prevStatus = getOr(S.lastStatusByMod, m, '');
            newStatus = '';
            newIdx = NaN;

            if ~hasKey(S.exists, m) || ~hasKey(S.hdrs, m) || ~hasKey(S.files, m)
                newStatus = 'missing-metadata';
                newIdx = NaN;
                if isequal(prevIdx, newIdx) && strcmp(prevStatus, newStatus)
                    continue;
                end
                cla(ax);
                axis(ax,'off');
                title(ax, '');
                frameLbl.Text = sprintf('%s — Missing metadata', m);
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            if ~getOr(S.exists, m, false)
                newStatus = 'missing-file';
                newIdx = NaN;
                if isequal(prevIdx, newIdx) && strcmp(prevStatus, newStatus)
                    continue;
                end
                cla(ax);
                axis(ax,'off');
                title(ax, '');
                frameLbl.Text = sprintf('%s — Missing file', m);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            hdrCandidate = getOr(S.hdrs, m, []);
            if isempty(hdrCandidate)
                newStatus = 'missing-header';
                newIdx = NaN;
                if isequal(prevIdx, newIdx) && strcmp(prevStatus, newStatus)
                    continue;
                end
                cla(ax);
                axis(ax,'off');
                title(ax, '');
                frameLbl.Text = sprintf('%s — Missing header', m);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            [fEff, status] = effectiveFrame(m, tCurrent);
            newStatus = status;
            newIdx = fEff;
            if isequal(prevIdx, newIdx) && strcmp(prevStatus, newStatus)
                continue;
            end

            cla(ax);
            title(ax, '');
            frameLbl.Text = 'Frame: -';
            fileLbl.Text  = '';

            maxF = getOr(S.maxFrames, m, NaN);
            if isnan(fEff)
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — No frame (max %g)', m, maxF);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            try
                img = fridge_read_frame(m, fEff, S.hdrs, S.files);
            catch ME
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — %s', m, ME.message);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            % VIS: show in color if 3-channel; otherwise grayscale
            if strcmp(m,'VIS-COLOR') && ndims(img)==3 && size(img,3)==3
                imgDisp = img;
                if ~isfloat(imgDisp)
                    maxVal = double(max(imgDisp(:)));
                    if maxVal > 0
                        imgDisp = double(imgDisp) / maxVal;
                    else
                        imgDisp = zeros(size(imgDisp),'double');
                    end
                else
                    mxv = max(imgDisp(:));
                    if mxv > 1
                        imgDisp = imgDisp ./ mxv;
                    end
                end
                hImg = imshow(imgDisp, 'Parent', ax);
            else
                hImg = imshow(img, [], 'Parent', ax);
            end


            hImg.ButtonDownFcn = @(src,evt) onPixelClick(src, evt, m, img);
            hImg.PickableParts = 'all';
            hImg.HitTest       = 'on';

            switch status
                case 'ok'
                    note = '';
                case 'held'
                    note = sprintf(' (held @%g of %g)', fEff, maxF);
                case 'looped'
                    note = sprintf(' (looped @%g of %g)', fEff, maxF);
                otherwise
                    note = '';
            end
            frameLbl.Text = sprintf('%s — Frame %d/%d%s', m, fEff, maxF, note);
            [~,fn,ext] = fileparts(getOr(S.files, m, ''));
            fileLbl.Text = [fn ext];
            S.lastFrameByMod(m) = newIdx;
            S.lastStatusByMod(m) = newStatus;
        end
    end

    function [fEff, status] = effectiveFrame(modality, targetTime)
        modality = keyify(modality);
        if nargin < 2
            targetTime = [];
        end

        maxF = getOr(S.maxFrames, modality, NaN);
        if isnan(maxF) || maxF < 1
            fEff   = NaN;
            status = 'missing';
            return;
        end

        tVec = getOr(S.fridgeTimesMap, modality, datetime.empty(0,1));
        if isempty(tVec) || ~isdatetime(tVec)
            if ~isnat(S.fridgeStartTime) && ~isnat(S.fridgeEndTime)
                durSec = seconds(S.fridgeEndTime - S.fridgeStartTime);
                if maxF <= 1
                    tVec = S.fridgeStartTime;
                else
                    offsets = linspace(0, durSec, maxF);
                    tVec = S.fridgeStartTime + seconds(offsets(:));
                end
            end
        end

        if isempty(tVec) || ~isdatetime(tVec)
            fEff   = NaN;
            status = 'missing';
            return;
        end

        tVec = tVec(~isnat(tVec));
        if isempty(tVec)
            fEff   = NaN;
            status = 'missing';
            return;
        end

        if isempty(targetTime) || ~isdatetime(targetTime)
            targetTime = S.tNow;
        end
        [fEff, status] = pickNearestIndex(tVec, targetTime, maxF);
    end

    %======================== PIXEL CLICKS =================================
    function onPixelClick(src, evt, modality, img) %#ok<INUSD>
        ax = ancestor(src, 'axes');
        cp = ax.CurrentPoint;
        x = round(cp(1,1));
        y = round(cp(1,2));
        sz = size(img);
        h  = sz(1);
        w  = sz(2);
        if x < 1 || x > w || y < 1 || y > h
            return;
        end
        v = img(y, x, :);
        lblPixel.Text = sprintf('Pixel: (%d, %d) — %s', x, y, modality);
        lblValue.Text = sprintf('Value: %s', viewer_format_pixel_value(v));
    end

    function onPixelClickCerb(src, evt, labelStr, img) %#ok<INUSD>
        ax = ancestor(src, 'axes');
        cp = ax.CurrentPoint;
        x = round(cp(1,1));
        y = round(cp(1,2));
        sz = size(img);
        h  = sz(1);
        w  = sz(2);
        if x < 1 || x > w || y < 1 || y > h
            return;
        end
        v = img(y, x, :);
        lblPixel.Text = sprintf('Pixel: (%d, %d) — %s', x, y, labelStr);
        lblValue.Text = sprintf('Value: %s', viewer_format_pixel_value(v));
    end

    %======================== TIME DISPLAY =================================
    function updateTimeDisplay()
        if isempty(S.tNow) || isnat(S.tNow)
            lblTime.Text = 'Time: (unavailable)';
            return;
        end
        lblTime.Text = sprintf('Time: %s', formatTimeCursor(S.tNow));
    end

    function s = formatTimeCursor(t)
        if isdatetime(t) && ~isnat(t)
            try
                tfmt = t;
                tfmt.Format = 'HH:mm:ss.SSS';
                s = char(tfmt);
                return;
            catch
            end
        end
        s = '(unavailable)';
    end

    function s = formatClockRange(t1, t2)
        s = sprintf('%s-%s', formatTimeCursor(t1), formatTimeCursor(t2));
    end

    function s = formatDatetimeSafe(t)
        if isdatetime(t) && ~isnat(t)
            try
                y = year(t);
                if any(~isfinite(y)) || any(y < 0 | y > 9999)
                    s = '(out of range)';
                    return;
                end
                tfmt = t;
                tfmt.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
                s = char(tfmt);
            catch
                try
                    s = datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF');
                catch
                    s = '(unavailable)';
                end
            end
        elseif isnumeric(t)
            try
                s = datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF');
            catch
                s = '(unavailable)';
            end
        else
            s = '(unavailable)';
        end
    end

    function updateReturnButtonState()
        if ~isempty(S.timelineFig) && isvalid(S.timelineFig)
            btnReturn.Enable = 'on';
        else
            btnReturn.Enable = 'off';
        end
    end

    function returnToTimeline()
        if isempty(S.timelineFig) || ~isvalid(S.timelineFig)
            btnReturn.Enable = 'off';
            return;
        end
        try
            S.timelineFig.Visible = 'on';
            if isprop(S.timelineFig,'WindowState')
                S.timelineFig.WindowState = 'normal';
            end
            figure(S.timelineFig);
        catch
            btnReturn.Enable = 'off';
        end
    end

    function onCloseRequest(figHandle)
        if isprop(figHandle, 'WindowState')
            try
                lastWindowState = figHandle.WindowState; %#ok<NASGU>
            catch
                % Ignore failures to capture window state.
            end
        end
        delete(figHandle);
    end

    %======================== RESET / INITIAL LOAD ========================
    function applyInitial(newInitial)
        % Allow timeline to reuse an existing viewer window instead of
        % recreating it; this resets and reloads using the latest selection.
        initial = newInitial;
        if nargin < 1 || isempty(initial)
            initial = struct();
        end

        enableHSI = true;
        if isfield(initial,'enableHSI') && ~isempty(initial.enableHSI)
            enableHSI = logical(initial.enableHSI);
        end

        targetStartTime = [];
        if isfield(initial,'initialTime')
            targetStartTime = initial.initialTime;
        end

        resetUI();

        if isfield(initial,'timelineFig') && ~isempty(initial.timelineFig) && isvalid(initial.timelineFig)
            S.timelineFig = initial.timelineFig;
        else
            S.timelineFig = [];
        end

        if enableHSI && isfield(initial,'hsiEvents')
            S.hsiEvents = initial.hsiEvents;
        else
            S.hsiEvents = struct('sensor', {}, 'time', {}, 'path', {}, 'modality', {});
        end
        S.enableHSI = enableHSI;
        if isfield(initial,'fridgeStartTime') && isfield(initial,'fridgeEndTime')
            S.fridgeStartTime = initial.fridgeStartTime;
            S.fridgeEndTime   = initial.fridgeEndTime;
        else
            S.fridgeStartTime = NaT;
            S.fridgeEndTime   = NaT;
        end
        [S.tStart, S.tEnd] = resolveTimeDomain(initial);
        if ~isnat(S.tStart)
            S.tNow = S.tStart;
        else
            S.tNow = NaT;
        end
        S.sliderStartTime = S.tStart;
        S.sliderEndTime   = S.tEnd;
        rebuildHsiGroups();
        updateHsiJumpAvailability();
        updateReturnButtonState();

        if isfield(initial,'rawFile') && ~isempty(initial.rawFile) && isfile(initial.rawFile)
            loadFromRawFile(initial.rawFile);
        end
        if enableHSI && isfield(initial,'cerbLWIR') && ~isempty(initial.cerbLWIR) && isfile(initial.cerbLWIR)
            loadCerbFromPath('LWIR', initial.cerbLWIR);
        end
        if enableHSI && isfield(initial,'cerbVNIR') && ~isempty(initial.cerbVNIR) && isfile(initial.cerbVNIR)
            loadCerbFromPath('VNIR', initial.cerbVNIR);
        end
        if enableHSI && isfield(initial,'mx20Hdr') && ~isempty(initial.mx20Hdr) && isfile(initial.mx20Hdr)
            loadMX20FromHdr(initial.mx20Hdr);
        end
        if enableHSI && isfield(initial,'fast') && ~isempty(initial.fast)
            mods = fieldnames(initial.fast);
            for ii = 1:numel(mods)
                hdrPath = initial.fast.(mods{ii});
                if isfile(hdrPath)
                    loadFastFromHdr(hdrPath, mods{ii});
                end
            end
        end

        recomputeSliderDataWindow();
        configureTimeSlider();
        if ~isnat(S.tNow)
            S.tNow = clampTime(S.tNow);
            updateAllPanesAtTime(S.tNow, true);
        elseif enableHSI && ~isempty(targetStartTime)
            syncHsiToTime(targetStartTime);
            updateTimeDisplay();
        else
            updateTimeDisplay();
            syncHsiToTime(timeForFrame(S.frame));
        end
    end

    function resetUI()
        S.files        = containers.Map(modalities, repmat({''},1,numel(modalities)));
        S.hdrs         = containers.Map(modalities, repmat({[]},1,numel(modalities)));
        S.exists       = containers.Map(modalities, num2cell(false(1,numel(modalities))));
        S.maxFrames    = containers.Map(modalities, num2cell(nan(1,numel(modalities))));
        S.frame        = 1;
        S.nFrames      = 0;
        S.frameCount   = 0;
        S.chosen       = '';
        S.fridgeTimes  = [];
        S.fridgeTimesMap = containers.Map(modalities, repmat({datetime.empty(0,1)},1,numel(modalities)));
        S.timelineTimes  = datetime.empty(0,1);
        S.sliderMode   = 'frame';
        S.tStart       = NaT;
        S.tEnd         = NaT;
        S.tNow         = NaT;
        S.sliderStartTime = NaT;
        S.sliderEndTime   = NaT;
        S.dataStartTime   = NaT;
        S.dataEndTime     = NaT;
        S.sliderRangeSec = NaN;
        S.sliderStepSec  = NaN;
        S.lastFrameByMod = containers.Map(modalities, num2cell(nan(1,numel(modalities))));
        S.lastStatusByMod = containers.Map(modalities, repmat({''},1,numel(modalities)));
        S.fridgeStartTime = NaT;
        S.fridgeEndTime   = NaT;
        S.hsiEvents    = struct('sensor', {}, 'time', {}, 'path', {}, 'modality', {});
        S.hsiGroupsMap = containers.Map('KeyType','char','ValueType','any');
        S.currentHsi   = struct('sensor','', 'time', NaT, 'effectiveTime', NaT);
        S.hsiPreciseCache = containers.Map('KeyType','char','ValueType','any');
        S.enableHSI = enableHSI;
        S.currentHsiMap = containers.Map('KeyType','char','ValueType','any');

        lblStatus.Text = 'Status: (no capture loaded)';
        lblFrames.Text = 'Frames: -';
        lblMem.Text    = '';
        lblPixel.Text  = 'Pixel: -';
        lblValue.Text  = 'Value: -';
        lblTime.Text   = 'Time: -';
        lblSelectionRange.Text = 'Selection: -';
        lblDataSpan.Text = 'Data span: -';

        disableAllScanControls();

        for i = 1:numel(modalities)
            modName = keyify(modalities{i});
            ax = getOr(axMap, modName, []);
            if ~isempty(ax)
                cla(ax);
                axis(ax,'off');
                title(ax, '');
            end

            % Use guarded map lookups in case upstream provided unexpected
            % modality keys; avoid containers.Map multi-level indexing errors.
            frameLbl = getOr(frameLabelMap, modName, []);
            if ~isempty(frameLbl)
                frameLbl.Text = 'Frame: -';
            end

            fileLbl = getOr(fileLabelMap, modName, []);
            if ~isempty(fileLbl)
                fileLbl.Text  = '';
            end
        end

        btnPrev.Enable     = 'off';
        btnNext.Enable     = 'off';
        btnSnapshot.Enable = 'off';
        btnExport.Enable   = 'off';
        btnJumpHsi.Enable  = 'off';

        frameSlider.Enable = 'off';
        frameSlider.Limits = [1 2];
        frameSlider.Value  = 1;

        cla(cerbAxLWIR); title(cerbAxLWIR,'CERB LWIR');
        cla(cerbAxVNIR); title(cerbAxVNIR,'CERB VNIR');
        cla(mxAx);       title(mxAx,'MX20 SW');
        if ~isempty(fastAxes)
            keys = fastAxes.keys;
            for ii = 1:numel(keys)
                ax = fastAxes(keys{ii});
                if isgraphics(ax)
                    cla(ax);
                    title(ax, sprintf('FAST %s', keys{ii}));
                end
            end
        end

        S.cerb = struct('LWIR',[],'VNIR',[]);
        S.mx20 = struct('hdr',[],'ctx',[]);
        S.fast = struct();

        refreshMontageLayout();
        updateReturnButtonState();
    end

    % Initial auto-load from timeline
    applyInitial(initial);
    setappdata(f, 'RMBV_Update', @applyInitial);

end
