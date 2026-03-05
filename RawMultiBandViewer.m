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

    % Header row with compact default controls + advanced drawer.
    headerWrap = uigridlayout(page,[3,1]);
    headerWrap.Layout.Row    = 1;
    headerWrap.Layout.Column = [1 3];
    headerWrap.RowHeight     = {'fit','fit','fit'};
    headerWrap.ColumnWidth   = {'1x'};
    headerWrap.Padding       = [8 8 8 8];

    headerRow = uigridlayout(headerWrap,[1,3]);
    headerRow.ColumnWidth   = {'1x','fit','fit'};

    header = makeLabel(headerRow, ...
        'Text','Multiband FRIDGE + HSI viewer (driven by timeline selection).', ...
        'FontWeight','bold','HorizontalAlignment','left');
    header.Layout.Row    = 1;
    header.Layout.Column = 1;

    btnAdvanced = uibutton(headerRow,'Text','[Advanced]', ...
        'ButtonPushedFcn',@(~,~)toggleAdvancedControls());
    btnAdvanced.Layout.Row = 1;
    btnAdvanced.Layout.Column = 2;

    btnReturn = uibutton(headerRow, 'Text','Return to Timeline', ...
        'Enable','off', ...
        'ButtonPushedFcn',@(~,~)returnToTimeline());
    btnReturn.Layout.Row    = 1;
    btnReturn.Layout.Column = 3;
    btnReturn.FontWeight    = 'bold';

    coreRow = uigridlayout(headerWrap,[1,6]);
    coreRow.ColumnWidth = {'fit','fit','fit','fit','1x','fit'};
    btnPrevOverlap = uibutton(coreRow,'Text','Prev Overlap','Enable','off', ...
        'ButtonPushedFcn',@(~,~)jumpOverlap(-1));
    btnNextOverlap = uibutton(coreRow,'Text','Next Overlap','Enable','off', ...
        'ButtonPushedFcn',@(~,~)jumpOverlap(1));
    btnResync = uibutton(coreRow,'Text','Re-sync','Enable','off','Visible','off', ...
        'ButtonPushedFcn',@(~,~)resyncPanels());
    lblDesync = makeLabel(coreRow,'Text','','FontWeight','bold');
    lblOverlapNote = makeLabel(coreRow,'Text','','HorizontalAlignment','left');
    makeLabel(coreRow,'Text','');

    advRow = uigridlayout(headerWrap,[1,10]);
    advRow.ColumnWidth = {'fit','fit','fit','fit','fit','fit','fit','fit','fit','1x'};
    advRow.Visible = 'off';
    makeLabel(advRow,'Text','Sync:');
    ddSyncMode = uidropdown(advRow,'Items',{'LOCKSTEP','FOLLOW_MASTER','FREE'},'Value','LOCKSTEP', ...
        'ValueChangedFcn',@(~,~)onSyncSettingsChanged());
    lblMaster = makeLabel(advRow,'Text','Master:');
    ddMasterSensor = uidropdown(advRow,'Items',{'LWIR'},'Value','LWIR','Enable','off', ...
        'ValueChangedFcn',@(~,~)onSyncSettingsChanged());
    makeLabel(advRow,'Text','Snap:');
    ddSnapMode = uidropdown(advRow,'Items',{'ANY','ALL','MASTER'},'Value','ANY', ...
        'ValueChangedFcn',@(~,~)onSyncSettingsChanged());
    makeLabel(advRow,'Text','Tolerance (ms):');
    fldTolerance = uieditfield(advRow,'numeric','Value',250,'Limits',[1 Inf], ...
        'RoundFractionalValues','on','ValueDisplayFormat','%.0f', ...
        'ValueChangedFcn',@(~,~)onSyncSettingsChanged());
    lblAdvancedHint = makeLabel(advRow,'Text','Advanced sync/snap settings', ...
        'HorizontalAlignment','left');

    % Axes names and grid positions
    modalities = {'LWIR','MWIR','SWIR','MONO','VIS-COLOR'};
    modalities = normalizeModalities(modalities);
    axMap          = containers.Map('KeyType','char','ValueType','any');
    frameLabelMap  = containers.Map('KeyType','char','ValueType','any');
    fileLabelMap   = containers.Map('KeyType','char','ValueType','any');
    panelMap       = containers.Map('KeyType','char','ValueType','any');
    hsiControlMap  = containers.Map('KeyType','char','ValueType','any');
    hsiFileLabelMap = containers.Map('KeyType','char','ValueType','any');

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

    function dispName = fridgeDisplayName(modName)
        if strcmp(modName,'VIS-COLOR')
            dispName = 'FRIDGE VIS';
        else
            dispName = ['FRIDGE ' modName];
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
    imgGrid.Padding       = [0 0 0 0];
    imgGrid.RowSpacing    = 0;
    imgGrid.ColumnSpacing = 0;

    % Off-screen parent used to hold inactive panels so the grid only sees
    % panes that are actually visible for the current selection.
    hiddenBin = uipanel(f, 'Visible','off', 'Position',[0 0 1 1]);

    %----------------------------------------------------------------------
    % FRIDGE panes: panel + inner grid (label row + axes row)
    %----------------------------------------------------------------------
    for i = 1:numel(modalities)
        pnl = uipanel(hiddenBin);

        pGrid = uigridlayout(pnl,[3,1]);
        pGrid.RowHeight   = {'fit','1x','fit'};
        pGrid.ColumnWidth = {'1x'};
        pGrid.Padding     = [0 0 0 0];
        pGrid.RowSpacing  = 0;
        pGrid.ColumnSpacing = 4;

        modName = keyify(modalities{i});
        dispName = fridgeDisplayName(modName);

        lblTop = makeLabel(pGrid, ...
            'Text', dispName, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');
        lblTop.Layout.Row    = 1;
        lblTop.Layout.Column = 1;

        ax = uiaxes(pGrid);
        ax.Layout.Row    = 2;
        ax.Layout.Column = 1;
        ax.DataAspectRatioMode = 'auto';
        ax.PlotBoxAspectRatioMode = 'auto';
        ax.XTick = [];
        ax.YTick = [];
        ax.XColor = 'none';
        ax.YColor = 'none';
        if isprop(ax, 'Toolbar')
            ax.Toolbar.Visible = 'off';
        end
        if isprop(ax, 'PositionConstraint')
            ax.PositionConstraint = 'innerposition';
        end
        if isprop(ax, 'LooseInset')
            ax.LooseInset = [0 0 0 0];
        end
        axis(ax,'off');
        title(ax, '');

        lblFilePane = makeLabel(pGrid, ...
            'Text','', ...
            'Interpreter','none', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top');
        lblFilePane.Layout.Row    = 3;
        lblFilePane.Layout.Column = 1;

        axMap(modName) = ax;
        frameLabelMap(modName) = lblTop;
        fileLabelMap(modName)  = lblFilePane;
        panelMap(modName)      = pnl;
    end

    % HSI panels (CERB LWIR/VNIR + MX20 + FAST)
    fastAxes = containers.Map('KeyType','char','ValueType','any');

    function [pnl, ax] = createHsiPanel(panelKey, displayName, sensorKey)
        pnl = uipanel(hiddenBin);
        pGrid = uigridlayout(pnl,[4,1]);
        pGrid.RowHeight   = {'fit','fit','1x','fit'};
        pGrid.ColumnWidth = {'1x'};
        pGrid.Padding     = [0 0 0 0];
        pGrid.RowSpacing  = 0;
        pGrid.ColumnSpacing = 4;

        lblTop = makeLabel(pGrid, ...
            'Text', displayName, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');
        lblTop.Layout.Row    = 1;
        lblTop.Layout.Column = 1;

        createHsiControlRow(pGrid, sensorKey, displayName, 2);

        ax = uiaxes(pGrid);
        ax.Layout.Row    = 3;
        ax.Layout.Column = 1;
        ax.DataAspectRatioMode = 'auto';
        ax.PlotBoxAspectRatioMode = 'auto';
        ax.XTick = [];
        ax.YTick = [];
        ax.XColor = 'none';
        ax.YColor = 'none';
        if isprop(ax, 'Toolbar')
            ax.Toolbar.Visible = 'off';
        end
        if isprop(ax, 'PositionConstraint')
            ax.PositionConstraint = 'innerposition';
        end
        if isprop(ax, 'LooseInset')
            ax.LooseInset = [0 0 0 0];
        end
        axis(ax,'off');
        title(ax, '');

        lblFilePane = makeLabel(pGrid, ...
            'Text','', ...
            'Interpreter','none', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top');
        lblFilePane.Layout.Row    = 4;
        lblFilePane.Layout.Column = 1;

        hsiFileLabelMap(panelKey) = lblFilePane;
        panelMap(panelKey) = pnl;
    end

    [~, cerbAxLWIR] = createHsiPanel('CERB_LWIR', 'CERBERUS LWIR', 'CERB_LWIR');
    [~, cerbAxVNIR] = createHsiPanel('CERB_VNIR', 'CERBERUS VNIR', 'CERB_VNIR');
    [~, mxAx] = createHsiPanel('MX20', 'MX20 SW', 'MX20');

    function modOut = normalizeFastModality(modIn)
        modOut = upper(string(modIn));
        if modOut == "FAST" || modOut == ""
            modOut = "LWIR";
        end
        if modOut ~= "LWIR" && modOut ~= "VNIR"
            modOut = "LWIR";
        end
        modOut = char(modOut);
    end

    function ax = ensureFastAxis(modality)
        key = normalizeFastModality(modality);
        if ~isKey(fastAxes, key) || isempty(fastAxes(key)) || ~isgraphics(fastAxes(key))
            [~, axNew] = createFastPanel(key);
            fastAxes(key) = axNew;
        end
        ax = fastAxes(key);
    end

    function [pnl, ax] = createFastPanel(modality)
        key = normalizeFastModality(modality);
        panelKey = sprintf('FAST_%s', key);
        displayName = sprintf('FAST %s', key);
        [pnl, ax] = createHsiPanel(panelKey, displayName, panelKey);
    end

    function createHsiControlRow(parent, sensorKey, displayName, layoutRow)
        if nargin < 4 || isempty(layoutRow)
            layoutRow = 1;
        end
        ctrlRow = uigridlayout(parent, [1 3]);
        ctrlRow.Layout.Row    = layoutRow;
        ctrlRow.Layout.Column = 1;
        ctrlRow.ColumnWidth   = {'fit','1x','fit'};
        ctrlRow.RowHeight     = {'fit'};
        ctrlRow.Padding       = [0 0 0 0];
        ctrlRow.RowSpacing    = 2;
        ctrlRow.ColumnSpacing = 4;

        btnPrevScan = uibutton(ctrlRow, 'Text','Prev Scan', 'Enable','off','Visible','off', ...
            'ButtonPushedFcn', @(~,~)stepScan(sensorKey, -1));
        btnPrevScan.Layout.Row    = 1;
        btnPrevScan.Layout.Column = 1;

        lblScan = makeLabel(ctrlRow, 'Text', sprintf('%s | Scan: -', displayName), ...
            'HorizontalAlignment','center');
        lblScan.Layout.Row    = 1;
        lblScan.Layout.Column = 2;

        btnNextScan = uibutton(ctrlRow, 'Text','Next Scan', 'Enable','off','Visible','off', ...
            'ButtonPushedFcn', @(~,~)stepScan(sensorKey, 1));
        btnNextScan.Layout.Row    = 1;
        btnNextScan.Layout.Column = 3;

        hsiControlMap(sensorKey) = struct('btnPrev', btnPrevScan, 'btnNext', btnNextScan, ...
            'label', lblScan, 'displayName', displayName);
    end

    % Create common FAST panels up front
    ensureFastAxis('LWIR');
    ensureFastAxis('VNIR');

    function setHsiFileLabel(panelKey, textValue)
        if isKey(hsiFileLabelMap, panelKey)
            lbl = hsiFileLabelMap(panelKey);
            if isgraphics(lbl)
                lbl.Text = textValue;
            end
        end
    end

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
    navCol = uigridlayout(ctrlWrapper,[4,1]);
    navCol.RowHeight   = {'fit','fit','fit','fit'};
    navCol.ColumnWidth = {'1x'};

    navTop = uigridlayout(navCol,[1,2]);
    navTop.ColumnWidth = {'fit','fit'};
    btnSnapshot = uibutton(navTop,'Text','Save Snapshot','Enable','off', ...
        'ButtonPushedFcn',@(~,~)saveSnapshot());
    btnExport = uibutton(navTop,'Text','Export Video...','Enable','off', ...
        'ButtonPushedFcn',@(~,~)exportVideo());

    instanceRow = uigridlayout(navCol,[1,3]);
    instanceRow.ColumnWidth = {'fit','fit','1x'};
    btnPrevInstance = uibutton(instanceRow,'Text','Prev Instance','Enable','off','Visible','on', ...
        'ButtonPushedFcn',@(~,~)stepInstance(-1));
    btnNextInstance = uibutton(instanceRow,'Text','Next Instance','Enable','off','Visible','on', ...
        'ButtonPushedFcn',@(~,~)stepInstance(+1));
    lblInstance = makeLabel(instanceRow,'Text','Instance: -','HorizontalAlignment','left','Visible','on');

    navBottom = uigridlayout(navCol,[1,3]);
    navBottom.ColumnWidth = {'fit','1x','fit'};
    lblSliderStart = makeLabel(navBottom, 'Text','Start: -', 'HorizontalAlignment','left');
    lblSliderStart.Layout.Row = 1;
    lblSliderStart.Layout.Column = 1;
    frameSlider = uislider(navBottom, ...
        'Limits',[1 2], ...
        'Value',1, ...
        'MajorTicks',[], ...
        'MinorTicks',[], ...
        'Enable','off', ...
        'ValueChangingFcn',@frameSliderChanging, ...
        'ValueChangedFcn',@frameSliderChanged);
    frameSlider.Layout.Row = 1;
    frameSlider.Layout.Column = 2;
    lblSliderEnd = makeLabel(navBottom, 'Text','End: -', 'HorizontalAlignment','right');
    lblSliderEnd.Layout.Row = 1;
    lblSliderEnd.Layout.Column = 3;

    fineRow = uigridlayout(navCol,[1,3]);
    fineRow.ColumnWidth = {'fit','1x','fit'};
    makeLabel(fineRow,'Text','FRIDGE frame:','HorizontalAlignment','right');
    fineSlider = uislider(fineRow, ...
        'Limits',[1 2], ...
        'Value',1, ...
        'MajorTicks',[], ...
        'MinorTicks',[], ...
        'Enable','off', ...
        'ValueChangingFcn',@fineSliderChanging, ...
        'ValueChangedFcn',@fineSliderChanged);
    lblFineFrame = makeLabel(fineRow,'Text','-','HorizontalAlignment','center');

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
    S.fridgeInstancesInRange = struct([]);
    S.activeFridgeInstanceIdx = NaN;
    S.isSwitchingFridge = false;
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
    S.instanceTimeline = struct('time', {}, 'type', {}, 'source', {}, 'ref', {});
    S.instanceTimes = datetime.empty(0,1);
    S.activeFineClip = struct('modality','', 'times', datetime.empty(0,1), 'startTime', NaT, 'endTime', NaT);
    S.perfEnabled = false;
    S.scrubPreviewFps = 20;
    S.lastScrubPreviewTic = tic;
    S.pendingScrubTime = NaT;
    S.renderJobId = 0;
    % Sync/snap model: playhead drives cross-sensor alignment.
    S.playheadTs = NaT;
    S.syncMode = 'LOCKSTEP';
    S.snapMode = 'ANY';
    S.masterSensor = '';
    S.selectedSensors = {};
    S.toleranceMs = 250;
    S.panelLocalTs = containers.Map('KeyType','char','ValueType','any');
    S.localScrubOverride = containers.Map('KeyType','char','ValueType','logical');
    S.showAdvanced = false;
    % Stale-while-revalidate frame cache: keep old image visible while loading new.
    S.frameCache = mv_lru_create(120);
    S.prefetchLoading = containers.Map('KeyType','char','ValueType','logical');
    S.abortedRequests = 0;
    S.perfSyncDebug = any(strcmp(getenv('RMBV_SYNC_DEBUG'), {'1','true','TRUE','on','ON'}));
    S.lastActivePanelsSig = '';

    targetStartTime = [];
    enableHSI = true;
    sliderInternalUpdate = false;  % prevent recursive slider callbacks
    fineSliderInternalUpdate = false;

    if isfield(initial,'perf') && logical(initial.perf)
        S.perfEnabled = true;
    end
    if ~S.perfEnabled
        perfEnv = getenv('RMBV_PERF');
        S.perfEnabled = any(strcmp(perfEnv, {'1','true','TRUE','on','ON'}));
    end

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

    function tf = shouldShowFridgePanel(modality)
        modality = keyify(modality);
        if ~getOr(S.exists, modality, false)
            tf = false;
            return;
        end
        if isempty(S.playheadTs) || ~isdatetime(S.playheadTs) || isnat(S.playheadTs)
            tf = true;
            return;
        end
        tVec = fridgeTimesForModality(modality, getOr(S.maxFrames, modality, NaN));
        if isempty(tVec)
            tf = false;
            return;
        end
        res = mv_resolve_panel_sample(tVec, S.playheadTs, struct('toleranceMs', S.toleranceMs, 'allowHold', false));
        tf = strcmp(res.status, 'OK');
    end

    function activePanels = getActivePanels()
        activePanels = {};
        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if shouldShowFridgePanel(m)
                activePanels{end+1} = m; %#ok<AGROW>
            end
        end

        if hasCerb('LWIR')
            activePanels{end+1} = 'CERB_LWIR'; %#ok<AGROW>
        end
        if hasCerb('VNIR')
            activePanels{end+1} = 'CERB_VNIR'; %#ok<AGROW>
        end
        if hasMx20()
            activePanels{end+1} = 'MX20'; %#ok<AGROW>
        end
        if isfield(S,'fast')
            fastMods = fieldnames(S.fast);
            for ii = 1:numel(fastMods)
                if hasFast(fastMods{ii})
                    activePanels{end+1} = sprintf('FAST_%s', normalizeFastModality(fastMods{ii})); %#ok<AGROW>
                end
            end
        end
    end

    function refreshMontageLayout()
        keysAll = panelMap.keys;
        for kk = 1:numel(keysAll)
            pnl = panelMap(keysAll{kk});
            pnl.Parent  = hiddenBin;
            pnl.Visible = 'off';
        end

        activePanels = getActivePanels();
        S.lastActivePanelsSig = strjoin(activePanels, '|');
        n = numel(activePanels);
        if n == 0
            imgGrid.RowHeight   = {'1x'};
            imgGrid.ColumnWidth = {'1x'};
            return;
        end

        maxCols = 4;
        cols = min(maxCols, max(1, ceil(sqrt(n))));
        rows = ceil(n / cols);
        imgGrid.RowHeight   = repmat({'1x'}, 1, rows);
        imgGrid.ColumnWidth = repmat({'1x'}, 1, cols);

        slot = 0;
        for idx = 1:n
            key = keyify(activePanels{idx});
            if ~isKey(panelMap, key)
                % Robustness: skip stale/unknown panel keys instead of erroring.
                continue;
            end
            slot = slot + 1;
            pnl = panelMap(key);
            pnl.Parent = imgGrid;
            pnl.Layout.Row    = ceil(slot / cols);
            pnl.Layout.Column = mod(slot-1, cols) + 1;
            pnl.Visible = 'on';
        end
    end


    function refreshMontageLayoutIfNeeded(force)
        if nargin < 1
            force = false;
        end
        activePanels = getActivePanels();
        sig = strjoin(activePanels, '|');
        if force || ~strcmp(S.lastActivePanelsSig, sig)
            refreshMontageLayout();
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
            'dataStartTime', S.dataStartTime, 'dataEndTime', S.dataEndTime, ...
            'fridgeInstancesInRange', S.fridgeInstancesInRange, ...
            'activeFridgeInstanceIdx', S.activeFridgeInstanceIdx, ...
            'isSwitchingFridge', S.isSwitchingFridge);
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
        S.fridgeInstancesInRange = savedTimeDomain.fridgeInstancesInRange;
        S.activeFridgeInstanceIdx = savedTimeDomain.activeFridgeInstanceIdx;
        S.isSwitchingFridge = savedTimeDomain.isSwitchingFridge;
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

        btnSnapshot.Enable = 'on';
        btnExport.Enable   = 'on';
        rebuildTimeline();
        recomputeSliderDataWindow();
        rebuildInstanceTimeline();
        updateSensorSelectionUi();
        onSyncSettingsChanged();
        setInstanceButtonsVisibility();
        configureTimeSlider();
        if ~isempty(targetStartTime) && (isnat(S.tNow) || isempty(S.tNow))
            updateAllPanesAtTime(targetStartTime, true);
        elseif ~isnat(S.tNow)
            updateAllPanesAtTime(S.tNow, true);
        end
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

    function toggleAdvancedControls()
        S.showAdvanced = ~S.showAdvanced;
        if S.showAdvanced
            advRow.Visible = 'on';
            btnAdvanced.Text = 'Hide Advanced';
        else
            advRow.Visible = 'off';
            btnAdvanced.Text = '[Advanced]';
        end
        setScanButtonsVisibility();
        setInstanceButtonsVisibility();
    end

    function onSyncSettingsChanged()
        S.syncMode = ddSyncMode.Value;
        S.snapMode = ddSnapMode.Value;
        S.toleranceMs = max(1, round(fldTolerance.Value));
        if strcmp(S.syncMode, 'FOLLOW_MASTER')
            ddMasterSensor.Enable = 'on';
            lblMaster.Visible = 'on';
            ddMasterSensor.Visible = 'on';
            S.masterSensor = ddMasterSensor.Value;
        else
            ddMasterSensor.Enable = 'off';
            lblMaster.Visible = 'off';
            ddMasterSensor.Visible = 'off';
        end
        updateDesyncUi();
    end

    function setInstanceButtonsVisibility()
        % Instance stepping is a core control and remains visible by default.
        btnPrevInstance.Visible = 'on';
        btnNextInstance.Visible = 'on';
        lblInstance.Visible = 'on';
    end

    function setScanButtonsVisibility()
        if isempty(hsiControlMap)
            return;
        end
        keys = hsiControlMap.keys;
        for ii = 1:numel(keys)
            ctrl = hsiControlMap(keys{ii});
            if isfield(ctrl,'btnPrev') && isgraphics(ctrl.btnPrev)
                ctrl.btnPrev.Visible = ternaryEnable(S.showAdvanced);
            end
            if isfield(ctrl,'btnNext') && isgraphics(ctrl.btnNext)
                ctrl.btnNext.Visible = ternaryEnable(S.showAdvanced);
            end
        end
    end

    function updateSensorSelectionUi()
        sensors = selectedSensorKeys();
        if isempty(sensors)
            ddMasterSensor.Items = {'LWIR'};
            ddMasterSensor.Value = 'LWIR';
            S.selectedSensors = {};
            S.masterSensor = '';
            btnPrevOverlap.Enable = 'off';
            btnNextOverlap.Enable = 'off';
            return;
        end
        S.selectedSensors = sensors;
        ddMasterSensor.Items = sensors;
        if isempty(S.masterSensor) || ~any(strcmp(sensors, S.masterSensor))
            if any(strcmp(sensors, 'LWIR'))
                S.masterSensor = 'LWIR';
            else
                S.masterSensor = sensors{1};
            end
        end
        ddMasterSensor.Value = S.masterSensor;
        sensorMap = buildOverlapSensorTimesMap();
        if sensorMap.Count > 0
            btnPrevOverlap.Enable = 'on';
            btnNextOverlap.Enable = 'on';
        else
            btnPrevOverlap.Enable = 'off';
            btnNextOverlap.Enable = 'off';
        end
    end

    function sensors = selectedSensorKeys()
        sensors = {};
        for si = 1:numel(modalities)
            mk = keyify(modalities{si});
            if getOr(S.exists, mk, false)
                sensors{end+1} = mk; %#ok<AGROW>
            end
        end
        hKeys = S.hsiGroupsMap.keys;
        for si = 1:numel(hKeys)
            sensors{end+1} = hKeys{si}; %#ok<AGROW>
        end
    end

    function jumpOverlap(direction)
        if isempty(S.playheadTs) || isnat(S.playheadTs)
            return;
        end
        sensorMap = buildOverlapSensorTimesMap();
        minSensors = minSensorsForOverlap(sensorMap);
        [tNext, info] = mv_find_next_overlap(S.playheadTs, direction, sensorMap, S.toleranceMs, minSensors, S.snapMode, S.masterSensor);
        if ~info.found || isempty(tNext)
            % Fallback: move to next/prev sample candidate even if strict overlap is absent.
            tNext = nextCandidateTime(S.playheadTs, direction, sensorMap, S.snapMode, S.masterSensor);
            if isempty(tNext)
                lblOverlapNote.Text = 'No further overlap found';
                return;
            end
            lblOverlapNote.Text = 'No strict overlap; moved to nearest next sample';
        elseif info.degraded
            lblOverlapNote.Text = 'No common overlap; snapped to nearest available';
        else
            lblOverlapNote.Text = '';
        end
        updateAllPanesAtTime(tNext, true);
    end

    function n = minSensorsForOverlap(sensorMap)
        keys = sensorMap.keys;
        countNonEmpty = 0;
        for ii = 1:numel(keys)
            t = sensorMap(keys{ii});
            if ~isempty(t)
                countNonEmpty = countNonEmpty + 1;
            end
        end
        n = min(2, max(1, countNonEmpty));
    end

    function t = nextCandidateTime(playheadTs, direction, sensorMap, snapMode, masterSensor)
        t = [];
        keys = sensorMap.keys;
        cands = datetime.empty(0,1);
        if strcmp(snapMode,'MASTER') && ~isempty(masterSensor) && isKey(sensorMap, masterSensor)
            cands = sensorMap(masterSensor);
        else
            for ii = 1:numel(keys)
                c = sensorMap(keys{ii});
                if isempty(c)
                    continue;
                end
                cands = [cands; c(:)]; %#ok<AGROW>
            end
        end
        if isempty(cands)
            return;
        end
        cands = unique(sort(cands));
        if direction >= 0
            idx = find(cands > playheadTs, 1, 'first');
        else
            idx = find(cands < playheadTs, 1, 'last');
        end
        if ~isempty(idx)
            t = cands(idx);
        end
    end

    function resyncPanels()
        if isempty(S.playheadTs) || isnat(S.playheadTs)
            return;
        end
        S.localScrubOverride = containers.Map('KeyType','char','ValueType','logical');
        S.lastActivePanelsSig = '';
        S.panelLocalTs = containers.Map('KeyType','char','ValueType','any');
        lblOverlapNote.Text = '';
        updateAllPanesAtTime(S.playheadTs, true);
    end

    function updateDesyncUi()
        panelTimes = {};
        if ~isempty(S.localScrubOverride)
            keys = S.localScrubOverride.keys;
            for ii = 1:numel(keys)
                k = keys{ii};
                if S.localScrubOverride(k) && isKey(S.panelLocalTs, k)
                    panelTimes{end+1} = S.panelLocalTs(k); %#ok<AGROW>
                end
            end
        end
        isDesync = mv_is_desynced(S.playheadTs, panelTimes, S.toleranceMs);
        if isDesync
            lblDesync.Text = 'DESYNC';
            btnResync.Visible = 'on';
            btnResync.Enable = 'on';
        else
            lblDesync.Text = '';
            btnResync.Visible = 'off';
            btnResync.Enable = 'off';
        end
    end


    function sensorMap = buildOverlapSensorTimesMap()
        % Overlap navigation treats all FRIDGE modalities as one sensor so
        % "Next Overlap" does not degenerate into per-frame stepping between
        % FRIDGE bands (e.g., SWIR vs MONO).
        sensorMap = containers.Map('KeyType','char','ValueType','any');

        fridgeCombined = datetime.empty(0,1);
        for si = 1:numel(modalities)
            m = keyify(modalities{si});
            if ~getOr(S.exists, m, false)
                continue;
            end
            t = fridgeTimesForModality(m, getOr(S.maxFrames, m, NaN));
            if isempty(t)
                continue;
            end
            fridgeCombined = [fridgeCombined; t(:)]; %#ok<AGROW>
        end
        if ~isempty(fridgeCombined)
            fridgeCombined = unique(sort(fridgeCombined(~isnat(fridgeCombined))));
            if ~isempty(fridgeCombined)
                sensorMap('FRIDGE') = fridgeCombined;
            end
        end

        hKeys = S.hsiGroupsMap.keys;
        for si = 1:numel(hKeys)
            t = S.hsiGroupsMap(hKeys{si}).timesUnique;
            if ~isempty(t)
                sensorMap(hKeys{si}) = t;
            end
        end
    end

    function sensorMap = buildSensorTimesMap()
        sensorMap = containers.Map('KeyType','char','ValueType','any');
        for si = 1:numel(modalities)
            m = keyify(modalities{si});
            if ~getOr(S.exists, m, false)
                continue;
            end
            t = fridgeTimesForModality(m, getOr(S.maxFrames, m, NaN));
            if ~isempty(t)
                sensorMap(m) = t;
            end
        end
        hKeys = S.hsiGroupsMap.keys;
        for si = 1:numel(hKeys)
            t = S.hsiGroupsMap(hKeys{si}).timesUnique;
            if ~isempty(t)
                sensorMap(hKeys{si}) = t;
            end
        end
    end

    function applySliderValue(val, isFinal)
        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime)
            return;
        end
        tTarget = S.sliderStartTime + seconds(val);
        if strcmp(S.snapMode,'MASTER') && ~isempty(S.masterSensor)
            sensorMap = buildSensorTimesMap();
            if isKey(sensorMap, S.masterSensor)
                [sampleTs, ~, ~] = mv_get_nearest_sample(sensorMap(S.masterSensor), tTarget);
                if ~isempty(sampleTs)
                    tTarget = sampleTs;
                end
            end
        end
        if ~isFinal
            S.pendingScrubTime = tTarget;
        end
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

    function tVec = syntheticFridgeTimes(maxF)
        tVec = datetime.empty(0,1);
        if isnat(S.fridgeStartTime) || isnat(S.fridgeEndTime)
            return;
        end
        if nargin < 1 || isempty(maxF) || isnan(maxF) || maxF < 1
            maxF = max(1, S.frameCount);
        end
        if maxF <= 1
            tVec = S.fridgeStartTime;
            return;
        end
        durSec = seconds(S.fridgeEndTime - S.fridgeStartTime);
        offsets = linspace(0, durSec, maxF);
        tVec = S.fridgeStartTime + seconds(offsets(:));
    end

    function tVec = fridgeTimesForModality(modality, maxF)
        modality = keyify(modality);
        tVec = getOr(S.fridgeTimesMap, modality, datetime.empty(0,1));
        if ~isdatetime(tVec)
            tVec = datetime.empty(0,1);
        end
        tVec = tVec(~isnat(tVec));

        tSynth = syntheticFridgeTimes(maxF);
        if isempty(tVec)
            tVec = tSynth;
            return;
        end

        % If parsed FRIDGE band times appear decoupled from the active
        % slider range (e.g., date mismatch), fall back to synthetic capture
        % times so scrubbing still advances FRIDGE frames synchronously.
        if ~isempty(tSynth) && ~isnat(S.sliderStartTime) && ~isnat(S.sliderEndTime)
            inRange = (tVec >= S.sliderStartTime) & (tVec <= S.sliderEndTime);
            if ~any(inRange)
                tVec = tSynth;
            end
        end
    end


    function rawPath = rawPathFromFridgeHdr(hdrPath)
        rawPath = '';
        if isempty(hdrPath)
            return;
        end
        [p,n,ext] = fileparts(hdrPath);
        if strcmpi(ext,'.hdr')
            rawPath = fullfile(p,[n '.raw']);
        else
            rawPath = strrep(hdrPath,'.hdr','.raw');
        end
    end

    function idxSel = bestFridgeInstanceIndexForTime(tNow)
        idxSel = NaN;
        insts = S.fridgeInstancesInRange;
        if isempty(insts)
            return;
        end
        nInst = numel(insts);
        starts = NaT(nInst,1);
        ends = NaT(nInst,1);
        for ii = 1:nInst
            starts(ii) = insts(ii).startTime;
            ends(ii) = insts(ii).endTime;
        end
        inRange = (tNow >= starts) & (tNow <= ends);
        idxIn = find(inRange, 1, 'first');
        if ~isempty(idxIn)
            idxSel = idxIn;
            return;
        end
        idxNext = find(starts >= tNow, 1, 'first');
        if ~isempty(idxNext)
            idxSel = idxNext;
            return;
        end
        validIdx = find(~isnat(starts), 1, 'last');
        if ~isempty(validIdx)
            idxSel = validIdx;
        end
    end

    function switched = switchFridgeInstanceForTime(tNow)
        switched = false;
        if S.isSwitchingFridge || isempty(S.fridgeInstancesInRange)
            return;
        end
        idxTarget = bestFridgeInstanceIndexForTime(tNow);
        if isempty(idxTarget) || isnan(idxTarget)
            return;
        end
        if ~isnan(S.activeFridgeInstanceIdx) && S.activeFridgeInstanceIdx == idxTarget
            return;
        end

        inst = S.fridgeInstancesInRange(idxTarget);
        if ~isfield(inst,'path') || isempty(inst.path)
            return;
        end
        rawPath = rawPathFromFridgeHdr(inst.path);
        if isempty(rawPath) || ~isfile(rawPath)
            return;
        end

        S.activeFridgeInstanceIdx = idxTarget;
        S.fridgeStartTime = inst.startTime;
        S.fridgeEndTime   = inst.endTime;

        S.isSwitchingFridge = true;
        c = onCleanup(@()setSwitchingFalse()); %#ok<NASGU>
        loadFromRawFile(rawPath);
        switched = true;
    end

    function setSwitchingFalse()
        S.isSwitchingFridge = false;
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
            if ~getOr(S.exists, m, false)
                continue;
            end
            maxF = getOr(S.maxFrames, m, NaN);
            tVec = fridgeTimesForModality(m, maxF);
            if isempty(tVec)
                continue;
            end
            tVec = tVec(tVec >= S.tStart & tVec <= S.tEnd);
            if ~isempty(tVec)
                tCollected = [tCollected; tVec(:)]; %#ok<AGROW>
            end
        end

        % Include FRIDGE instance intervals in-range so slider span can
        % cover multiple FRIDGE captures inside one timeline selection.
        if ~isempty(S.fridgeInstancesInRange)
            insts = S.fridgeInstancesInRange;
            for ii = 1:numel(insts)
                if isempty(insts(ii).startTime) || isempty(insts(ii).endTime)
                    continue;
                end
                t1 = max(insts(ii).startTime, S.tStart);
                t2 = min(insts(ii).endTime,   S.tEnd);
                if t2 >= t1
                    tCollected = [tCollected; t1; t2]; %#ok<AGROW>
                end
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
            lblSliderStart.Text = 'Start: -';
            lblSliderEnd.Text = 'End: -';
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

        lblSliderStart.Text = sprintf('Start: %s', formatTimeCursor(S.sliderStartTime));
        lblSliderEnd.Text = sprintf('End: %s', formatTimeCursor(S.sliderEndTime));
        lblSelectionRange.Text = sprintf('Selection: %s', formatClockRange(S.tStart, S.tEnd));
        if ~isnat(S.dataStartTime) && ~isnat(S.dataEndTime)
            lblDataSpan.Text = sprintf('Data span: %s', formatClockRange(S.dataStartTime, S.dataEndTime));
        else
            lblDataSpan.Text = 'Data span: (no data in selection)';
        end
        if ~isnat(S.tNow)
            updateSliderValueFromTime(S.tNow);
        end
    end

    function updateSliderValueFromTime(tNow)
        if isnat(S.sliderStartTime) || isnat(S.sliderEndTime) || isempty(tNow) || isnat(tNow)
            return;
        end
        val = seconds(tNow - S.sliderStartTime);
        if isempty(val) || ~isfinite(val)
            return;
        end
        lims = frameSlider.Limits;
        if numel(lims) < 2 || any(~isfinite(lims))
            return;
        end
        lo = lims(1);
        hi = lims(2);
        if hi < lo
            return;
        end
        val = min(max(val, lo), hi);
        sliderInternalUpdate = true;
        try
            frameSlider.Value = val;
        catch
        end
        sliderInternalUpdate = false;
    end

    function rebuildInstanceTimeline()
        S.instanceTimeline = struct('time', {}, 'type', {}, 'source', {}, 'ref', {});
        S.instanceTimes = datetime.empty(0,1);
        if isnat(S.tStart) || isnat(S.tEnd)
            updateInstanceControls();
            return;
        end

        entries = struct('time', {}, 'type', {}, 'source', {}, 'ref', {});

        if ~isempty(S.fridgeInstancesInRange)
            insts = S.fridgeInstancesInRange;
            for ii = 1:numel(insts)
                if ~isfield(insts(ii),'startTime') || isempty(insts(ii).startTime)
                    continue;
                end
                tVal = insts(ii).startTime;
                if isempty(tVal) || isnat(tVal)
                    continue;
                end
                if tVal < S.tStart || tVal > S.tEnd
                    continue;
                end
                if isfield(insts(ii),'modality') && ~isempty(insts(ii).modality)
                    source = keyify(insts(ii).modality);
                else
                    source = 'FRIDGE';
                end
                entries(end+1) = struct('time', tVal, 'type', 'FRIDGE', ...
                    'source', source, 'ref', ii); %#ok<AGROW>
            end
        end

        if S.enableHSI && ~isempty(S.hsiGroupsMap)
            keys = S.hsiGroupsMap.keys;
            for kk = 1:numel(keys)
                key = keys{kk};
                data = S.hsiGroupsMap(key);
                if ~isfield(data,'timesUnique') || isempty(data.timesUnique)
                    continue;
                end
                tVec = data.timesUnique(:);
                tVec = tVec(~isnat(tVec));
                tVec = tVec(tVec >= S.tStart & tVec <= S.tEnd);
                for jj = 1:numel(tVec)
                    entries(end+1) = struct('time', tVec(jj), 'type', 'HSI', ...
                        'source', key, 'ref', jj); %#ok<AGROW>
                end
            end
        end

        if isempty(entries)
            updateInstanceControls();
            return;
        end

        allTimes = [entries.time]';
        allTimes = unique(allTimes);
        allTimes = sort(allTimes);
        S.instanceTimes = allTimes;
        S.instanceTimeline = entries;
        updateInstanceControls();
    end

    function updateInstanceControls()
        if isempty(S.instanceTimes)
            btnPrevInstance.Enable = 'off';
            btnNextInstance.Enable = 'off';
            lblInstance.Text = 'Instance: -';
            return;
        end

        tNow = S.tNow;
        if isempty(tNow) || isnat(tNow)
            tNow = S.instanceTimes(1);
        end
        hasPrev = any(S.instanceTimes < tNow);
        hasNext = any(S.instanceTimes > tNow);
        btnPrevInstance.Enable = ternaryEnable(hasPrev);
        btnNextInstance.Enable = ternaryEnable(hasNext);

        [~, idx] = min(abs(S.instanceTimes - tNow));
        tNear = S.instanceTimes(idx);
        lblInstance.Text = sprintf('Instance: %s (%d/%d)', formatTimeCursor(tNear), idx, numel(S.instanceTimes));
    end

    function stepInstance(direction)
        tNow = S.playheadTs;
        if isempty(tNow) || isnat(tNow)
            tNow = S.tNow;
        end
        if isempty(tNow) || isnat(tNow)
            if ~isempty(S.instanceTimes)
                tNow = S.instanceTimes(1);
            else
                return;
            end
        end

        sensorMap = buildSensorTimesMap();
        if strcmp(S.syncMode,'FOLLOW_MASTER') || strcmp(S.snapMode,'MASTER')
            if isKey(sensorMap, S.masterSensor)
                cands = sort(sensorMap(S.masterSensor));
            else
                cands = S.instanceTimes;
            end
        elseif strcmp(S.snapMode,'ALL')
            [tNext, info] = mv_find_next_overlap(tNow, direction, sensorMap, S.toleranceMs, 2, 'ALL', S.masterSensor);
            if info.found && ~isempty(tNext)
                updateAllPanesAtTime(tNext, true);
            end
            return;
        else
            cands = S.instanceTimes;
        end

        if isempty(cands)
            return;
        end
        if direction > 0
            idx = find(cands > tNow, 1, 'first');
        else
            idx = find(cands < tNow, 1, 'last');
        end
        if isempty(idx)
            return;
        end
        updateAllPanesAtTime(cands(idx), true);
    end

    function e = ternaryEnable(tf)
        if tf
            e = 'on';
        else
            e = 'off';
        end
    end

    function updateFineSliderForTime(tNow)
        S.activeFineClip = struct('modality','', 'times', datetime.empty(0,1), 'startTime', NaT, 'endTime', NaT);
        fineSliderInternalUpdate = true;
        fineSlider.Enable = 'off';
        fineSlider.Limits = [1 2];
        fineSlider.Value = 1;
        lblFineFrame.Text = '-';

        if isempty(tNow) || isnat(tNow)
            fineSliderInternalUpdate = false;
            return;
        end

        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if ~getOr(S.exists, m, false)
                continue;
            end
            maxF = getOr(S.maxFrames, m, NaN);
            tVec = fridgeTimesForModality(m, maxF);
            if isempty(tVec)
                continue;
            end
            tVec = tVec(:);
            tVec = tVec(~isnat(tVec));
            if isempty(tVec)
                continue;
            end
            t1 = tVec(1);
            t2 = tVec(end);
            if tNow < t1 || tNow > t2
                continue;
            end
            [~, idx] = min(abs(tVec - tNow));
            nF = numel(tVec);
            fineSlider.Enable = 'on';
            fineSlider.Limits = [1 max(2,nF)];
            if isprop(fineSlider, 'SliderStep')
                fineSlider.SliderStep = [min(1/max(nF,1),0.25) min(10/max(nF,1),0.8)];
            end
            fineSlider.Value = idx;
            lblFineFrame.Text = sprintf('%s %d/%d', m, idx, nF);
            S.activeFineClip = struct('modality',m, 'times', tVec, 'startTime', t1, 'endTime', t2);
            break;
        end
        fineSliderInternalUpdate = false;
    end

    function fineSliderChanging(~, evt)
        if fineSliderInternalUpdate
            return;
        end
        applyFineSliderValue(evt.Value, false);
    end

    function fineSliderChanged(src, ~)
        if fineSliderInternalUpdate
            return;
        end
        applyFineSliderValue(src.Value, true);
    end

    function applyFineSliderValue(val, isFinal)
        if isempty(S.activeFineClip.times)
            return;
        end
        tVecActive = S.activeFineClip.times;
        idx = round(val);
        idx = min(max(1, idx), numel(tVecActive));
        tFrame = tVecActive(idx);

        % FRIDGE frame slider is local-first and propagates by timestamp
        % across available FRIDGE modalities so they scrub together.
        nUpdated = 0;
        for ii = 1:numel(modalities)
            m = keyify(modalities{ii});
            if ~getOr(S.exists, m, false)
                continue;
            end
            maxF = getOr(S.maxFrames, m, NaN);
            tVec = fridgeTimesForModality(m, maxF);
            if isempty(tVec)
                continue;
            end
            tVec = tVec(:);
            tVec = tVec(~isnat(tVec));
            if isempty(tVec)
                continue;
            end
            [sampleTs, ~, ~] = mv_get_nearest_sample(tVec, tFrame);
            if isempty(sampleTs) || isnat(sampleTs)
                continue;
            end
            S.panelLocalTs(m) = sampleTs;
            S.localScrubOverride(m) = true;
            nUpdated = nUpdated + 1;
        end

        if nUpdated > 0
            lblOverlapNote.Text = 'FRIDGE local scrub active (DESYNC)';
        end

        if isFinal
            S.renderJobId = S.renderJobId + 1;
        end
        tDraw = S.playheadTs;
        if isempty(tDraw) || isnat(tDraw)
            tDraw = S.tNow;
        end
        if isempty(tDraw) || isnat(tDraw)
            tDraw = tFrame;
        end
        drawAll(tDraw, S.renderJobId);
        updateFineSliderForTime(tFrame);
        updateDesyncUi();
    end

    function updateAllPanesAtTime(tNow, isFinal, syncPanels)
        if nargin < 2
            isFinal = true;
        end
        if nargin < 3
            syncPanels = true;
        end
        if isempty(tNow) || isnat(tNow) || isnat(S.sliderStartTime) || isnat(S.sliderEndTime)
            return;
        end
        tNow = clampTime(tNow);
        S.tNow = tNow;
        S.playheadTs = tNow;
        cancelObsoletePreloads();
        if syncPanels
            S.panelLocalTs = containers.Map('KeyType','char','ValueType','any');
            S.localScrubOverride = containers.Map('KeyType','char','ValueType','logical');
        end

        switched = switchFridgeInstanceForTime(tNow);
        if switched
            tNow = clampTime(S.tNow);
            S.tNow = tNow;
        end

        updateSliderValueFromTime(tNow);

        if hasFridgeTimes()
            idx = frameForTime(tNow);
            if ~isempty(idx) && ~isnan(idx)
                S.frame = idx;
            end
        end

        shouldDraw = true;
        if ~isFinal
            dt = toc(S.lastScrubPreviewTic);
            shouldDraw = dt >= (1 / max(1, S.scrubPreviewFps));
        end

        if shouldDraw || isFinal
            refreshMontageLayoutIfNeeded();
            S.renderJobId = S.renderJobId + 1;
            S.lastScrubPreviewTic = tic;
            if S.perfEnabled
                tDrawStart = tic;
            end
            drawAll(tNow, S.renderJobId);
            if S.perfEnabled
                fprintf('[perf] drawAll %.1f ms (final=%d)\n', toc(tDrawStart)*1000, isFinal);
            end
        end

        updateTimeDisplay();
        if isFinal || shouldDraw
            syncHsiToTime(tNow);
            updateInstanceControls();
            updateFineSliderForTime(tNow);
        end
        updateDesyncUi();
        preloadAroundTime(tNow, ~isFinal);

        if S.perfSyncDebug && isFinal
            req = S.frameCache.hits + S.frameCache.misses;
            hitRate = 0;
            if req > 0
                hitRate = 100 * (S.frameCache.hits / req);
            end
            fprintf('[sync/perf] cache hit %.1f%% (h=%d m=%d ev=%d aborted=%d\n', ...
                hitRate, S.frameCache.hits, S.frameCache.misses, S.frameCache.evictions, S.abortedRequests);
        end

        if isFinal && ~isnat(S.pendingScrubTime)
            S.pendingScrubTime = NaT;
        end
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
                    modality = normalizeFastModality(evt.modality);
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
            modality   = normalizeFastModality(extractAfter(sensorKey, 'FAST_'));
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
            dispName = sprintf('FAST %s', normalizeFastModality(extractAfter(sensorKey, 'FAST_')));
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
        setScanButtonsVisibility();
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

    function txt = setHsiTimestampText(sensorKey)
        txt = displayNameForSensor(sensorKey);
        if isKey(S.currentHsiMap, sensorKey)
            st = S.currentHsiMap(sensorKey);
            if isfield(st,'time') && isdatetime(st.time) && ~isnat(st.time)
                txt = sprintf('%s @ %s', txt, formatHsiTimestamp(st.time));
            end
        end
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
            allowHold = strcmp(S.syncMode,'FREE');
            if strcmp(S.syncMode,'FREE') && isKey(S.panelLocalTs, sensorKey)
                tResolve = S.panelLocalTs(sensorKey);
            else
                tResolve = tTarget;
            end
            res = mv_resolve_panel_sample(data.timesUnique, tResolve, struct('toleranceMs', S.toleranceMs, 'allowHold', allowHold));
            groupIdx = res.idx;
            if isnan(groupIdx) || strcmp(res.status,'NO_SAMPLE_NEARBY')
                setScanControlsState(sensorKey, NaN, NaN);
                setHsiFileLabel(sensorKey, sprintf('No data near playhead (Δt > %dms)', S.toleranceMs));
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
            if ~isempty(res.sample) && isdatetime(res.sample) && ~isnat(res.sample)
                S.panelLocalTs(sensorKey) = res.sample;
                setHsiFileLabel(sensorKey, sprintf('%s | Δt %+0.0fms', setHsiTimestampText(sensorKey), res.deltaMs));
            end
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
                hImg = imshow(ctx, [], 'Parent', cerbAxLWIR, 'InitialMagnification','fit');
                hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, 'CERB LWIR', ctx);
                hImg.PickableParts = 'all';
                hImg.HitTest       = 'on';
                [~,fnOnly,ext] = fileparts(fullpath);
                title(cerbAxLWIR, '');
                setHsiFileLabel('CERB_LWIR', [fnOnly ext]);
                S.cerb.LWIR = struct('path',fullpath,'ctx',ctx);
            case 'VNIR'
                hImg = imshow(ctx, [], 'Parent', cerbAxVNIR, 'InitialMagnification','fit');
                hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, 'CERB VNIR', ctx);
                hImg.PickableParts = 'all';
                hImg.HitTest       = 'on';
                [~,fnOnly,ext] = fileparts(fullpath);
                title(cerbAxVNIR, '');
                setHsiFileLabel('CERB_VNIR', [fnOnly ext]);
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
        hImg = imshow(ctx, [], 'Parent', mxAx, 'InitialMagnification','fit');
        hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, 'MX20 SW', ctx);
        hImg.PickableParts = 'all';
        hImg.HitTest       = 'on';
        [~,fnOnly,ext] = fileparts(hsicPath);
        title(mxAx, '');
        setHsiFileLabel('MX20', [fnOnly ext]);
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
        hImg = imshow(ctx, [], 'Parent', ax, 'InitialMagnification','fit');
        hImg.ButtonDownFcn = @(src,evt) onPixelClickCerb(src, evt, ['FAST ' key], ctx);
        hImg.PickableParts = 'all';
        hImg.HitTest       = 'on';
        [~,fnOnly,ext] = fileparts(hsicPath);
        title(ax, '');
        setHsiFileLabel(sprintf('FAST_%s', key), [fnOnly ext]);
        S.fast.(key) = struct('hdr',hdrOrHsicPath,'ctx',ctx);

        if doLayout
            refreshMontageLayout();
        end
    end

    function key = frameCacheKey(modality, frameIdx)
        key = sprintf('%s#%d', keyify(modality), frameIdx);
    end

    function [img, cacheHit, ok, errText] = getFrameWithCache(modality, frameIdx)
        % SWrV cache lookup: render previous image until this returns.
        key = frameCacheKey(modality, frameIdx);
        [S.frameCache, img, cacheHit] = mv_lru_get(S.frameCache, key);
        if cacheHit
            ok = true;
            errText = '';
            return;
        end

        ok = false;
        errText = '';
        tLoad = tic;
        try
            img = fridge_read_frame(modality, frameIdx, S.hdrs, S.files);
            S.frameCache = mv_lru_put(S.frameCache, key, img);
            ok = true;
        catch ME
            img = [];
            errText = ME.message;
        end
        if S.perfSyncDebug
            fprintf('[sync/perf] frame load %s f=%d cacheHit=%d %.1fms\n', modality, frameIdx, cacheHit, toc(tLoad)*1000);
        end
    end

    function cancelObsoletePreloads()
        keys = S.prefetchLoading.keys;
        for ii = 1:numel(keys)
            if S.prefetchLoading(keys{ii})
                S.abortedRequests = S.abortedRequests + 1;
            end
        end
        S.prefetchLoading = containers.Map('KeyType','char','ValueType','logical');
    end

    function preloadAroundTime(tNow, isScrubbing)
        if isempty(tNow) || isnat(tNow)
            return;
        end
        preloadSpan = 1;
        if ~isScrubbing
            preloadSpan = 2;
        end
        for i = 1:numel(modalities)
            m = keyify(modalities{i});
            if ~getOr(S.exists, m, false)
                continue;
            end
            [idxNow, statusNow] = effectiveFrame(m, tNow);
            if isnan(idxNow) || strcmp(statusNow, 'no-nearby')
                continue;
            end
            maxF = getOr(S.maxFrames, m, NaN);
            if isnan(maxF)
                continue;
            end
            offsets = [-preloadSpan:-1 1:preloadSpan];
            for oi = 1:numel(offsets)
                idx = idxNow + offsets(oi);
                if idx < 1 || idx > maxF
                    continue;
                end
                key = frameCacheKey(m, idx);
                if isKey(S.frameCache.map, key)
                    continue;
                end
                if isKey(S.prefetchLoading, key) && S.prefetchLoading(key)
                    continue;
                end
                S.prefetchLoading(key) = true;
                [~,~,ok,~] = getFrameWithCache(m, idx);
                S.prefetchLoading(key) = false;
                if ~ok && isKey(S.prefetchLoading, key)
                    remove(S.prefetchLoading, key);
                end
            end
        end
    end

    %======================== DRAWING / IO (FRIDGE) ========================
    function drawAll(tCurrent, jobId)
        if nargin < 1 || isempty(tCurrent) || isnat(tCurrent)
            tCurrent = S.tNow;
        end
        if nargin < 2
            jobId = S.renderJobId;
        end
        for i = 1:numel(modalities)
            if jobId ~= S.renderJobId
                return;
            end
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
                frameLbl.Text = sprintf('%s — Missing metadata', fridgeDisplayName(m));
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
                frameLbl.Text = sprintf('%s — Missing file', fridgeDisplayName(m));
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
                frameLbl.Text = sprintf('%s — Missing header', fridgeDisplayName(m));
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            [fEff, status, deltaMs, sampleTs] = effectiveFrame(m, tCurrent);
            newStatus = status;
            newIdx = fEff;
            if isequal(prevIdx, newIdx) && strcmp(prevStatus, newStatus)
                continue;
            end

            maxF = getOr(S.maxFrames, m, NaN);
            if isnan(fEff) || strcmp(status,'no-nearby')
                cla(ax);
                axis(ax,'off');
                title(ax, '');
                if strcmp(status,'no-nearby')
                    frameLbl.Text = sprintf('%s — No data near playhead (Δt > %dms)', fridgeDisplayName(m), S.toleranceMs);
                else
                    frameLbl.Text = sprintf('%s — No frame (max %g)', fridgeDisplayName(m), maxF);
                end
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                S.lastFrameByMod(m) = newIdx;
                S.lastStatusByMod(m) = newStatus;
                continue;
            end

            [img, cacheHit, loadOk, errText] = getFrameWithCache(m, fEff);
            if ~loadOk
                % Stale-while-revalidate: keep prior panel image on load errors.
                frameLbl.Text = sprintf('%s — load failed (showing previous)', fridgeDisplayName(m));
                if S.perfSyncDebug
                    fprintf('[sync/perf] load error %s f=%d: %s\n', m, fEff, errText);
                end
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                continue;
            end

            frameLbl.Text = sprintf('%s — Loading...', fridgeDisplayName(m));
            drawnow limitrate nocallbacks;
            cla(ax);
            title(ax, '');
            fileLbl.Text  = '';

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
                hImg = imshow(imgDisp, 'Parent', ax, 'InitialMagnification','fit');
            else
                hImg = imshow(img, [], 'Parent', ax, 'InitialMagnification','fit');
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
            if ~isempty(sampleTs) && isdatetime(sampleTs) && ~isnat(sampleTs)
                S.panelLocalTs(m) = sampleTs;
            end
            cacheNote = '';
            if cacheHit
                cacheNote = ' [cache]';
            end
            frameLbl.Text = sprintf('%s — %d/%d%s | Δt %+0.0fms%s', fridgeDisplayName(m), fEff, maxF, note, deltaMs, cacheNote);
            [~,fn,ext] = fileparts(getOr(S.files, m, ''));
            fileLbl.Text = [fn ext];
            S.lastFrameByMod(m) = newIdx;
            S.lastStatusByMod(m) = newStatus;
        end
    end

    function [fEff, status, deltaMs, sampleTs] = effectiveFrame(modality, targetTime)
        modality = keyify(modality);
        if nargin < 2
            targetTime = [];
        end

        maxF = getOr(S.maxFrames, modality, NaN);
        if isnan(maxF) || maxF < 1
            fEff   = NaN;
            status = 'missing';
            deltaMs = NaN;
            sampleTs = NaT;
            return;
        end

        tVec = fridgeTimesForModality(modality, maxF);
        if isempty(tVec) || ~isdatetime(tVec)
            fEff   = NaN;
            status = 'missing';
            deltaMs = NaN;
            sampleTs = NaT;
            return;
        end

        if isempty(targetTime) || ~isdatetime(targetTime)
            targetTime = S.playheadTs;
        end

        allowHold = strcmp(S.syncMode,'FREE');
        useLocalOverride = isKey(S.localScrubOverride, modality) && S.localScrubOverride(modality) && isKey(S.panelLocalTs, modality);
        if useLocalOverride
            targetResolve = S.panelLocalTs(modality);
            allowHold = true;
        elseif strcmp(S.syncMode,'FREE') && isKey(S.panelLocalTs, modality)
            targetResolve = S.panelLocalTs(modality);
            allowHold = true;
        else
            targetResolve = targetTime;
        end

        res = mv_resolve_panel_sample(tVec, targetResolve, struct('toleranceMs', S.toleranceMs, 'allowHold', allowHold));
        fEff = res.idx;
        sampleTs = res.sample;
        if useLocalOverride && ~isempty(S.playheadTs) && ~isnat(S.playheadTs) && ~isempty(sampleTs) && ~isnat(sampleTs)
            deltaMs = milliseconds(sampleTs - S.playheadTs);
        else
            deltaMs = res.deltaMs;
        end
        if strcmp(res.status, 'OK')
            status = 'ok';
        elseif strcmp(res.status, 'HOLDING_LAST')
            status = 'held';
        else
            status = 'no-nearby';
        end
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
        if isfield(initial,'fridgeInstancesInRange') && ~isempty(initial.fridgeInstancesInRange)
            S.fridgeInstancesInRange = initial.fridgeInstancesInRange;
        else
            S.fridgeInstancesInRange = struct([]);
        end
        S.activeFridgeInstanceIdx = NaN;
        [S.tStart, S.tEnd] = resolveTimeDomain(initial);
        if ~isnat(S.tStart)
            S.tNow = S.tStart;
        else
            S.tNow = NaT;
        end
        S.sliderStartTime = S.tStart;
        S.sliderEndTime   = S.tEnd;
        rebuildHsiGroups();
        rebuildInstanceTimeline();
        updateReturnButtonState();

        if isfield(initial,'rawFile') && ~isempty(initial.rawFile) && isfile(initial.rawFile)
            loadFromRawFile(initial.rawFile);
            if ~isempty(S.fridgeInstancesInRange)
                for ii = 1:numel(S.fridgeInstancesInRange)
                    candRaw = rawPathFromFridgeHdr(S.fridgeInstancesInRange(ii).path);
                    if strcmpi(candRaw, initial.rawFile)
                        S.activeFridgeInstanceIdx = ii;
                        break;
                    end
                end
            end
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
        rebuildInstanceTimeline();
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
        S.fridgeInstancesInRange = struct([]);
        S.activeFridgeInstanceIdx = NaN;
        S.isSwitchingFridge = false;
        S.hsiEvents    = struct('sensor', {}, 'time', {}, 'path', {}, 'modality', {});
        S.hsiGroupsMap = containers.Map('KeyType','char','ValueType','any');
        S.currentHsi   = struct('sensor','', 'time', NaT, 'effectiveTime', NaT);
        S.hsiPreciseCache = containers.Map('KeyType','char','ValueType','any');
        S.enableHSI = enableHSI;
        S.currentHsiMap = containers.Map('KeyType','char','ValueType','any');
        S.instanceTimeline = struct('time', {}, 'type', {}, 'source', {}, 'ref', {});
        S.instanceTimes = datetime.empty(0,1);
        S.activeFineClip = struct('modality','', 'times', datetime.empty(0,1), 'startTime', NaT, 'endTime', NaT);

        lblStatus.Text = 'Status: (no capture loaded)';
        lblFrames.Text = 'Frames: -';
        lblMem.Text    = '';
        lblPixel.Text  = 'Pixel: -';
        lblValue.Text  = 'Value: -';
        lblTime.Text   = 'Time: -';
        lblSelectionRange.Text = 'Selection: -';
        lblDataSpan.Text = 'Data span: -';
        lblSliderStart.Text = 'Start: -';
        lblSliderEnd.Text = 'End: -';

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
                frameLbl.Text = sprintf('%s — -', fridgeDisplayName(modName));
            end

            fileLbl = getOr(fileLabelMap, modName, []);
            if ~isempty(fileLbl)
                fileLbl.Text  = '';
            end
        end

        btnSnapshot.Enable = 'off';
        btnExport.Enable   = 'off';
        btnPrevInstance.Enable = 'off';
        btnNextInstance.Enable = 'off';
        lblInstance.Text = 'Instance: -';

        frameSlider.Enable = 'off';
        frameSlider.Limits = [1 2];
        frameSlider.Value  = 1;

        fineSlider.Enable = 'off';
        fineSlider.Limits = [1 2];
        fineSlider.Value  = 1;
        lblFineFrame.Text = '-';

        cla(cerbAxLWIR); title(cerbAxLWIR,'');
        cla(cerbAxVNIR); title(cerbAxVNIR,'');
        cla(mxAx);       title(mxAx,'');
        setHsiFileLabel('CERB_LWIR', '');
        setHsiFileLabel('CERB_VNIR', '');
        setHsiFileLabel('MX20', '');
        if ~isempty(fastAxes)
            keys = fastAxes.keys;
            for ii = 1:numel(keys)
                ax = fastAxes(keys{ii});
                if isgraphics(ax)
                    cla(ax);
                    title(ax, '');
                    setHsiFileLabel(sprintf('FAST_%s', keys{ii}), '');
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
