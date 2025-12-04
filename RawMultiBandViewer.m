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

    %======================== UI LAYOUT ===================================
    f = uifigure('Name','AARO Multi-Band Viewer','Position',[80 80 1280 900]);

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

    % Normalize map keys so callers can provide either char or string
    % modality names without triggering containers.Map indexing errors.
    labelCtor = pickLabelCtor();
    makeLabel = @(parent, varargin) labelCtor(parent, varargin{:});

    function ctor = pickLabelCtor()
        % Prefer the class constructor, then the helper function, with a
        % final fallback to uicontrol text for legacy runtimes.
        if exist('matlab.ui.control.Label','class')
            ctor = @(parent, varargin) matlab.ui.control.Label(parent, varargin{:});
            return;
        end
        if exist('uilabel','builtin') || exist('uilabel','file')
            ctor = @uilabel;
            return;
        end
        ctor = @fallbackLabel;
    end

    function lbl = fallbackLabel(parent, varargin)
        lbl = uicontrol(parent, 'Style','text');
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
    hsiPanel = uipanel(hiddenBin);
    hsiGrid = uigridlayout(hsiPanel,[2,1]);
    hsiGrid.RowHeight   = {'fit','1x'};
    hsiGrid.ColumnWidth = {'1x'};

    makeLabel(hsiGrid, 'Text','HSI Context', 'FontWeight','bold', ...
        'HorizontalAlignment','center').Layout.Row = 1;

    cerbTabs = uitabgroup(hsiGrid);
    cerbTabs.Layout.Row    = 2;
    cerbTabs.Layout.Column = 1;

    tabLWIR = uitab(cerbTabs,'Title','CERB LWIR');
    tabVNIR = uitab(cerbTabs,'Title','CERB VNIR');
    tabMX20 = uitab(cerbTabs,'Title','MX20 SW');
    
    % --- CERB LWIR tab: 1x1 grid, axes fills whole tab ---
    tabLWIRGrid = uigridlayout(tabLWIR,[1 1]);
    tabLWIRGrid.RowHeight   = {'1x'};
    tabLWIRGrid.ColumnWidth = {'1x'};
    
    cerbAxLWIR = uiaxes(tabLWIRGrid);
    cerbAxLWIR.Layout.Row    = 1;
    cerbAxLWIR.Layout.Column = 1;
    axis(cerbAxLWIR,'off');
    title(cerbAxLWIR,'CERB LWIR');
    
    % --- CERB VNIR tab ---
    tabVNIRGrid = uigridlayout(tabVNIR,[1 1]);
    tabVNIRGrid.RowHeight   = {'1x'};
    tabVNIRGrid.ColumnWidth = {'1x'};
    
    cerbAxVNIR = uiaxes(tabVNIRGrid);
    cerbAxVNIR.Layout.Row    = 1;
    cerbAxVNIR.Layout.Column = 1;
    axis(cerbAxVNIR,'off');
    title(cerbAxVNIR,'CERB VNIR');
    
    % --- MX20 tab ---
    tabMX20Grid = uigridlayout(tabMX20,[1 1]);
    tabMX20Grid.RowHeight   = {'1x'};
    tabMX20Grid.ColumnWidth = {'1x'};
    
    mxAx = uiaxes(tabMX20Grid);
    mxAx.Layout.Row    = 1;
    mxAx.Layout.Column = 1;
    axis(mxAx,'off');
    title(mxAx,'MX20 SW');

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

    navTop = uigridlayout(navCol,[1,1]);
    navTop.ColumnWidth = {'fit'};
    btnSave = uibutton(navTop,'Text','Save Montage PNG','Enable','off', ...
        'ButtonPushedFcn',@(~,~)saveMontage());

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
    timeCol = uigridlayout(ctrlWrapper,[2,1]);
    timeCol.RowHeight   = {'fit','fit'};
    timeCol.ColumnWidth = {'1x'};
    makeLabel(timeCol,'Text','Current time','FontWeight','bold', ...
        'HorizontalAlignment','left');
    lblTime = makeLabel(timeCol,'Text','Time: -','HorizontalAlignment','left');

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
    S.hsiEvents    = struct('sensor', {}, 'time', {}, 'path', {});
    S.currentHsi   = struct('sensor','', 'time', NaT, 'effectiveTime', NaT);
    S.hsiPreciseCache = containers.Map('KeyType','char','ValueType','any');
    S.timelineFig  = [];

    S.cerb = struct('LWIR',[],'VNIR',[]);
    S.mx20 = struct('hdr',[],'ctx',[]);

    targetStartTime = [];
    if isfield(initial,'initialTime')
        targetStartTime = initial.initialTime;
    end
    if isfield(initial,'hsiEvents')
        S.hsiEvents = initial.hsiEvents;
    end
    if isfield(initial,'timelineFig') && ~isempty(initial.timelineFig) && isvalid(initial.timelineFig)
        S.timelineFig = initial.timelineFig;
    end
    updateHsiJumpAvailability();
    sliderInternalUpdate = false;  % prevent recursive slider callbacks
    updateReturnButtonState();

    %======================== LAYOUT HELPERS ==============================
    function tf = hasCerb(whichMod)
        tf = isfield(S, 'cerb') && isfield(S.cerb, whichMod) && ...
             ~isempty(S.cerb.(whichMod));
    end

    function tf = hasMx20()
        tf = isfield(S, 'mx20') && isfield(S.mx20, 'hdr') && ...
             ~isempty(S.mx20.hdr);
    end

    function tf = hasAnyHsi()
        tf = hasCerb('LWIR') || hasCerb('VNIR') || hasMx20();
    end

    function updateHsiTabVisibility()
        if hasCerb('LWIR')
            tabLWIR.Visible = 'on';
        else
            tabLWIR.Visible = 'off';
        end

        if hasCerb('VNIR')
            tabVNIR.Visible = 'on';
        else
            tabVNIR.Visible = 'off';
        end

        if hasMx20()
            tabMX20.Visible = 'on';
        else
            tabMX20.Visible = 'off';
        end

        if hasAnyHsi()
            if strcmp(tabLWIR.Visible,'on')
                cerbTabs.SelectedTab = tabLWIR;
            elseif strcmp(tabVNIR.Visible,'on')
                cerbTabs.SelectedTab = tabVNIR;
            else
                cerbTabs.SelectedTab = tabMX20;
            end
        else
            cerbTabs.SelectedTab = [];
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
        resetUI();
        S.hsiEvents = existingEvents;

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
        btnSave.Enable     = 'on';
        updateHsiJumpAvailability();

        rebuildTimeline();

        if ~isempty(targetStartTime)
            jumpToTime(targetStartTime);
        else
            drawAll();
            updateTimeDisplay();
            syncHsiToTime(timeForFrame(S.frame));
        end
    end

    function step(delta)
        if S.nFrames < 1
            return;
        end
        S.frame = min(max(1, S.frame + delta), S.nFrames);
        setSliderFromFrame();
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
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
        % Build a combined time axis so the slider spans the earliest to
        % latest FRIDGE timestamps across all modalities.
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
                isfield(initial,'fridgeStartTime') && ...
                isfield(initial,'fridgeEndTime') && ...
                ~isempty(initial.fridgeStartTime) && ...
                ~isempty(initial.fridgeEndTime)
            % Synthesize evenly spaced times when only capture bounds exist.
            durSec = seconds(initial.fridgeEndTime - initial.fridgeStartTime);
            steps  = max(2, max(1, S.frameCount));
            offsets = linspace(0, durSec, steps);
            allTimes = initial.fridgeStartTime + seconds(offsets(:));
        end

        if ~isempty(allTimes)
            allTimes       = unique(allTimes);
            allTimes       = sort(allTimes);
            S.timelineTimes = allTimes;
            S.nFrames      = numel(S.timelineTimes);

            S.sliderMode = 'time';
            S.sliderOrigin = S.timelineTimes(1);
            sliderLimits = seconds(S.timelineTimes([1 end]) - S.sliderOrigin);
            if sliderLimits(1) == sliderLimits(2)
                sliderLimits(2) = sliderLimits(2) + eps;
            end
            frameSlider.Limits = sliderLimits;
            frameSlider.Value  = sliderLimits(1);
            frameSlider.Enable = 'on';
        else
            S.sliderMode = 'frame';
            S.sliderOrigin = NaT;
            S.nFrames    = max(1, S.frameCount);

            if S.nFrames <= 1
                frameSlider.Limits = [1 2];
            else
                frameSlider.Limits = [1 S.nFrames];
            end
            frameSlider.Value  = 1;
            frameSlider.Enable = 'on';
        end

        setSliderFromFrame();
    end

    function setSliderFromFrame()
        sliderInternalUpdate = true;
        if strcmp(S.sliderMode,'time') && hasFridgeTimes()
            frameSlider.Value = seconds(timeForFrame(S.frame) - S.sliderOrigin);
        else
            frameSlider.Value = S.frame;
        end
        sliderInternalUpdate = false;
    end

    function jumpToTime(tTarget)
        if isempty(tTarget)
            return;
        end

        idx = frameForTime(tTarget);
        if isempty(idx) || isnan(idx)
            % No FRIDGE time data: still try to sync HSI to the target time
            syncHsiToTime(tTarget);
            updateTimeDisplay();
            return;
        end

        S.frame = idx;
        setSliderFromFrame();
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
    end

    function applySliderValue(val, isFinal)
        if strcmp(S.sliderMode,'time') && hasFridgeTimes()
            targetTime = S.sliderOrigin + seconds(val);
            if targetTime < S.timelineTimes(1)
                targetTime = S.timelineTimes(1);
            elseif targetTime > S.timelineTimes(end)
                targetTime = S.timelineTimes(end);
            end
            idx = frameForTime(targetTime);
            if isempty(idx) || isnan(idx)
                return;
            end
            if idx == S.frame && ~isFinal
                return;
            end
            S.frame = idx;
            sliderInternalUpdate = true;
            frameSlider.Value = seconds(timeForFrame(S.frame) - S.sliderOrigin);
            sliderInternalUpdate = false;
            drawAll();
            updateTimeDisplay();
            syncHsiToTime(timeForFrame(S.frame));
            return;
        end

        if S.nFrames < 1
            return;
        end
        newFrame = round(val);
        newFrame = max(1, min(S.nFrames, newFrame));
        if newFrame == S.frame && ~isFinal
            return;
        end
        S.frame = newFrame;
        sliderInternalUpdate = true;
        frameSlider.Value = S.frame;
        sliderInternalUpdate = false;
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
    end

    function syncHsiToTime(tTarget)
        if isempty(S.hsiEvents)
            return;
        end

        if nargin < 1 || isempty(tTarget)
            effTimesTmp = arrayfun(@(e) effectiveHsiTime(e), S.hsiEvents);
            effTimesTmp = effTimesTmp(~isnat(effTimesTmp));
            if isempty(effTimesTmp)
                return;
            end
            tTarget = min(effTimesTmp);
        end

        effTimes = arrayfun(@(e) effectiveHsiTime(e), S.hsiEvents);
        diffs    = abs(effTimes - tTarget);
        [~, idxEvt] = min(diffs);
        evt = S.hsiEvents(idxEvt);
        evtEff = effTimes(idxEvt);

        if strcmp(S.currentHsi.sensor, evt.sensor) && ...
           isdatetime(S.currentHsi.effectiveTime) && S.currentHsi.effectiveTime == evtEff
            return;
        end

        switch evt.sensor
            case 'CERB'
                loadCerbFromPath('LWIR', evt.path);
            case 'MX20'
                loadMX20FromHdr(evt.path);
        end

        S.currentHsi = struct('sensor', evt.sensor, 'time', evt.time, ...
                              'effectiveTime', evtEff);
    end

    function tEff = effectiveHsiTime(evt)
        % Use per-file unixtime when present to distinguish closely spaced
        % HSI scans that share the same coarse timestamp.
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

    function updateHsiJumpAvailability()
        if isempty(S.hsiEvents)
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

    function saveMontage()
        % Uses fridge_read_frame + your toUint8/padTo helpers
        names = {'LWIR','MWIR','SWIR'; 'MONO','VIS-COLOR',''};
        tiles = cell(2,3);
        maxH  = 0;
        maxW  = 0;

        for r = 1:2
            for c = 1:3
                name = keyify(names{r,c});
                if isempty(name)
                    tiles{r,c} = uint8(255);
                    continue;
                end
                if ~getOr(S.exists, name, false)
                    tile = uint8(255*ones(100,160,'uint8'));
                else
                    [fEff, ~] = effectiveFrame(name, S.frame, timeForFrame(S.frame));
                    if isnan(fEff)
                        tile = uint8(255*ones(100,160,'uint8'));
                    else
                        tile = toUint8(fridge_read_frame(name, fEff, S.hdrs, S.files));
                    end
                end
                tiles{r,c} = tile;
                [h,w] = size(tile);
                maxH = max(maxH,h);
                maxW = max(maxW,w);
            end
        end

        for r = 1:2
            for c = 1:3
                tiles{r,c} = padTo(tiles{r,c}, [maxH,maxW]);
            end
        end

        row1 = [tiles{1,1}, tiles{1,2}, tiles{1,3}];
        row2 = [tiles{2,1}, tiles{2,2}];
        M    = [row1; row2];

        [file, path] = uiputfile({'*.png'}, 'Save Montage as PNG', ...
                                 fullfile(S.dir, sprintf('montage_frame_%d.png', S.frame)));
        if isequal(file,0)
            return;
        end
        imwrite(M, fullfile(path,file));
        uialert(f,'Montage saved.','Done');
    end

    %======================== CERBERUS / MX20 ==============================
    function loadCerbFromPath(whichMod, fullpath)
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

        refreshMontageLayout();
    end

    function loadMX20FromHdr(hdrOrHsicPath)
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

        refreshMontageLayout();
    end

    %======================== DRAWING / IO (FRIDGE) ========================
    function drawAll()
        tCurrent = timeForFrame(S.frame);
        for i = 1:numel(modalities)
            m  = keyify(modalities{i});
            ax = axMap(m);
            frameLbl = frameLabelMap(m);
            fileLbl  = fileLabelMap(m);
            cla(ax);
            title(ax, '');
            frameLbl.Text = 'Frame: -';
            fileLbl.Text  = '';

            if ~hasKey(S.exists, m) || ~hasKey(S.hdrs, m) || ~hasKey(S.files, m)
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — Missing metadata', m);
                continue;
            end

            if ~getOr(S.exists, m, false)
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — Missing file', m);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                continue;
            end

            hdrCandidate = getOr(S.hdrs, m, []);
            if isempty(hdrCandidate)
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — Missing header', m);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                continue;
            end

            [fEff, status] = effectiveFrame(m, S.frame, tCurrent);
            maxF = getOr(S.maxFrames, m, NaN);
            if isnan(fEff)
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — No frame (max %g)', m, maxF);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
                continue;
            end

            try
                img = fridge_read_frame(m, fEff, S.hdrs, S.files);
            catch ME
                axis(ax,'off');
                frameLbl.Text = sprintf('%s — %s', m, ME.message);
                [~,fn,ext] = fileparts(getOr(S.files, m, ''));
                fileLbl.Text = [fn ext];
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
        end
    end

    function [fEff, status] = effectiveFrame(modality, fReq, targetTime)
        modality = keyify(modality);
        if nargin < 3
            targetTime = [];
        end

        maxF = getOr(S.maxFrames, modality, NaN);
        if isnan(maxF) || maxF < 1
            fEff   = NaN;
            status = 'missing';
            return;
        end

        % Time-driven: pick the nearest timestamp for this modality and hold
        % the first/last frame outside its bounds so all panes stay aligned.
        if strcmp(S.sliderMode,'time') && hasFridgeTimes() && ...
                hasKey(S.fridgeTimesMap, modality)
            tVec = getOr(S.fridgeTimesMap, modality, datetime.empty(0,1));
            if isdatetime(tVec) && ~isempty(tVec)
                tVec = tVec(~isnat(tVec));
                if isempty(targetTime)
                    targetTime = timeForFrame(fReq);
                end
                if ~isempty(targetTime) && isdatetime(targetTime)
                    if targetTime <= tVec(1)
                        idxSel = 1;
                        status = 'held';
                    elseif targetTime >= tVec(end)
                        idxSel = numel(tVec);
                        status = 'held';
                    else
                        [~, idxSel] = min(abs(tVec - targetTime));
                        status = 'ok';
                    end
                    idxSel = min(max(1, idxSel), numel(tVec));
                    if ~isnan(maxF)
                        idxSel = min(idxSel, maxF);
                    end
                    fEff = idxSel;
                    return;
                end
            end
        end

        switch S.behavior
            case 'hold'
                if fReq > maxF
                    fEff   = maxF;
                    status = 'held';
                else
                    fEff   = fReq;
                    status = 'ok';
                end
            case 'missing'
                if fReq > maxF
                    fEff   = NaN;
                    status = 'missing';
                else
                    fEff   = fReq;
                    status = 'ok';
                end
            case 'loop'
                if fReq <= maxF
                    fEff   = fReq;
                    status = 'ok';
                else
                    fEff   = mod(fReq-1, maxF) + 1;
                    status = 'looped';
                end
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
        if ~hasFridgeTimes() || numel(S.timelineTimes) < S.frame
            lblTime.Text = 'Time: (no FRIDGE time data)';
            return;
        end
        t = timeForFrame(S.frame);
        lblTime.Text = sprintf('Time: %s', datestr(t,'yyyy-mm-dd HH:MM:SS.FFF'));
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

    %======================== RESET / INITIAL LOAD ========================
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
        S.hsiEvents    = struct('sensor', {}, 'time', {}, 'path', {});
        S.currentHsi   = struct('sensor','', 'time', NaT, 'effectiveTime', NaT);
        S.hsiPreciseCache = containers.Map('KeyType','char','ValueType','any');

        lblStatus.Text = 'Status: (no capture loaded)';
        lblFrames.Text = 'Frames: -';
        lblMem.Text    = '';
        lblPixel.Text  = 'Pixel: -';
        lblValue.Text  = 'Value: -';
        lblTime.Text   = 'Time: -';

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
        btnSave.Enable     = 'off';
        btnJumpHsi.Enable  = 'off';

        frameSlider.Enable = 'off';
        frameSlider.Limits = [1 2];
        frameSlider.Value  = 1;

        cla(cerbAxLWIR); title(cerbAxLWIR,'CERB LWIR');
        cla(cerbAxVNIR); title(cerbAxVNIR,'CERB VNIR');
        cla(mxAx);       title(mxAx,'MX20 SW');

        S.cerb = struct('LWIR',[],'VNIR',[]);
        S.mx20 = struct('hdr',[],'ctx',[]);

        refreshMontageLayout();
        updateReturnButtonState();
    end

    % Initial auto-load from timeline
    if isfield(initial,'rawFile') && ~isempty(initial.rawFile) && isfile(initial.rawFile)
        loadFromRawFile(initial.rawFile);
    end
    if isfield(initial,'cerbLWIR') && ~isempty(initial.cerbLWIR) && isfile(initial.cerbLWIR)
        loadCerbFromPath('LWIR', initial.cerbLWIR);
    end
    if isfield(initial,'cerbVNIR') && ~isempty(initial.cerbVNIR) && isfile(initial.cerbVNIR)
        loadCerbFromPath('VNIR', initial.cerbVNIR);
    end
    if isfield(initial,'mx20Hdr') && ~isempty(initial.mx20Hdr) && isfile(initial.mx20Hdr)
        loadMX20FromHdr(initial.mx20Hdr);
    end

    updateTimeDisplay();
    if ~isempty(targetStartTime)
        syncHsiToTime(targetStartTime);
    else
        syncHsiToTime(timeForFrame(S.frame));
    end

end
