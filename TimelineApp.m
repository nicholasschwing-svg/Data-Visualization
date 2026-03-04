function TimelineApp()
    % TimelineApp
    % Timeline UI with:
    % - Date dropdown
    % - CERBERUS dots + MX20 dots + FRIDGE bars
    % - Click on/near a dot or bar to see a popup
    % - Adaptive HH:MM ticks based on current zoom
    % - Drag-rectangle selection to choose a time range.
    %
    % Vertical position of the selection box controls which HSI sensor:
    %   - near CERB row -> CERBERUS only
    %   - near MX20 row -> MX20 only
    %   - otherwise     -> both allowed
    %
    % Sensor visibility is driven by available data (no checkboxes for now).

    %% CONFIGURATION

    % Date list is derived dynamically from the scanned data; start empty.
    dateList    = datetime.empty(0,1);
    dateStrings = {};

    % CERBERUS filenames, e.g. 2024-11-19_15-11-38_LWIR_Scan_00198_cal_hsi
    CERB_PATTERN      = '*cal_hsi*';
    CERB_TIME_PATTERN = ...
        '(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_(?<hour>\d{2})-(?<min>\d{2})-(?<sec>\d{2})';

    % FAST filenames, e.g. 2024-11-20_1817_17_point-00_LWIR_cal_hsi.hdr
    FAST_PATTERN      = '*cal_hsi.hdr';
    FAST_TIME_PATTERN = ...
        '(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_(?<hour>\d{2})(?<min>\d{2})_(?<sec>\d{2})';

    % FRIDGE header filenames, e.g. AARO_core_7_LWIR.hdr
    FRIDGE_PATTERN = 'AARO*.hdr';
    FRIDGE_DEFAULT_DURATION_SEC = 2;

    perfEnv = getenv('RMBV_PERF');
    perfEnabled = any(strcmp(perfEnv, {'1','true','TRUE','on','ON'}));

    %% STATE STORAGE

    nDays = 0;

    % CERBERUS
    cerbTimesByDay = {};
    cerbMetaByDay  = {};

    % MX20
    mxTimesByDay = {};
    mxMetaByDay  = {};

    % FAST
    fastTimesByDayMap = containers.Map('KeyType','char','ValueType','any');
    fastMetaByDayMap  = containers.Map('KeyType','char','ValueType','any');
    fastModalities    = {};

    % FRIDGE
    fridgeInstancesByDay = {};

    fridgeRootDir          = '';
    hsiRootDir             = '';
    activeWorkspace        = WorkspaceManager('default');
    currentDayIndex        = 0;
    currentFridgeInstances = struct( ...
        'startTime', datetime.empty(0,1), ...
        'endTime',   datetime.empty(0,1), ...
        'wavelength', {{}}, ...
        'path',      {{}} );

    % Rectangle-selection state
    selectionRect = gobjects(1,1);  % handle to selection rectangle patch
    isDragging    = false;
    dragStart     = [NaN NaN];      % [x0 y0] in axes data units
    clickThresh   = 0.02;           % hours; small drag = click

    % Sensor enable flags (tied to checkboxes)
    fridgeEnabled = true;
    hsiCerbEnabled  = true;
    hsiMxEnabled    = true;
    fastEnabledMap  = containers.Map('KeyType','char','ValueType','logical');

    % Sensor availability flags (computed after scanning)
    hasCerbAny   = false;
    hasMxAny     = false;
    hasFridgeAny = false;
    hasFastAny   = false;

    %% UI FIGURE & AXES

    f = uifigure('Name', 'Timeline App', ...
                 'Position', [100 100 900 500]);

    ax = uiaxes('Parent', f, ...
                'Position', [75 170 800 280]);

    % No date yet
    ax.Title.String  = 'Timeline (no data loaded)';
    ax.XLabel.String = 'Time of Day';
    ax.YLabel.String = '';

    ax.YLim   = [0 1];
    ax.YTick  = [];
    ax.YColor = 'none';
    ax.XLim   = [0 24];
    ax.XGrid  = 'on';

    hold(ax, 'on');

    baselineY = 0.5;
    baseLineHandle = plot(ax, [0 24], [baselineY baselineY], '-', 'LineWidth', 2);
    set(baseLineHandle, 'HandleVisibility', 'off');

    % CERBERUS points (blue-ish)
    cerbY = 0.8;
    cerbScatter = scatter(ax, nan, nan, 36, 'filled', ...
        'HitTest', 'off', ...
        'PickableParts', 'none', ...
        'MarkerFaceColor', [0 0.4470 0.7410]);

    % MX20 points (orange-ish, slightly higher)
    mxY = 0.9;
    mxScatter = scatter(ax, nan, nan, 36, 'filled', ...
        'HitTest', 'off', ...
        'PickableParts', 'none', ...
        'MarkerFaceColor', [0.8500 0.3250 0.0980]);

    % FAST points (green-ish, slightly lower than CERBERUS)
    fastY = 0.7;
    fastScatterMap = containers.Map('KeyType','char','ValueType','any');
    fastMarkerCycle = {'^','v','s','d','o','<','>'};
    fastColor       = [0.2 0.6 0.2];

    % Keep a hidden patch handle for FRIDGE so drawFridgeBars can reuse its
    % styling when needed. This handle is also used for the native legend.
    fridgeLegendPatch = patch(ax, [nan nan nan nan], [nan nan nan nan], ...
        [0.20 0.20 0.20], 'EdgeColor', 'none', 'FaceAlpha', 0.85, ...
        'Visible', 'off', 'HandleVisibility', 'on');

    lgd = legend(ax, 'off');

    displayPanel      = uipanel(f, 'Title', 'Display?', 'Position', [700 95 170 160], 'Visible', 'off');
    fridgeCheckbox    = [];
    cerbCheckbox      = [];
    mxCheckbox        = [];
    fastCheckboxMap   = containers.Map('KeyType','char','ValueType','any');

    fridgePatches = gobjects(0);

    % One click handler for the whole axes
    ax.ButtonDownFcn = @axesMouseDown;

    % Listener for XLim changes (zoom / pan) to update ticks
    addlistener(ax, 'XLim', 'PostSet', @(~,~)updateTimeTicks());

    % Enable legacy exploration so pan/zoom objects can coexist with a
    % customized toolbar (avoids pan/set Motion errors on uifigure).
    enableLegacyExplorationModes(f);

    % Attach a toolbar with a restore button that calls the reset helper so
    % the view returns to the full-day window without duplicating ticks.
    tb = axtoolbar(ax, {'pan', 'zoomin', 'zoomout', 'restoreview'});
    restoreBtn = findobj(tb.Children, 'Tooltip', 'Restore View');
    if ~isempty(restoreBtn)
        restoreBtn.ButtonPushedFcn = @(~,~)resetViewLimits();
    end

    % Mouse-move and mouse-up for drag selection
    f.WindowButtonMotionFcn = @mouseMoved;
    f.WindowButtonUpFcn     = @mouseReleased;

    %% DATE DROPDOWN

    uilabel(f, ...
        'Position', [75 105 200 20], ...
        'Text', 'Select Date:', ...
        'HorizontalAlignment', 'left');

    % Start blank & disabled
    dateDropdown = uidropdown(f, ...
        'Position', [75 75 200 30], ...
        'Items', {''}, ...
        'Value', '', ...
        'Enable', 'off', ...
        'ValueChangedFcn', @(~,~)dateChangedCallback());

    %% WORKSPACE / DATA SETUP
    fridgeRootLabel = uilabel(f, ...
        'Position', [300 105 550 20], ...
        'Text', 'Workspace: (none loaded)', ...
        'HorizontalAlignment', 'left');

    hsiRootLabel = uilabel(f, ...
        'Position', [300 45 550 20], ...
        'Text', 'Campaign root: (none selected)', ...
        'HorizontalAlignment', 'left');

    uibutton(f, ...
        'Position', [300 75 150 30], ...
        'Text', 'Data Setup...', ...
        'ButtonPushedFcn', @(~,~)openDataSetupDialog());

    uibutton(f, ...
        'Position', [460 75 150 30], ...
        'Text', 'Refresh Index', ...
        'ButtonPushedFcn', @(~,~)startIndexRefresh());


    % Auto-load last workspace if available.
    try
        lastWsPath = WorkspaceManager('getlast');
        if ~isempty(lastWsPath) && isfile(lastWsPath)
            activeWorkspace = WorkspaceManager('load', lastWsPath);
            fridgeRootLabel.Text = ['Workspace: ' activeWorkspace.name];
            hsiRootLabel.Text = ['Campaign root: ' activeWorkspace.campaignRoot];
            rescanDataAndRefresh();
        end
    catch me
        warning('TimelineApp:WorkspaceLoadFailed', 'Failed to load last workspace: %s', me.message);
    end

    %% PAN & ZOOM
    p = pan(f);   p.Motion = 'horizontal'; p.Enable = 'on';
    z = zoom(f);  z.Motion = 'horizontal'; z.Enable = 'on';

    %======================================================================
    %% DATE CHANGE + ROOT SELECTION
    %======================================================================

    function dateChangedCallback()
        % If dropdown is disabled or blank, do nothing
        if strcmp(dateDropdown.Enable,'off') || isempty(dateDropdown.Value)
            return;
        end
        tPopulate = tic;
        tFirstItems = NaN;

        if isempty(dateList) || isempty(dateStrings)
            ax.Title.String = 'Timeline (no data loaded)';
            return;
        end

        % Switch to selected day and redraw everything
        idx = find(strcmp(dateDropdown.Value, dateStrings));
        if isempty(idx)
            idx = 1;
        end
        if idx > numel(dateList)
            return;
        end
        currentDayIndex = idx;
        currentDate     = dateList(idx);

        ax.Title.String = sprintf('Timeline for %s', datestr(currentDate, 'mm/dd'));

        % Reset X-limits to full day on date change
        ax.XLim = [0 24];

        % Draw a minimal skeleton first so timeline responds quickly.
        cerbScatter.XData = nan; cerbScatter.YData = nan;
        mxScatter.XData = nan;   mxScatter.YData = nan;
        drawnow limitrate nocallbacks;
        if perfEnabled
            fprintf('[perf] timeline initial skeleton: %.1f ms\n', toc(tPopulate)*1000);
        end

        % ----- CERBERUS -----
        timesToday = cerbTimesByDay{idx};
        metaToday  = cerbMetaByDay{idx};

        if isempty(timesToday)
            cerbScatter.XData = nan;
            cerbScatter.YData = nan;
        else
            hoursOfDay = hour(timesToday) + minute(timesToday)/60 + second(timesToday)/3600;
            cerbScatter.XData = hoursOfDay;
            cerbScatter.YData = cerbY * ones(size(hoursOfDay));
        end

        cerbScatter.UserData = struct('times', timesToday, 'meta', metaToday);
        cerbScatter.Visible  = ternary(hsiCerbEnabled,'on','off');
        if isnan(tFirstItems) && ~isempty(timesToday)
            tFirstItems = toc(tPopulate);
            if perfEnabled
                fprintf('[perf] timeline first items (CERBERUS): %.1f ms\n', tFirstItems*1000);
            end
        end

        % ----- MX20 -----
        mxTimesToday = mxTimesByDay{idx};
        mxMetaToday  = mxMetaByDay{idx};

        if isempty(mxTimesToday)
            mxScatter.XData = nan;
            mxScatter.YData = nan;
        else
            mxHoursOfDay = hour(mxTimesToday) + minute(mxTimesToday)/60 + second(mxTimesToday)/3600;
            mxScatter.XData = mxHoursOfDay;
            mxScatter.YData = mxY * ones(size(mxHoursOfDay));
        end

        mxScatter.UserData = struct('times', mxTimesToday, 'meta', mxMetaToday);
        mxScatter.Visible  = ternary(hsiMxEnabled,'on','off');
        if isnan(tFirstItems) && ~isempty(mxTimesToday)
            tFirstItems = toc(tPopulate);
            if perfEnabled
                fprintf('[perf] timeline first items (MX20): %.1f ms\n', tFirstItems*1000);
            end
        end

        % ----- FAST -----
        for fm = 1:numel(fastModalities)
            key = fastModalities{fm};
            sc = ensureFastScatter(key, fm);

            [timesToday, metaToday] = getFastDay(key, idx);

            if isempty(timesToday)
                sc.XData = nan;
                sc.YData = nan;
            else
                fastHours = hour(timesToday) + minute(timesToday)/60 + second(timesToday)/3600;
                sc.XData = fastHours;
                sc.YData = fastY * ones(size(fastHours));
            end

            sc.UserData = struct('times', timesToday, 'meta', metaToday, 'modality', key);
            sc.Visible  = ternary(getOrFastEnabled(key) && hasFastAny, 'on', 'off');
            if isnan(tFirstItems) && ~isempty(timesToday)
                tFirstItems = toc(tPopulate);
                if perfEnabled
                    fprintf('[perf] timeline first items (FAST): %.1f ms\n', tFirstItems*1000);
                end
            end
        end
        if isempty(fastModalities)
            keys = fastScatterMap.keys;
            for fk = 1:numel(keys)
                sc = fastScatterMap(keys{fk});
                if isgraphics(sc)
                    sc.XData = nan;
                    sc.YData = nan;
                    sc.Visible = 'off';
                end
            end
        end

        % ----- FRIDGE -----
        if ~isempty(fridgePatches) && all(isgraphics(fridgePatches))
            delete(fridgePatches);
        end
        fridgePatches = gobjects(0);

        instancesToday         = fridgeInstancesByDay{idx};
        currentFridgeInstances = instancesToday;

        if ~isempty(instancesToday) && ~isempty(instancesToday(1).startTime)
            fridgePatches = drawFridgeBars(ax, instancesToday);
            vis = ternary(fridgeEnabled,'on','off');
            set(fridgePatches, 'Visible', vis);
            if isnan(tFirstItems)
                tFirstItems = toc(tPopulate);
                if perfEnabled
                    fprintf('[perf] timeline first items (FRIDGE): %.1f ms\n', tFirstItems*1000);
                end
            end
        end

        % Initial tick layout for this date
        updateTimeTicks();
        drawnow limitrate nocallbacks;
        if perfEnabled
            fprintf('[perf] timeline full population: %.1f ms\n', toc(tPopulate)*1000);
        end
    end

    function selectFridgeRootCallback()
        newRoot = uigetdir(pwd, 'Select FRIDGE root folder');
        if isequal(newRoot,0)
            return;
        end

        fridgeRootDir = newRoot;
        fridgeRootLabel.Text = ['FRIDGE root: ' fridgeRootDir];

        % Refresh scan results using the latest FRIDGE/HSI roots
        rescanDataAndRefresh();
    end

    function selectHsiRootCallback()
        newRoot = uigetdir(pwd, 'Select HSI root folder');
        if isequal(newRoot,0)
            return;
        end

        hsiRootDir = newRoot;
        hsiRootLabel.Text = ['HSI root: ' hsiRootDir];

        % Refresh scan results using the latest FRIDGE/HSI roots
        rescanDataAndRefresh();
    end

    function rescanDataAndRefresh()
        % Load timeline data from SQLite index for enabled workspace sources.
        if isempty(activeWorkspace) || isempty(activeWorkspace.campaignRoot)
            dateDropdown.Items  = {''};
            dateDropdown.Value  = '';
            dateDropdown.Enable = 'off';
            ax.Title.String     = 'Timeline (no workspace loaded)';
            return;
        end

        loadDlg = uiprogressdlg(f, 'Title', 'Loading timeline', ...
            'Message', 'Reading index summary...', 'Indeterminate', 'on');
        loadCleanup = onCleanup(@() closeProgressDlg(loadDlg)); %#ok<NASGU>

        try
            summary = IndexStore('summaries', activeWorkspace.indexDbPath);
        catch
            summary = table();
        end

        if isempty(summary)
            updateDateList(datetime.empty(0,1));
            resetDataArrays();
            dateDropdown.Items  = {''};
            dateDropdown.Value  = '';
            dateDropdown.Enable = 'off';
            ax.Title.String     = 'Timeline (index empty - run Refresh Index)';
            return;
        end

        enabledIds = {};
        for i = 1:height(summary)
            enabledVal = tableScalar(summary, 'enabled', i);
            if isEnabledValue(enabledVal)
                sidVal = tableScalar(summary, 'source_id', i);
                enabledIds{end+1} = char(string(sidVal)); %#ok<AGROW>
            end
        end
        if isempty(enabledIds)
            return;
        end

        loadDlg.Message = 'Querying indexed timeline items...';
        rows = IndexStore('queryrange', activeWorkspace.indexDbPath, enabledIds, int64(0), intmax('int64'));
        if isempty(rows)
            updateDateList(datetime.empty(0,1));
            resetDataArrays();
            return;
        end

        validTs = [];
        for ri = 1:height(rows)
            tsVal = tableScalar(rows, 'timestamp_utc', ri);
            if ~isempty(tsVal)
                validTs(end+1,1) = double(tsVal); %#ok<AGROW>
            end
        end
        if isempty(validTs)
            updateDateList(datetime.empty(0,1));
            resetDataArrays();
            return;
        end
        dt = datetime(validTs/1000, 'ConvertFrom','posixtime');
        updateDateList(unique(dateshift(dt, 'start', 'day')));
        resetDataArrays();

        loadDlg.Indeterminate = 'off';
        loadDlg.Value = 0;
        for r = 1:height(rows)
            if mod(r, 5000) == 0 || r == height(rows)
                loadDlg.Value = r / max(height(rows),1);
                loadDlg.Message = sprintf('Loading %d / %d indexed items...', r, height(rows));
                drawnow limitrate;
            end
            tsVal = tableScalar(rows, 'timestamp_utc', r);
            if isempty(tsVal), continue; end
            t = datetime(double(tsVal)/1000, 'ConvertFrom','posixtime');
            dayIdx = find(dateList == dateshift(t,'start','day'), 1);
            if isempty(dayIdx), continue; end
            sid = upper(string(tableScalar(rows, 'source_id', r)));
            kind = upper(string(tableScalar(rows, 'kind', r)));
            fpath = char(string(tableScalar(rows, 'file_path', r)));
            if contains(sid, 'FAST')
                key = 'FAST';
                if ~isKey(fastTimesByDayMap,key)
                    fastTimesByDayMap(key) = repmat({datetime.empty(0,1)}, nDays, 1);
                    fastMetaByDayMap(key) = repmat({struct('time', datetime.empty(0,1), 'paths', {{}})}, nDays,1);
                    fastModalities = unique([fastModalities {'FAST'}]);
                end
                timesCell = fastTimesByDayMap(key); metaCell = fastMetaByDayMap(key);
                timesCell{dayIdx}(end+1,1) = t;
                m = metaCell{dayIdx}; m.time(end+1,1)=t; m.paths{end+1,1}=fpath; metaCell{dayIdx}=m;
                fastTimesByDayMap(key)=timesCell; fastMetaByDayMap(key)=metaCell;
            elseif contains(sid,'FRIDGE') || strcmp(kind,'CUBE')
                inst = fridgeInstancesByDay{dayIdx};
                if isempty(inst) || isempty(inst(1).startTime)
                    inst = struct('startTime', datetime.empty(0,1), 'endTime', datetime.empty(0,1), 'wavelength', {{}}, 'path', {{}});
                end
                n = numel(inst) + (~isempty(inst) && ~isempty(inst(1).startTime));
                if isempty(inst) || isempty(inst(1).startTime)
                    n = 1;
                end
                inst(n).startTime = t;
                inst(n).endTime = t + seconds(FRIDGE_DEFAULT_DURATION_SEC);
                inst(n).wavelength = 'UNKNOWN';
                inst(n).path = fpath;
                fridgeInstancesByDay{dayIdx}=inst;
            elseif contains(fpath, 'mx20', 'IgnoreCase', true)
                mxTimesByDay{dayIdx}(end+1,1)=t;
                mm = mxMetaByDay{dayIdx}; mm.time(end+1,1)=t; mm.paths{end+1,1}=fpath; mxMetaByDay{dayIdx}=mm;
            else
                cerbTimesByDay{dayIdx}(end+1,1)=t;
                cm = cerbMetaByDay{dayIdx}; cm.time(end+1,1)=t; cm.paths{end+1,1}=fpath; cerbMetaByDay{dayIdx}=cm;
            end
        end

        hasCerbAny = any(cellfun(@(c) ~isempty(c), cerbTimesByDay(:)'));
        hasMxAny     = any(cellfun(@(c) ~isempty(c), mxTimesByDay(:)'));
        hasFridgeAny = any(cellfun(@(inst) ~isempty(inst) && ~isempty(inst(1).startTime), fridgeInstancesByDay(:)'));
        hasFastAny   = ~isempty(fastModalities);
        if ~isempty(dateStrings)
            dateDropdown.Items  = dateStrings;
            dateDropdown.Value  = dateStrings{1};
            dateDropdown.Enable = 'on';
            updateLegendAndFilters();
            dateChangedCallback();
        end
    end

    function startIndexRefresh()
        if isempty(activeWorkspace) || isempty(activeWorkspace.campaignRoot)
            uialert(f, 'Open Data Setup and select a campaign root first.', 'No Workspace');
            return;
        end
        cancelFlag = false;
        dlg = uiprogressdlg(f, 'Title', 'Indexing', 'Message', 'Preparing index...', 'Cancelable', true, 'Indeterminate', 'off', 'Value', 0);
        cleanup = onCleanup(@() closeProgressDlg(dlg)); %#ok<NASGU>
        progressFcn = @(ev) setProgress(ev);
        cancelFcn = @() cancelRequested();
        updateWorkspaceIndex(activeWorkspace, progressFcn, cancelFcn);
        if cancelFlag
            uialert(f, 'Indexing canceled.', 'Index');
        end
        rescanDataAndRefresh();

        function tf = cancelRequested()
            cancelFlag = dlg.CancelRequested;
            tf = cancelFlag;
        end
        function setProgress(ev)
            dlg.Message = sprintf('Indexing %s\nDirs: %d  Files visited: %d  Indexed: %d\n%s', ev.source, ev.dirsScanned, ev.filesVisited, ev.filesIndexed, ev.currentPath);
            dlg.Value = min(0.95, dlg.Value + 0.001);
            drawnow;
        end
    end

    function openDataSetupDialog()
        d = uifigure('Name','Data Setup','Position',[200 200 900 450]);
        uilabel(d,'Position',[20 410 120 20],'Text','Campaign Root:');
        rootField = uieditfield(d,'text','Position',[140 410 560 22],'Value',activeWorkspace.campaignRoot);
        uibutton(d,'Position',[710 410 120 24],'Text','Change...','ButtonPushedFcn',@(~,~)changeRoot());
        srcTable = uitable(d,'Position',[20 120 860 270],'ColumnName',{'Enabled','Label','Type','Path'},'ColumnEditable',[true false false false]);
        refreshTable();
        uibutton(d,'Position',[20 80 150 28],'Text','Auto-Discover','ButtonPushedFcn',@(~,~)discover());
        uibutton(d,'Position',[180 80 150 28],'Text','Start / Update Index','ButtonPushedFcn',@(~,~)startIndexRefresh());
        uibutton(d,'Position',[340 80 120 28],'Text','Save Workspace','ButtonPushedFcn',@(~,~)saveWs());
        uibutton(d,'Position',[470 80 120 28],'Text','Load Workspace','ButtonPushedFcn',@(~,~)loadWs());
        uibutton(d,'Position',[600 80 120 28],'Text','Close','ButtonPushedFcn',@(~,~)close(d));

        function changeRoot()
            r = uigetdir(pwd, 'Select Campaign Root');
            if isequal(r,0), return; end
            rootField.Value = r;
            activeWorkspace.campaignRoot = r;
            hsiRootLabel.Text = ['Campaign root: ' r];
        end
        function discover()
            activeWorkspace.campaignRoot = rootField.Value;
            dDlg = uiprogressdlg(d, 'Title', 'Discovering sources', ...
                'Message', 'Scanning folders...', 'Cancelable', true, 'Indeterminate', 'on');
            dCleanup = onCleanup(@() closeProgressDlg(dDlg)); %#ok<NASGU>
            prog = @(ev)setDiscoverProgress(ev);
            stop = @() dDlg.CancelRequested;
            activeWorkspace.sources = discoverDataSources(activeWorkspace.campaignRoot, activeWorkspace.excludePatterns, 3, prog, stop);
            refreshTable();

            function setDiscoverProgress(ev)
                dDlg.Message = sprintf('Visited %d folders\n%s', ev.dirsVisited, ev.currentPath);
                drawnow limitrate;
            end
        end
        function refreshTable()
            data = cell(numel(activeWorkspace.sources),4);
            for ii=1:numel(activeWorkspace.sources)
                s = activeWorkspace.sources(ii);
                data{ii,1} = logical(s.enabled);
                data{ii,2} = s.label;
                data{ii,3} = s.type;
                data{ii,4} = s.rootPath;
            end
            srcTable.Data = data;
            srcTable.CellEditCallback = @(~,evt)setEnabled(evt);
        end
        function setEnabled(evt)
            activeWorkspace.sources(evt.Indices(1)).enabled = logical(evt.NewData);
        end
        function saveWs()
            [fn,fp] = uiputfile('*.json', 'Save Workspace', fullfile(prefdir, 'timeline_workspace.json'));
            if isequal(fn,0), return; end
            activeWorkspace.lastOpenedAt = posixtime(datetime('now'));
            WorkspaceManager('save', activeWorkspace, fullfile(fp,fn));
            WorkspaceManager('setlast', fullfile(fp,fn));
            fridgeRootLabel.Text = ['Workspace: ' activeWorkspace.name];
        end
        function loadWs()
            [fn,fp] = uigetfile('*.json', 'Load Workspace');
            if isequal(fn,0), return; end
            activeWorkspace = WorkspaceManager('load', fullfile(fp,fn));
            WorkspaceManager('setlast', fullfile(fp,fn));
            rootField.Value = activeWorkspace.campaignRoot;
            fridgeRootLabel.Text = ['Workspace: ' activeWorkspace.name];
            hsiRootLabel.Text = ['Campaign root: ' activeWorkspace.campaignRoot];
            refreshTable();
            rescanDataAndRefresh();
        end
    end

    function out = tableScalar(tbl, varName, rowIdx)
        % Robust scalar extraction across table variable storage layouts.
        out = [];
        if ~istable(tbl) || rowIdx < 1 || rowIdx > height(tbl) || ~ismember(varName, tbl.Properties.VariableNames)
            return;
        end
        col = tbl.(varName);
        if iscell(col)
            out = col{rowIdx};
        elseif isstring(col) || ischar(col)
            out = col(rowIdx,:);
        else
            out = col(rowIdx);
        end
    end

    function tf = isEnabledValue(v)
        tf = false;
        if isempty(v), return; end
        if islogical(v)
            tf = any(v);
        elseif isnumeric(v)
            tf = any(v ~= 0);
        else
            vv = lower(strtrim(char(string(v))));
            tf = any(strcmp(vv, {'1','true','on','yes'}));
        end
    end

    function updateDateList(newDates)
        % Normalize and store the list of dates used throughout the UI
        newDates = unique(dateshift(newDates, 'start', 'day'));
        dateList    = newDates(:);
        dateStrings = cellstr(datestr(dateList, 'mm/dd'));
        nDays       = numel(dateList);

        if nDays == 0
            currentDayIndex = 0;
        else
            currentDayIndex = 1;
        end
    end

    function resetDataArrays()
        cerbTimesByDay = cell(nDays, 1);
        cerbMetaByDay  = cell(nDays, 1);
        mxTimesByDay   = cell(nDays, 1);
        mxMetaByDay    = cell(nDays, 1);
        fridgeInstancesByDay = cell(nDays, 1);
        fastTimesByDayMap = containers.Map('KeyType','char','ValueType','any');
        fastMetaByDayMap  = containers.Map('KeyType','char','ValueType','any');

        for k = 1:nDays
            cerbTimesByDay{k} = datetime.empty(0,1);
            cerbMetaByDay{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
            mxTimesByDay{k}   = datetime.empty(0,1);
            mxMetaByDay{k}    = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
            fridgeInstancesByDay{k} = struct( ...
                'startTime', datetime.empty(0,1), ...
                'endTime',   datetime.empty(0,1), ...
                'wavelength', {{}}, ...
                'path',      {{}} );
        end

        for fm = 1:numel(fastModalities)
            key = fastModalities{fm};

            timesCell = cell(nDays,1);
            metaCell  = cell(nDays,1);
            for k = 1:nDays
                timesCell{k} = datetime.empty(0,1);
                metaCell{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
            end

            fastTimesByDayMap(key) = timesCell;
            fastMetaByDayMap(key)  = metaCell;
        end
    end

    function sc = ensureFastScatter(modality, idx)
        key = upper(modality);
        if nargin < 2
            idx = 1;
        end

        if ~isKey(fastScatterMap, key) || ~isgraphics(fastScatterMap(key))
            marker = fastMarkerCycle{mod(idx-1, numel(fastMarkerCycle))+1};
            sc = scatter(ax, nan, nan, 60, 'filled', ...
                'HitTest', 'off', 'PickableParts', 'none', ...
                'MarkerFaceColor', fastColor, 'Marker', marker);
            fastScatterMap(key) = sc;
        else
            sc = fastScatterMap(key);
        end
    end

    function val = getOrFastEnabled(modality)
        key = upper(modality);
        if isKey(fastEnabledMap, key)
            val = fastEnabledMap(key);
        else
            val = true;
        end
    end

    function anyData = fastHasAny(modality)
        % True when any day contains FAST data for the modality
        key = upper(modality);
        anyData = false;
        if ~isKey(fastTimesByDayMap, key)
            return;
        end

        cellArr = fastTimesByDayMap(key);
        anyData = any(cellfun(@(c) ~isempty(c), cellArr(:)));
    end

    function [times, meta] = getFastDay(modality, dayIndex)
        % Safely fetch FAST per-day data, tolerating size mismatches
        % between dateList and the stored cell arrays.
        if nargin < 2 || isempty(dayIndex)
            dayIndex = 1;
        end

        times = datetime.empty(0,1);
        meta  = struct('time', datetime.empty(0,1), 'paths', {{}});

        key = upper(modality);
        if ~isKey(fastTimesByDayMap, key) || ~isKey(fastMetaByDayMap, key)
            return;
        end

        timesCell = fastTimesByDayMap(key);
        metaCell  = fastMetaByDayMap(key);

        if dayIndex > numel(timesCell) || dayIndex > numel(metaCell)
            return;
        end

        times = timesCell{dayIndex};
        meta  = metaCell{dayIndex};
    end

    function updateLegendAndFilters()
        legendHandles = [];
        legendNames   = {};

        if hasFridgeAny
            if ~isempty(fridgePatches) && all(isgraphics(fridgePatches))
                set(fridgePatches, 'Visible', 'on');
            end
            fridgeLegendPatch.Visible = 'on';
            set(fridgeLegendPatch, 'HandleVisibility', 'on');
            legendHandles(end+1) = fridgeLegendPatch; %#ok<AGROW>
            legendNames{end+1} = 'FRIDGE'; %#ok<AGROW>
            fridgeEnabled = ensureCheckboxVisible('fridge');
        else
            fridgeEnabled           = false;
            fridgeLegendPatch.Visible = 'off';
            set(fridgeLegendPatch, 'HandleVisibility', 'off');
            hideCheckbox('fridge');
        end

        if hasCerbAny
            cerbScatter.Visible = 'on';
            cerbScatter.HandleVisibility = 'on';
            legendHandles(end+1) = cerbScatter; %#ok<AGROW>
            legendNames{end+1} = 'CERBERUS'; %#ok<AGROW>
            hsiCerbEnabled = ensureCheckboxVisible('cerb');
        else
            hsiCerbEnabled     = false;
            cerbScatter.Visible = 'off';
            cerbScatter.HandleVisibility = 'off';
            hideCheckbox('cerb');
        end

        if hasMxAny
            mxScatter.Visible = 'on';
            mxScatter.HandleVisibility = 'on';
            legendHandles(end+1) = mxScatter; %#ok<AGROW>
            legendNames{end+1} = 'MX-20'; %#ok<AGROW>
            hsiMxEnabled = ensureCheckboxVisible('mx');
        else
            hsiMxEnabled     = false;
            mxScatter.Visible = 'off';
            mxScatter.HandleVisibility = 'off';
            hideCheckbox('mx');
        end

        if hasFastAny && ~isempty(fastModalities)
            for fm = 1:numel(fastModalities)
                key = fastModalities{fm};
                sc = ensureFastScatter(key, fm);
                hasModality = fastHasAny(key);
                if hasModality
                    sc.Visible = 'on';
                    sc.HandleVisibility = 'on';
                    legendHandles(end+1) = sc; %#ok<AGROW>
                    legendNames{end+1} = formatFastLabel(key); %#ok<AGROW>
                    fastEnabledMap(key) = ensureFastCheckbox(key);
                else
                    sc.Visible = 'off';
                    sc.HandleVisibility = 'off';
                    fastEnabledMap(key) = false;
                    hideCheckbox(key);
                end
            end
        elseif ~isempty(fastModalities)
            for fm = 1:numel(fastModalities)
                key = fastModalities{fm};
                if isKey(fastScatterMap, key) && isgraphics(fastScatterMap(key))
                    sc = fastScatterMap(key);
                    sc.Visible = 'off';
                    sc.HandleVisibility = 'off';
                end
                fastEnabledMap(key) = false;
                hideCheckbox(key);
            end
        end

        if isempty(legendHandles)
            legend(ax, 'off');
            displayPanel.Visible = 'off';
            return;
        end

        lgd = legend(ax, legendHandles, legendNames, ...
            'Location', 'southoutside', ...
            'Orientation', 'horizontal'); %#ok<NASGU>
        lgd.AutoUpdate = 'off';

        updateCheckboxLayout();
        applyCheckboxVisibility();
    end

    function closeProgressDlg(dlg)
        if ~isempty(dlg) && isvalid(dlg)
            close(dlg);
        end
    end

    function val = ensureCheckboxVisible(kind)
        ensureDisplayPanelExists();

        switch kind
            case 'fridge'
                if isempty(fridgeCheckbox) || ~isgraphics(fridgeCheckbox)
                    fridgeCheckbox = uicheckbox(displayPanel, ...
                        'Text', 'FRIDGE', ...
                        'Value', true, ...
                        'ValueChangedFcn', @(src,~)onFridgeToggle(src.Value));
                else
                    fridgeCheckbox.Value = true;
                    fridgeCheckbox.Visible = 'on';
                end
                val = logical(fridgeCheckbox.Value);
            case 'cerb'
                if isempty(cerbCheckbox) || ~isgraphics(cerbCheckbox)
                    cerbCheckbox = uicheckbox(displayPanel, ...
                        'Text', 'CERBERUS', ...
                        'Value', true, ...
                        'ValueChangedFcn', @(src,~)onCerbToggle(src.Value));
                else
                    cerbCheckbox.Value = true;
                    cerbCheckbox.Visible = 'on';
                end
                val = logical(cerbCheckbox.Value);
            case 'mx'
                if isempty(mxCheckbox) || ~isgraphics(mxCheckbox)
                    mxCheckbox = uicheckbox(displayPanel, ...
                        'Text', 'MX-20', ...
                        'Value', true, ...
                        'ValueChangedFcn', @(src,~)onMxToggle(src.Value));
                else
                    mxCheckbox.Value = true;
                    mxCheckbox.Visible = 'on';
                end
                val = logical(mxCheckbox.Value);
            otherwise
                val = true;
        end
        displayPanel.Visible = 'on';
    end

    function val = ensureFastCheckbox(modality)
        ensureDisplayPanelExists();
        key = upper(modality);
        if ~isKey(fastCheckboxMap, key) || ~isgraphics(fastCheckboxMap(key))
            cb = uicheckbox(displayPanel, ...
                'Text', formatFastLabel(key), ...
                'Value', true, ...
                'ValueChangedFcn', @(src,~)onFastToggle(key, src.Value));
            fastCheckboxMap(key) = cb;
        else
            cb = fastCheckboxMap(key);
            cb.Value   = true;
            cb.Visible = 'on';
        end
        fastEnabledMap(key) = logical(cb.Value);
        val = fastEnabledMap(key);
        displayPanel.Visible = 'on';
    end

    function hideCheckbox(kind)
        switch kind
            case 'fridge'
                if isgraphics(fridgeCheckbox)
                    fridgeCheckbox.Value   = false;
                    fridgeCheckbox.Visible = 'off';
                end
            case 'cerb'
                if isgraphics(cerbCheckbox)
                    cerbCheckbox.Value   = false;
                    cerbCheckbox.Visible = 'off';
                end
            case 'mx'
                if isgraphics(mxCheckbox)
                    mxCheckbox.Value   = false;
                    mxCheckbox.Visible = 'off';
                end
            otherwise
                key = upper(kind);
                if isKey(fastCheckboxMap, key) && isgraphics(fastCheckboxMap(key))
                    cb = fastCheckboxMap(key);
                    cb.Value   = false;
                    cb.Visible = 'off';
                end
        end
    end

    function updateCheckboxLayout()
        if ~isgraphics(displayPanel)
            return;
        end

        yStart = displayPanel.Position(4) - 35;
        step   = 25;
        nextY  = yStart;

        handles = {};
        if ~isempty(fridgeCheckbox) && isgraphics(fridgeCheckbox) && strcmp(fridgeCheckbox.Visible,'on')
            handles{end+1} = fridgeCheckbox; %#ok<AGROW>
        end
        if ~isempty(cerbCheckbox) && isgraphics(cerbCheckbox) && strcmp(cerbCheckbox.Visible,'on')
            handles{end+1} = cerbCheckbox; %#ok<AGROW>
        end
        if ~isempty(mxCheckbox) && isgraphics(mxCheckbox) && strcmp(mxCheckbox.Visible,'on')
            handles{end+1} = mxCheckbox; %#ok<AGROW>
        end
        if ~isempty(fastCheckboxMap)
            keys = fastCheckboxMap.keys;
            for kk = 1:numel(keys)
                cb = fastCheckboxMap(keys{kk});
                if isgraphics(cb) && strcmp(cb.Visible,'on')
                    handles{end+1} = cb; %#ok<AGROW>
                end
            end
        end

        for hh = 1:numel(handles)
            handles{hh}.Position = [10 nextY 150 20];
            nextY = nextY - step;
        end

        anyVisible = any(cellfun(@(h) isgraphics(h) && strcmp(h.Visible,'on'), handles));

        displayPanel.Visible = ternary(anyVisible, 'on', 'off');
    end

    function label = formatFastLabel(modality)
        key = upper(string(modality));
        if key == "FAST"
            label = 'FAST';
        else
            label = sprintf('FAST %s', key);
        end
    end

    function applyCheckboxVisibility()
        if ~isempty(cerbCheckbox) && isgraphics(cerbCheckbox)
            cerbScatter.Visible = ternary(logical(cerbCheckbox.Value) && hasCerbAny, 'on', 'off');
            hsiCerbEnabled      = strcmp(cerbScatter.Visible, 'on');
        end

        if ~isempty(mxCheckbox) && isgraphics(mxCheckbox)
            mxScatter.Visible = ternary(logical(mxCheckbox.Value) && hasMxAny, 'on', 'off');
            hsiMxEnabled      = strcmp(mxScatter.Visible, 'on');
        end

        if ~isempty(fridgeCheckbox) && isgraphics(fridgeCheckbox) && ~isempty(fridgePatches) && all(isgraphics(fridgePatches))
            vis = ternary(logical(fridgeCheckbox.Value) && hasFridgeAny, 'on', 'off');
            set(fridgePatches, 'Visible', vis);
            fridgeEnabled = strcmp(vis, 'on');
        end

        if ~isempty(fastCheckboxMap)
            keys = fastCheckboxMap.keys;
            for kk = 1:numel(keys)
                cb = fastCheckboxMap(keys{kk});
                if ~isgraphics(cb)
                    continue;
                end
                fastEnabledMap(keys{kk}) = logical(cb.Value) && hasFastAny;
                if isKey(fastScatterMap, keys{kk})
                    sc = fastScatterMap(keys{kk});
                    if isgraphics(sc)
                        visVal = ternary(fastEnabledMap(keys{kk}), 'on', 'off');
                        sc.Visible = visVal;
                        fastScatterMap(keys{kk}) = sc; % ensure stored handle remains graphics
                    end
                end
            end
        end
    end

    function ensureDisplayPanelExists()
        if isempty(displayPanel) || ~isgraphics(displayPanel)
            displayPanel = uipanel(f, 'Title', 'Display?', 'Position', [700 95 170 160]);
        end
    end

    function onFridgeToggle(val)
        fridgeEnabled = logical(val);
        if ~isempty(fridgePatches) && all(isgraphics(fridgePatches))
            set(fridgePatches, 'Visible', ternary(fridgeEnabled && hasFridgeAny, 'on', 'off'));
        end
    end

    function onCerbToggle(val)
        hsiCerbEnabled = logical(val);
        cerbScatter.Visible = ternary(hsiCerbEnabled && hasCerbAny, 'on', 'off');
    end

    function onMxToggle(val)
        hsiMxEnabled = logical(val);
        mxScatter.Visible = ternary(hsiMxEnabled && hasMxAny, 'on', 'off');
    end

    function onFastToggle(modality, val)
        key = upper(modality);
        fastEnabledMap(key) = logical(val);
        if isKey(fastScatterMap, key) && isgraphics(fastScatterMap(key))
            sc = fastScatterMap(key);
            sc.Visible = ternary(fastEnabledMap(key), 'on', 'off');
        end
    end

    %======================================================================
    %% MOUSE INTERACTION
    %======================================================================

    function axesMouseDown(~, event)
        % Start drag-rectangle selection OR single-click selection

        % Convert click to data coordinates
        if isfield(event, 'IntersectionPoint') && ~isempty(event.IntersectionPoint)
            x0 = event.IntersectionPoint(1);
            y0 = event.IntersectionPoint(2);
        else
            cp = ax.CurrentPoint;
            x0 = cp(1,1);
            y0 = cp(1,2);
        end

        % Start drag
        isDragging = true;
        dragStart  = [x0 y0];

        % Remove old rectangle if exists
        if ~isempty(selectionRect) && isgraphics(selectionRect)
            delete(selectionRect);
        end

        % Create a zero-size rectangle to be updated on drag
        selectionRect = patch(ax, ...
            [x0 x0 x0 x0], [y0 y0 y0 y0], ...
            'k', 'FaceAlpha', 0.1, 'EdgeColor', 'k', ...
            'LineStyle', '--', ...
            'HitTest', 'off', 'PickableParts', 'none');
    end

    function mouseMoved(~, ~)
        % Update rectangle while dragging
        if ~isDragging || ~isgraphics(selectionRect)
            return;
        end

        cp = ax.CurrentPoint;
        x1 = cp(1,1);
        y1 = cp(1,2);

        x0 = dragStart(1);
        y0 = dragStart(2);

        xMin = min(x0, x1);
        xMax = max(x0, x1);
        yMin = min(y0, y1);
        yMax = max(y0, y1);

        set(selectionRect, 'XData', [xMin xMax xMax xMin], ...
                           'YData', [yMin yMin yMax yMax]);
    end

    function mouseReleased(~, ~)
        if ~isDragging
            return;
        end
        isDragging = false;

        if ~isgraphics(selectionRect)
            return;
        end

        cp = ax.CurrentPoint;
        x1 = cp(1,1);
        y1 = cp(1,2);

        x0 = dragStart(1);
        y0 = dragStart(2);

        dx = abs(x1 - x0);
        dy = abs(y1 - y0);

        % Save rectangle bounds before deleting the patch
        xMin = min(x0, x1);
        xMax = max(x0, x1);
        yMin = min(y0, y1);
        yMax = max(y0, y1);

        % Remove rectangle visualization
        if isgraphics(selectionRect)
            delete(selectionRect);
        end
        selectionRect = gobjects(1,1);

        % Small movement: treat as click (original behavior)
        if dx < clickThresh && dy < clickThresh
            handlePointClick(x0, y0);
            return;
        end

        % Larger drag: treat as range selection in time
        handleRangeSelection(xMin, xMax, yMin, yMax);
    end

    %======================================================================
    %% POINT CLICK (INFO POPUPS)
    %======================================================================

    function handlePointClick(xClick, yClick)
        % ----- CERBERUS (near cerbY) -----
        cerbUD = cerbScatter.UserData;
        if hsiCerbEnabled && ~isempty(cerbUD) && isfield(cerbUD, 'times') && ~isempty(cerbUD.times)
            timesToday = cerbUD.times;
            metaToday  = cerbUD.meta;

            hoursOfDay = hour(timesToday) + minute(timesToday)/60 + second(timesToday)/3600;

            if abs(yClick - cerbY) < 0.05
                [~, idxPoint] = min(abs(hoursOfDay - xClick));
                if ~isempty(idxPoint) && ~isnan(hoursOfDay(idxPoint))
                    t   = timesToday(idxPoint);
                    pth = metaToday.paths{idxPoint};

                    msg = sprintf(['CERBERUS point\n\nTime: %s\nFile:\n%s'], ...
                        datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF'), pth);
                    uialert(f, msg, 'CERBERUS');
                    return;
                end
            end
        end

        % ----- MX20 (near mxY) -----
        mxUD = mxScatter.UserData;
        if hsiMxEnabled && ~isempty(mxUD) && isfield(mxUD,'times') && ~isempty(mxUD.times)
            mxTimesToday = mxUD.times;
            mxMetaToday  = mxUD.meta;

            mxHours = hour(mxTimesToday) + minute(mxTimesToday)/60 + second(mxTimesToday)/3600;

            if abs(yClick - mxY) < 0.05
                [~, idxPoint] = min(abs(mxHours - xClick));
                if ~isempty(idxPoint) && ~isnan(mxHours(idxPoint))
                    t   = mxTimesToday(idxPoint);
                    pth = mxMetaToday.paths{idxPoint};

                    msg = sprintf(['MX20 point\n\nTime: %s\nFile:\n%s'], ...
                        datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF'), pth);
                    uialert(f, msg, 'MX20');
                    return;
                end
            end
        end

        % ----- FAST (near fastY) -----
        fastKeys = fastScatterMap.keys;
        for fk = 1:numel(fastKeys)
            key = fastKeys{fk};
            sc = fastScatterMap(key);
            if ~getOrFastEnabled(key) || ~isgraphics(sc)
                continue;
            end
            fastUD = sc.UserData;
            if isempty(fastUD) || ~isfield(fastUD,'times') || isempty(fastUD.times)
                continue;
            end
            fastTimes = fastUD.times;
            fastMeta  = fastUD.meta;
            fhours    = hour(fastTimes) + minute(fastTimes)/60 + second(fastTimes)/3600;
            if abs(yClick - fastY) < 0.05
                [~, idxPoint] = min(abs(fhours - xClick));
                if ~isempty(idxPoint) && ~isnan(fhours(idxPoint))
                    t   = fastTimes(idxPoint);
                    pth = fastMeta.paths{idxPoint};
                    msg = sprintf(['FAST %s point\n\nTime: %s\nFile:\n%s'], ...
                        key, datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF'), pth);
                    uialert(f, msg, sprintf('FAST %s', key));
                    return;
                end
            end
        end

        % ----- FRIDGE (bars around y = 0.2) -----
        yBottom = 0.15;
        yTop    = 0.30;
        tolY    = 0.05;

        if ~fridgeEnabled || isempty(currentFridgeInstances) || isempty(currentFridgeInstances(1).startTime)
            return;
        end

        if yClick < (yBottom - tolY) || yClick > (yTop + tolY)
            return;
        end

        nInst = numel(currentFridgeInstances);
        x1    = zeros(nInst,1);
        x2    = zeros(nInst,1);

        for kk = 1:nInst
            t1 = currentFridgeInstances(kk).startTime;
            t2 = currentFridgeInstances(kk).endTime;

            x1(kk) = hour(t1) + minute(t1)/60 + second(t1)/3600;
            x2(kk) = hour(t2) + minute(t2)/60 + second(t2)/3600;
        end

        contains = (xClick >= x1) & (xClick <= x2);

        if any(contains)
            kk = find(contains, 1, 'first');
        else
            centers = (x1 + x2) / 2;
            [~, kk] = min(abs(centers - xClick));
        end

        inst = currentFridgeInstances(kk);

        msg = sprintf(['FRIDGE capture\n\nStart: %s\nEnd:   %s\n\nHeader:\n%s'], ...
            datestr(inst.startTime, 'yyyy-mm-dd HH:MM:SS.FFF'), ...
            datestr(inst.endTime,   'yyyy-mm-dd HH:MM:SS.FFF'), ...
            inst.path);

        uialert(f, msg, 'FRIDGE');
    end

    %======================================================================
    %% RANGE SELECTION (DRAG BOX → VIEWER)
    %======================================================================

    function handleRangeSelection(xMin, xMax, yMin, yMax)
        % Gather HSI selections via helper
        cerbUD = cerbScatter.UserData;
        mxUD   = mxScatter.UserData;
        fastUDMap = containers.Map('KeyType','char','ValueType','any');
        fastKeys = fastScatterMap.keys;
        for fk = 1:numel(fastKeys)
            key = fastKeys{fk};
            if isKey(fastScatterMap, key) && isgraphics(fastScatterMap(key))
                fastUDMap(key) = fastScatterMap(key).UserData;
            end
        end

        [cerbSel, mxSel, fastSel] = selectHSIEventsInRange( ...
            xMin, xMax, yMin, yMax, ...
            cerbUD, mxUD, fastUDMap, cerbY, mxY, fastY, ...
            hsiCerbEnabled, hsiMxEnabled, fastEnabledMap);

        % Geometry-based selection flags so FRIDGE-only drags do not trigger
        % HSI panes when the rectangle never touched an HSI lane.
        hsiThresh = 0.10;
        overlapsCerb = (yMax >= (cerbY - hsiThresh)) && (yMin <= (cerbY + hsiThresh));
        overlapsMx   = (yMax >= (mxY   - hsiThresh)) && (yMin <= (mxY   + hsiThresh));
        overlapsFast = (yMax >= (fastY - hsiThresh)) && (yMin <= (fastY + hsiThresh));

        fridgeBandY = [0.15 0.30];
        overlapsFridge = (yMax >= fridgeBandY(1)) && (yMin <= fridgeBandY(2));

        selection = struct();
        selection.hasHSI    = overlapsCerb || overlapsMx || overlapsFast;
        selection.hasFRIDGE = fridgeEnabled && overlapsFridge;
        selection.tStart    = xMin;
        selection.tEnd      = xMax;

        % Translate the selected hour range into absolute datetimes for the
        % viewer time cursor.
        if currentDayIndex >= 1 && currentDayIndex <= numel(dateList)
            baseDate = dateList(currentDayIndex);
            selection.tStart = baseDate + hours(xMin);
            selection.tEnd   = baseDate + hours(xMax);
        end

        % FRIDGE selection (prefer instance closest to the anchor HSI time)
        anchorTime = [];
        hsiTimes = datetime.empty(0,1);
        if cerbSel.has
            hsiTimes(end+1,1) = cerbSel.time; %#ok<AGROW>
        end
        if mxSel.has
            hsiTimes(end+1,1) = mxSel.time; %#ok<AGROW>
        end
        if fastSel.has
            hsiTimes(end+1,1) = fastSel.time; %#ok<AGROW>
        end
        if ~isempty(hsiTimes)
            anchorTime = min(hsiTimes);
        end

        fridgeSel = struct('has', false, 'instance', []);

        % Only include FRIDGE when enabled and when the selection clearly
        % targets that band (or no HSI was picked, in which case FRIDGE is a
        % reasonable fallback).
        if selection.hasFRIDGE || (fridgeEnabled && ~selection.hasHSI && ...
                ~cerbSel.has && ~mxSel.has)
            fridgeSel = selectFridgeInstanceInRange(xMin, xMax, currentFridgeInstances, anchorTime);
        end

        % Pass all overlapping FRIDGE instances so the viewer can switch
        % captures as the absolute-time slider moves.
        if ~isempty(currentFridgeInstances) && ~isempty(currentFridgeInstances(1).startTime)
            nInst = numel(currentFridgeInstances);
            xs = zeros(nInst,1);
            xe = zeros(nInst,1);
            for kk = 1:nInst
                t1 = currentFridgeInstances(kk).startTime;
                t2 = currentFridgeInstances(kk).endTime;
                xs(kk) = hour(t1) + minute(t1)/60 + second(t1)/3600;
                xe(kk) = hour(t2) + minute(t2)/60 + second(t2)/3600;
            end
            overlapMask = (xe >= xMin) & (xs <= xMax);
            selection.fridgeInstancesInRange = currentFridgeInstances(overlapMask);
        else
            selection.fridgeInstancesInRange = struct([]);
        end

        % Launch viewer via helper
        launchViewerFromSelection(cerbSel, mxSel, fastSel, fridgeSel, xMin, xMax, f, selection);
    end

    %======================================================================
    %% ADAPTIVE TICKS
    %======================================================================

    function updateTimeTicks()
        % Adjust XTick / XTickLabel based on current XLim span

        xl   = ax.XLim;
        span = diff(xl);   % hours

        % Choose step in minutes based on span
        if span >= 12
            stepMin = 60;   % 1 hour
        elseif span >= 4
            stepMin = 30;
        elseif span >= 1
            stepMin = 10;
        elseif span >= 0.5
            stepMin = 5;
        else
            stepMin = 1;
        end

        stepH = stepMin / 60;

        % Align ticks nicely to the grid
        startH = ceil(xl(1) / stepH) * stepH;
        endH   = floor(xl(2) / stepH) * stepH;

        if endH < startH
            ticks = mean(xl);
        else
            ticks = startH:stepH:endH;
        end

        ax.XTick      = ticks;
        ax.XTickLabel = arrayfun(@fmtHourMinute, ticks, 'UniformOutput', false);
    end

    function resetViewLimits()
        ax.XLim      = [0 24];
        ax.XLimMode  = 'manual';
        ax.XTickMode = 'manual';
        updateTimeTicks();
    end

    function s = fmtHourMinute(x)
        % Format a numeric hour value as 'HH:MM'
        h = floor(x);
        m = round((x - h)*60);
        if m == 60
            h = h + 1;
            m = 0;
        end
        s = sprintf('%02d:%02d', h, m);
    end

    % Small utility (ternary for strings)
    function out = ternary(cond, a, b)
        if cond
            out = a;
        else
            out = b;
        end
    end

end
