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
        [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
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

    %% DATA ROOT BUTTONS
    % Allow independent selection of FRIDGE and HSI roots so the user can
    % point FRIDGE scanning and HSI scanning at different locations.

    fridgeRootLabel = uilabel(f, ...
        'Position', [300 105 550 20], ...
        'Text', 'FRIDGE root: (none selected)', ...
        'HorizontalAlignment', 'left');

    uibutton(f, ...
        'Position', [300 75 150 30], ...
        'Text', 'Select FRIDGE Directory', ...
        'ButtonPushedFcn', @(~,~)selectFridgeRootCallback());

    hsiRootLabel = uilabel(f, ...
        'Position', [300 45 550 20], ...
        'Text', 'HSI root: (none selected)', ...
        'HorizontalAlignment', 'left');

    uibutton(f, ...
        'Position', [300 15 150 30], ...
        'Text', 'Select HSI Directory', ...
        'ButtonPushedFcn', @(~,~)selectHsiRootCallback());

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
        end

        % Initial tick layout for this date
        updateTimeTicks();
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
        % Re-scan FRIDGE and HSI data using independent roots. This clears
        % stale state so switching directories cannot leave old events on
        % the plot.

        dlg = uiprogressdlg(f, 'Title','Scanning', ...
            'Message','Scanning directories...', ...
            'Indeterminate','on', ...
            'Cancelable', true);
        dlgCleanup = onCleanup(@() closeProgressDlg(dlg)); %#ok<NASGU>

        if isempty(fridgeRootDir) && isempty(hsiRootDir)
            updateDateList(datetime.empty(0,1));
            resetDataArrays();
            dateDropdown.Items  = {''};
            dateDropdown.Value  = '';
            dateDropdown.Enable = 'off';
            ax.Title.String     = 'Timeline (no data loaded)';
            updateLegendAndFilters();
            uialert(f, ...
                'Select at least one FRIDGE or HSI directory to scan for data.', ...
                'No Data Roots');
            return;
        end

        % First pass: scan each sensor independently to learn which dates
        % exist. This allows the dropdown to adapt to whatever files are
        % present without any hard-coded date ranges.
        dateCandidates = datetime.empty(0,1);
        mxRoot = '';
        fastRoot = '';

        % --- CERBERUS HSI ---
        cerbRoot = '';
        if ~isempty(hsiRootDir)
            if dlg.CancelRequested, return; end
            dlg.Message = 'Scanning CERBERUS...'; drawnow;
            % Prefer nested layout HSI/CERBERUS, then CERBERUS/, then root.
            cerbCandidates = { ...
                fullfile(hsiRootDir, 'HSI', 'CERBERUS'), ...
                fullfile(hsiRootDir, 'CERBERUS'), ...
                hsiRootDir};

            for cc = 1:numel(cerbCandidates)
                if isfolder(cerbCandidates{cc})
                    cerbRoot = cerbCandidates{cc};
                    break;
                end
            end

            if isempty(cerbRoot)
                fprintf('No CERBERUS folder found under configured HSI root (%s).\n', hsiRootDir);
            else
                [~, ~, cerbDates] = scanCerberusFiles( ...
                    cerbRoot, datetime.empty(0,1), CERB_PATTERN, CERB_TIME_PATTERN);
                dateCandidates = [dateCandidates; cerbDates(:)]; %#ok<AGROW>
            end

            % --- MX20 HSI ---
            % Prefer MX20 subfolders so CERBERUS headers do not get
            % misinterpreted as MX20. Only fall back to the HSI root if we
            % explicitly find MX20 in the path.
            if dlg.CancelRequested, return; end
            dlg.Message = 'Scanning MX20...'; drawnow;
            mxCandidates = { ...
                fullfile(hsiRootDir, 'HSI', 'MX20'), ...
                fullfile(hsiRootDir, 'MX20')};

            mxRoot = '';
            for mc = 1:numel(mxCandidates)
                if isfolder(mxCandidates{mc})
                    mxRoot = mxCandidates{mc};
                    break;
                end
            end

            if ~isempty(mxRoot)
                [~, ~, mxDates] = scanMX20Files( ...
                    mxRoot, datetime.empty(0,1), CERB_TIME_PATTERN);
                dateCandidates = [dateCandidates; mxDates(:)]; %#ok<AGROW>
            else
                fprintf('No MX20 folder found under configured HSI root (%s).\n', hsiRootDir);
            end

            % --- FAST HSI ---
            if dlg.CancelRequested, return; end
            dlg.Message = 'Scanning FAST...'; drawnow;
            fastCandidates = { ...
                fullfile(hsiRootDir, 'HSI', 'FAST'), ...
                fullfile(hsiRootDir, 'FAST')};

            for fc = 1:numel(fastCandidates)
                if isfolder(fastCandidates{fc})
                    fastRoot = fastCandidates{fc};
                    break;
                end
            end

            if ~isempty(fastRoot)
                [~, ~, fastDates, fastMods] = scanFastFiles( ...
                    fastRoot, datetime.empty(0,1));
                dateCandidates = [dateCandidates; fastDates(:)]; %#ok<AGROW>
                fastModalities = fastMods;
            else
                fprintf('No FAST folder found under configured HSI root (%s).\n', hsiRootDir);
                fastModalities = {};
            end
        else
            fprintf('HSI root not set; skipping CERBERUS/MX20 scanning.\n');
        end

        % --- FRIDGE ---
        if ~isempty(fridgeRootDir)
            if dlg.CancelRequested, return; end
            dlg.Message = 'Scanning FRIDGE headers...'; drawnow;
            [~, fridgeDates] = scanFridgeHeaders( ...
                fridgeRootDir, datetime.empty(0,1), FRIDGE_PATTERN, FRIDGE_DEFAULT_DURATION_SEC);
            dateCandidates = [dateCandidates; fridgeDates(:)]; %#ok<AGROW>
        else
            fprintf('FRIDGE root not set; skipping FRIDGE scanning.\n');
        end

        dateCandidates = unique(dateshift(dateCandidates,'start','day'));
        updateDateList(dateCandidates);
        resetDataArrays();

        % Second pass: with the unified date list, populate per-day buckets
        % for each sensor.
        if dlg.CancelRequested, return; end
        dlg.Message = 'Loading CERBERUS events...'; drawnow;
        if ~isempty(cerbRoot)
            [cerbTimesByDay, cerbMetaByDay] = scanCerberusFiles( ...
                cerbRoot, dateList, CERB_PATTERN, CERB_TIME_PATTERN);

            fprintf('\nCERBERUS event counts per day:\n');
            for di = 1:numel(dateList)
                fprintf('  %s: %d events\n', datestr(dateList(di), 'mm/dd'), ...
                        numel(cerbTimesByDay{di}));
            end
        end

        if dlg.CancelRequested, return; end
        dlg.Message = 'Loading MX20 events...'; drawnow;
        if ~isempty(mxRoot)
            [mxTimesByDay, mxMetaByDay] = scanMX20Files( ...
                mxRoot, dateList, CERB_TIME_PATTERN);

            fprintf('\nMX20 event counts per day:\n');
            for di = 1:numel(dateList)
                fprintf('  %s: %d events\n', datestr(dateList(di), 'mm/dd'), ...
                        numel(mxTimesByDay{di}));
            end
        end

        if dlg.CancelRequested, return; end
        dlg.Message = 'Loading FAST events...'; drawnow;
        if ~isempty(fastRoot)
            [fastTimesByDayMap, fastMetaByDayMap, ~, fastMods] = scanFastFiles( ...
                fastRoot, dateList);
            if ~isempty(fastMods)
                fastModalities = fastMods;
            end

            fprintf('\nFAST event counts per day:\n');
            for di = 1:numel(dateList)
                counts = arrayfun(@(m) numel(getFastDay(m{1}, di)), fastModalities);
                fprintf('  %s: %d events\n', datestr(dateList(di), 'mm/dd'), sum(counts));
            end
        end

        if dlg.CancelRequested, return; end
        dlg.Message = 'Loading FRIDGE instances...'; drawnow;
        if ~isempty(fridgeRootDir)
            [fridgeInstancesByDay, ~] = scanFridgeHeaders( ...
                fridgeRootDir, dateList, FRIDGE_PATTERN, FRIDGE_DEFAULT_DURATION_SEC);

            fprintf('\nFRIDGE instance counts per day:\n');
            for di = 1:numel(dateList)
                nInst = 0;
                insts = fridgeInstancesByDay{di};
                if ~isempty(insts) && ~isempty(insts(1).startTime)
                    nInst = numel(insts);
                end
                fprintf('  %s: %d instances\n', datestr(dateList(di), 'mm/dd'), nInst);
            end
        end

        % ---------- filter dropdown to days that have any data ----------
        nDaysLocal = numel(dateList);
        hasData = false(nDaysLocal,1);

        for di = 1:nDaysLocal
            hasCerb   = ~isempty(cerbTimesByDay{di});
            hasMx     = ~isempty(mxTimesByDay{di});
            insts     = fridgeInstancesByDay{di};
            hasFridge = ~isempty(insts) && ~isempty(insts(1).startTime);
            hasFast   = false;
            for fm = 1:numel(fastModalities)
                key = fastModalities{fm};
                if ~isempty(getFastDay(key, di))
                    hasFast = true;
                    break;
                end
            end

            % If you ONLY want days with HSI images, use:
            %   hasData(di) = hasCerb || hasMx;
            % Currently: any CERB, MX20, or FRIDGE counts as "has data"
            hasData(di) = hasCerb || hasMx || hasFridge || hasFast;
        end

        hasCerbAny   = any(cellfun(@(c) ~isempty(c), cerbTimesByDay(:)'));
        hasMxAny     = any(cellfun(@(c) ~isempty(c), mxTimesByDay(:)'));
        hasFridgeAny = any(cellfun(@(inst) ~isempty(inst) && ~isempty(inst(1).startTime), fridgeInstancesByDay(:)'));
        hasFastAny   = false;
        for fm = 1:numel(fastModalities)
            key = fastModalities{fm};
            if fastHasAny(key)
                hasFastAny = true; %#ok<AGROW>
            end
        end

        validIdx = find(hasData);

        if isempty(validIdx)
            % No data for any configured dates
            dateDropdown.Items  = {''};
            dateDropdown.Value  = '';
            dateDropdown.Enable = 'off';
            ax.Title.String     = 'Timeline (no data found)';
            hasCerbAny   = false;
            hasMxAny     = false;
            hasFridgeAny = false;
            updateLegendAndFilters();
            uialert(f, ...
                'No CERBERUS, MX20, FAST, or FRIDGE data found for any configured dates.', ...
                'No Data');
            return;  % important: do NOT call dateChangedCallback()
        else
            % Restrict dropdown items to only dates that have data
            newItems = dateStrings(validIdx);
            dateDropdown.Items  = newItems;
            dateDropdown.Value  = newItems{1};
            dateDropdown.Enable = 'on';

            % Update title to first valid date
            ax.Title.String = sprintf('Timeline for %s', newItems{1});
            updateLegendAndFilters();
        end
        % -----------------------------------------------------------------

        % Redraw for the (possibly new) current date
        dateChangedCallback();
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
            legendNames{end+1} = 'MX20'; %#ok<AGROW>
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
                    legendNames{end+1} = sprintf('FAST %s', key); %#ok<AGROW>
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
                        'Text', 'MX20', ...
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
                'Text', sprintf('FAST %s HSI', key), ...
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

        yStart = 50;
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
            displayPanel = uipanel(f, 'Title', 'Display?', 'Position', [700 95 170 100]);
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
