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
    % Checkboxes below the plot control which HSI sensors are
    % visible and selectable.

    %% CONFIGURATION

    % Date list is derived dynamically from the scanned data; start empty.
    dateList    = datetime.empty(0,1);
    dateStrings = {};

    % CERBERUS filenames, e.g. 2024-11-19_15-11-38_LWIR_Scan_00198_cal_hsi
    CERB_PATTERN      = '*cal_hsi*';
    CERB_TIME_PATTERN = ...
        '(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_(?<hour>\d{2})-(?<min>\d{2})-(?<sec>\d{2})';

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

    % Sensor availability flags (computed after scanning)
    hasCerbAny   = false;
    hasMxAny     = false;
    hasFridgeAny = false;

    %% UI FIGURE & AXES

    f = uifigure('Name', 'Timeline App', ...
                 'Position', [100 100 900 500]);

    ax = uiaxes('Parent', f, ...
                'Position', [75 150 800 320]);

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
    plot(ax, [0 24], [baselineY baselineY], '-', 'LineWidth', 2);

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

    % Legend handles (populated dynamically based on available data)
    fridgeLegendPatch = patch(ax, [nan nan nan nan], [nan nan nan nan], ...
        [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
        'Visible', 'off', 'HandleVisibility', 'on');

    lgd = legend(ax, 'off');

    % Sensor enable/disable checkboxes (placed near the legend) – hidden
    % until data is loaded so sensors do not appear before roots are
    % selected.
    displayLabel = uilabel(f, ...
        'Text', 'Display?', ...
        'Position', [650 115 80 20], ...
        'HorizontalAlignment', 'center'); %#ok<NASGU>

    cbFridge = uicheckbox(f, ...
        'Text', 'FRIDGE', ...
        'Value', false, ...
        'Position', [650 95 80 20], ...
        'Visible', 'off', ...
        'Enable', 'off', ...
        'ValueChangedFcn', @(cb,~)toggleFridge(cb.Value)); %#ok<NASGU>

    cbCerb = uicheckbox(f, ...
        'Text', 'CERBERUS', ...
        'Value', false, ...
        'Position', [650 75 90 20], ...
        'Visible', 'off', ...
        'Enable', 'off', ...
        'ValueChangedFcn', @(cb,~)toggleCerb(cb.Value)); %#ok<NASGU>

    cbMx = uicheckbox(f, ...
        'Text', 'MX20', ...
        'Value', false, ...
        'Position', [650 55 80 20], ...
        'Visible', 'off', ...
        'Enable', 'off', ...
        'ValueChangedFcn', @(cb,~)toggleMx(cb.Value)); %#ok<NASGU>

    fridgePatches = gobjects(0);

    % One click handler for the whole axes
    ax.ButtonDownFcn = @axesMouseDown;

    % Listener for XLim changes (zoom / pan) to update ticks
    addlistener(ax, 'XLim', 'PostSet', @(~,~)updateTimeTicks());

    % Enable legacy exploration so pan/zoom objects can coexist with a
    % customized toolbar (avoids pan/set Motion errors on uifigure).
    enableLegacyExplorationModes(ax);

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
    %% CHECKBOX CALLBACKS
    %======================================================================

    function toggleFridge(val)
        if ~hasFridgeAny
            cbFridge.Value = false;
            fridgeEnabled  = false;
            set(fridgePatches, 'Visible', 'off');
            return;
        end
        fridgeEnabled = logical(val);
        vis = ternary(fridgeEnabled, 'on', 'off');
        if ~isempty(fridgePatches) && all(isgraphics(fridgePatches))
            set(fridgePatches, 'Visible', vis);
        end
        updateLegendAndFilters();
    end

    function toggleCerb(val)
        if ~hasCerbAny
            cbCerb.Value = false;
            cerbScatter.Visible = 'off';
            return;
        end
        hsiCerbEnabled = logical(val);
        if hsiCerbEnabled
            cerbScatter.Visible = 'on';
        else
            cerbScatter.Visible = 'off';
        end
    end

    function toggleMx(val)
        if ~hasMxAny
            cbMx.Value = false;
            mxScatter.Visible = 'off';
            return;
        end
        hsiMxEnabled = logical(val);
        if hsiMxEnabled
            mxScatter.Visible = 'on';
        else
            mxScatter.Visible = 'off';
        end
    end

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

        % --- CERBERUS HSI ---
        cerbRoot = '';
        if ~isempty(hsiRootDir)
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
        else
            fprintf('HSI root not set; skipping CERBERUS/MX20 scanning.\n');
        end

        % --- FRIDGE ---
        if ~isempty(fridgeRootDir)
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
        if ~isempty(cerbRoot)
            [cerbTimesByDay, cerbMetaByDay] = scanCerberusFiles( ...
                cerbRoot, dateList, CERB_PATTERN, CERB_TIME_PATTERN);

            fprintf('\nCERBERUS event counts per day:\n');
            for di = 1:numel(dateList)
                fprintf('  %s: %d events\n', datestr(dateList(di), 'mm/dd'), ...
                        numel(cerbTimesByDay{di}));
            end
        end

        if ~isempty(mxRoot)
            [mxTimesByDay, mxMetaByDay] = scanMX20Files( ...
                mxRoot, dateList, CERB_TIME_PATTERN);

            fprintf('\nMX20 event counts per day:\n');
            for di = 1:numel(dateList)
                fprintf('  %s: %d events\n', datestr(dateList(di), 'mm/dd'), ...
                        numel(mxTimesByDay{di}));
            end
        end

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

            % If you ONLY want days with HSI images, use:
            %   hasData(di) = hasCerb || hasMx;
            % Currently: any CERB, MX20, or FRIDGE counts as "has data"
            hasData(di) = hasCerb || hasMx || hasFridge;
        end

        hasCerbAny   = any(cellfun(@(c) ~isempty(c), cerbTimesByDay(:)'));
        hasMxAny     = any(cellfun(@(c) ~isempty(c), mxTimesByDay(:)'));
        hasFridgeAny = any(cellfun(@(inst) ~isempty(inst) && ~isempty(inst(1).startTime), fridgeInstancesByDay(:)'));

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
                'No CERBERUS, MX20, or FRIDGE data found for any configured dates.', ...
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
    end

    function updateLegendAndFilters()
        handles = [];
        labels  = {};

        if hasFridgeAny
            if strcmp(cbFridge.Enable, 'off')
                fridgeEnabled = true; % default to checked when data appears
            end
            cbFridge.Visible = 'on';
            cbFridge.Enable  = 'on';
            cbFridge.Value   = fridgeEnabled;
            if ~isempty(fridgePatches) && all(isgraphics(fridgePatches))
                set(fridgePatches, 'Visible', ternary(fridgeEnabled,'on','off'));
            end
            fridgeLegendPatch.Visible = 'on';
            handles(end+1) = fridgeLegendPatch; %#ok<AGROW>
            labels{end+1}  = 'FRIDGE'; %#ok<AGROW>
        else
            fridgeEnabled      = false;
            cbFridge.Value     = false;
            cbFridge.Enable    = 'off';
            cbFridge.Visible   = 'off';
            fridgeLegendPatch.Visible = 'off';
        end

        if hasCerbAny
            if strcmp(cbCerb.Enable, 'off')
                hsiCerbEnabled = true;  % default to checked when data appears
            end
            cbCerb.Visible = 'on';
            cbCerb.Enable  = 'on';
            cbCerb.Value   = hsiCerbEnabled;
            cerbScatter.Visible = ternary(hsiCerbEnabled,'on','off');
            handles(end+1) = cerbScatter; %#ok<AGROW>
            labels{end+1}  = 'CERBERUS'; %#ok<AGROW>
        else
            hsiCerbEnabled   = false;
            cbCerb.Value     = false;
            cbCerb.Enable    = 'off';
            cbCerb.Visible   = 'off';
            cerbScatter.Visible = 'off';
        end

        if hasMxAny
            if strcmp(cbMx.Enable, 'off')
                hsiMxEnabled = true;  % default to checked when data appears
            end
            cbMx.Visible = 'on';
            cbMx.Enable  = 'on';
            cbMx.Value   = hsiMxEnabled;
            mxScatter.Visible = ternary(hsiMxEnabled,'on','off');
            handles(end+1) = mxScatter; %#ok<AGROW>
            labels{end+1}  = 'MX20'; %#ok<AGROW>
        else
            hsiMxEnabled   = false;
            cbMx.Value     = false;
            cbMx.Enable    = 'off';
            cbMx.Visible   = 'off';
            mxScatter.Visible = 'off';
        end

        if isempty(handles)
            lgd = legend(ax, 'off'); %#ok<NASGU>
        else
            lgd = legend(ax, handles, labels, 'Location', 'southoutside'); %#ok<NASGU>
            lgd.AutoUpdate = 'off';
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

        [cerbSel, mxSel] = selectHSIEventsInRange( ...
            xMin, xMax, yMin, yMax, ...
            cerbUD, mxUD, cerbY, mxY, ...
            hsiCerbEnabled, hsiMxEnabled);

        % FRIDGE selection (prefer instance closest to the anchor HSI time)
        anchorTime = [];
        if cerbSel.has && mxSel.has
            anchorTime = min(cerbSel.time, mxSel.time); % earliest HSI in range
        elseif cerbSel.has
            anchorTime = cerbSel.time;
        elseif mxSel.has
            anchorTime = mxSel.time;
        end

        fridgeSel = struct('has', false, 'instance', []);

        % Only include FRIDGE when enabled and when the selection clearly
        % targets that band (or no HSI was picked, in which case FRIDGE is a
        % reasonable fallback).
        fridgeBandY = [0.15 0.30];
        overlapsFridge = (yMax >= fridgeBandY(1)) && (yMin <= fridgeBandY(2));

        if fridgeEnabled && (overlapsFridge || (~cerbSel.has && ~mxSel.has))
            fridgeSel = selectFridgeInstanceInRange(xMin, xMax, currentFridgeInstances, anchorTime);
        end

        % Launch viewer via helper
        launchViewerFromSelection(cerbSel, mxSel, fridgeSel, xMin, xMax, f);
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
