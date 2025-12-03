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
    page.RowHeight   = {34, '1x', 'fit'};
    page.ColumnWidth = {'1x','1x','1x'};

    % Header
    header = uilabel(page, ...
        'Text','Multiband FRIDGE + HSI viewer (driven by timeline selection).', ...
        'FontWeight','bold','HorizontalAlignment','center');
    header.Layout.Row    = 1;
    header.Layout.Column = [1 3];

    % Axes names and grid positions
    modalities = {'LWIR','MWIR','SWIR','MONO','VIS-COLOR'};
    axPos      = [1,1; 1,2; 1,3; 2,1; 2,2];
    axMap      = containers.Map;

    % Image grid (2x3)
    imgGrid = uigridlayout(page,[2,3]);
    imgGrid.Layout.Row    = 2;
    imgGrid.Layout.Column = [1 3];
    imgGrid.RowHeight     = {'1x','1x'};
    imgGrid.ColumnWidth   = {'1x','1x','1x'};

    %----------------------------------------------------------------------
    % FRIDGE panes: panel + inner grid (label row + axes row)
    %----------------------------------------------------------------------
    for i = 1:numel(modalities)
        pnl = uipanel(imgGrid);
        pnl.Layout.Row    = axPos(i,1);
        pnl.Layout.Column = axPos(i,2);

        pGrid = uigridlayout(pnl,[2,1]);
        pGrid.RowHeight   = {'fit','1x'};
        pGrid.ColumnWidth = {'1x'};

        modName = modalities{i};
        if strcmp(modName,'VIS-COLOR')
            dispName = 'FRIDGE VIS';
        else
            dispName = ['FRIDGE ' modName];
        end

        lblTop = uilabel(pGrid, ...
            'Text', dispName, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center');
        lblTop.Layout.Row    = 1;
        lblTop.Layout.Column = 1;

        ax = uiaxes(pGrid);
        ax.Layout.Row    = 2;
        ax.Layout.Column = 1;
        axis(ax,'off');
        title(ax, modName);

        axMap(modalities{i}) = ax;
    end

    % HSI tab group (CERB LWIR/VNIR + MX20)
    cerbTabs = uitabgroup(imgGrid);
    cerbTabs.Layout.Row    = 2;
    cerbTabs.Layout.Column = 3;
    
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

    fileRow = uigridlayout(infoCol,[1,3]);
    fileRow.ColumnWidth = {'1x','fit','fit'};
    lblFile   = uilabel(fileRow,'Text','Filename: -','HorizontalAlignment','left');
    lblFrames = uilabel(fileRow,'Text','Frames: -');
    lblMem    = uilabel(fileRow,'Text','');

    pixelRow = uigridlayout(infoCol,[1,2]);
    pixelRow.ColumnWidth = {'1x','1x'};
    lblPixel = uilabel(pixelRow, ...
        'Text','Pixel: -', ...
        'HorizontalAlignment','left');
    lblValue = uilabel(pixelRow, ...
        'Text','Value: -', ...
        'HorizontalAlignment','left');

    % Center column: navigation + slider kept wide
    navCol = uigridlayout(ctrlWrapper,[2,1]);
    navCol.RowHeight   = {'fit','fit'};
    navCol.ColumnWidth = {'1x'};

    navTop = uigridlayout(navCol,[1,6]);
    navTop.ColumnWidth = {'fit','fit','fit','fit','fit','fit'};
    uilabel(navTop,'Text','Out-of-range:','HorizontalAlignment','right');
    ddBehavior = uidropdown(navTop,'Items',{'Hold last','Show missing','Loop'}, ...
        'Value','Hold last', ...
        'Tooltip','When a modality runs out of frames', ...
        'ValueChangedFcn',@(dd,~) setBehavior(dd.Value)); %#ok<NASGU>
    uilabel(navTop,'Text','Go to frame:','HorizontalAlignment','right');
    goField = uieditfield(navTop,'numeric','Limits',[1 Inf], ...
        'RoundFractionalValues','on','Value',1,'Enable','off');
    btnGo   = uibutton(navTop,'Text','Go','Enable','off', ...
        'ButtonPushedFcn',@(~,~)gotoFrame());
    btnSave = uibutton(navTop,'Text','Save Montage PNG','Enable','off', ...
        'ButtonPushedFcn',@(~,~)saveMontage());

    navBottom = uigridlayout(navCol,[1,6]);
    navBottom.ColumnWidth = {'fit','fit','fit','1x','fit','fit'};
    btnPrev = uibutton(navBottom,'Text','Previous','Enable','off', ...
        'ButtonPushedFcn',@(~,~)step(-1));
    btnNext = uibutton(navBottom,'Text','Next','Enable','off', ...
        'ButtonPushedFcn',@(~,~)step(+1));
    uilabel(navBottom,'Text','Time slider:','HorizontalAlignment','right');
    frameSlider = uislider(navBottom, ...
        'Limits',[1 2], ...
        'Value',1, ...
        'MajorTicks',[], ...
        'MinorTicks',[], ...
        'Enable','off', ...
        'ValueChangedFcn',@frameSliderChanged);
    % Spacer to keep right-side items from crowding the slider
    uilabel(navBottom,'Text','');
    % Timestamp label sits in the right column

    % Right column: dedicated timestamp display
    timeCol = uigridlayout(ctrlWrapper,[2,1]);
    timeCol.RowHeight   = {'fit','fit'};
    timeCol.ColumnWidth = {'1x'};
    uilabel(timeCol,'Text','Current time','FontWeight','bold', ...
        'HorizontalAlignment','left');
    lblTime = uilabel(timeCol,'Text','Time: -','HorizontalAlignment','left');

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
    S.behavior     = 'hold';   % 'hold'|'missing'|'loop'
    S.fridgeTimes  = [];       % legacy datetime vector
    S.fridgeTimesMap = containers.Map(modalities, repmat({datetime.empty(0,1)},1,numel(modalities)));
    S.timelineTimes  = datetime.empty(0,1);   % union of all modality times
    S.sliderMode   = 'frame';  % 'frame' (fallback) or 'time'
    S.sliderOrigin = NaT;      % reference time for slider (time mode)
    S.hsiEvents    = struct('sensor', {}, 'time', {}, 'path', {});
    S.currentHsi   = struct('sensor','', 'time', NaT);

    S.cerb = struct('LWIR',[],'VNIR',[]);
    S.mx20 = struct('hdr',[],'ctx',[]);

    targetStartTime = [];
    if isfield(initial,'initialTime')
        targetStartTime = initial.initialTime;
    end
    if isfield(initial,'hsiEvents')
        S.hsiEvents = initial.hsiEvents;
    end
    sliderInternalUpdate = false;  % prevent recursive slider callbacks

    %======================== CORE CALLBACKS ===============================
    function loadFromRawFile(fullRawPath)
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
        [S.files, S.hdrs, S.exists, S.maxFrames, ...
         S.nFrames, S.fridgeTimes, S.fridgeTimesMap] = ...
            fridge_init_from_raw(path, prefix, modalities);
        S.frameCount = S.nFrames;

        % If the timeline already passed per-frame timestamps, prefer them so
        % alignment works even when headers omit band_names.
        if isfield(initial, 'fridgeTimes') && ~isempty(initial.fridgeTimes)
            S.fridgeTimes = initial.fridgeTimes(:);
            for ii = 1:numel(modalities)
                S.fridgeTimesMap(modalities{ii}) = S.fridgeTimes;
            end
        end
        if isempty(S.fridgeTimes) && ~isempty(S.fridgeTimesMap)
            % Legacy callers may have only populated the map
            mods = S.fridgeTimesMap.keys;
            for ii = 1:numel(mods)
                if ~isempty(S.fridgeTimesMap(mods{ii}))
                    S.fridgeTimes = S.fridgeTimesMap(mods{ii});
                    break;
                end
            end
        end
        % -----------------------------------------------------------------

        lblFile.Text   = ['Filename: ', [file ext]];
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
        goField.Enable     = 'on';
        btnGo.Enable       = 'on';

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
        goField.Value = S.frame;
        setSliderFromFrame();
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
    end

    function gotoFrame()
        val = round(goField.Value);
        if ~isfinite(val) || val < 1
            val = 1;
        elseif val > S.nFrames
            val = S.nFrames;
        end
        S.frame = val;
        goField.Value = S.frame;
        setSliderFromFrame();
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
    end

    function frameSliderChanged(src, ~)
        if sliderInternalUpdate
            return;
        end

        if strcmp(S.sliderMode,'time') && hasFridgeTimes()
            targetTime = S.sliderOrigin + seconds(src.Value);
            jumpToTime(targetTime);
            return;
        end

        if S.nFrames < 1
            return;
        end
        newFrame = round(src.Value);
        newFrame = max(1, min(S.nFrames, newFrame));
        if newFrame == S.frame
            return;
        end
        S.frame = newFrame;
        goField.Value = S.frame;
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
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
            m = modalities{ii};
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
            goField.Limits = [1 S.nFrames];

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
            goField.Limits = [1 S.nFrames];

            if S.nFrames <= 1
                frameSlider.Limits = [1 2];
            else
                frameSlider.Limits = [1 S.nFrames];
            end
            frameSlider.Value  = 1;
            frameSlider.Enable = 'on';
        end

        goField.Value = 1;
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
        goField.Value = S.frame;
        setSliderFromFrame();
        drawAll();
        updateTimeDisplay();
        syncHsiToTime(timeForFrame(S.frame));
    end

    function syncHsiToTime(tTarget)
        if isempty(S.hsiEvents)
            return;
        end

        if nargin < 1 || isempty(tTarget)
            tTarget = S.hsiEvents(1).time;
        end

        diffs = abs([S.hsiEvents.time] - tTarget);
        [~, idxEvt] = min(diffs);
        evt = S.hsiEvents(idxEvt);

        if strcmp(S.currentHsi.sensor, evt.sensor) && ...
           isdatetime(S.currentHsi.time) && S.currentHsi.time == evt.time
            return;
        end

        switch evt.sensor
            case 'CERB'
                loadCerbFromPath('LWIR', evt.path);
            case 'MX20'
                loadMX20FromHdr(evt.path);
        end

        S.currentHsi = struct('sensor', evt.sensor, 'time', evt.time);
    end

    function setBehavior(val)
        switch val
            case 'Hold last'
                S.behavior = 'hold';
            case 'Show missing'
                S.behavior = 'missing';
            case 'Loop'
                S.behavior = 'loop';
        end
        drawAll();
        updateTimeDisplay();
    end

    function saveMontage()
        % Uses fridge_read_frame + your toUint8/padTo helpers
        names = {'LWIR','MWIR','SWIR'; 'MONO','VIS-COLOR',''};
        tiles = cell(2,3);
        maxH  = 0;
        maxW  = 0;

        for r = 1:2
            for c = 1:3
                name = names{r,c};
                if isempty(name)
                    tiles{r,c} = uint8(255);
                    continue;
                end
                if ~S.exists(name)
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
    end

    %======================== DRAWING / IO (FRIDGE) ========================
    function drawAll()
        tCurrent = timeForFrame(S.frame);
        for i = 1:numel(modalities)
            m  = modalities{i};
            ax = axMap(m);
            cla(ax);

            if ~S.exists(m)
                axis(ax,'off');
                title(ax, sprintf('%s — Missing file', m));
                continue;
            end

            [fEff, status] = effectiveFrame(m, S.frame, tCurrent);
            maxF = S.maxFrames(m);
            if isnan(fEff)
                axis(ax,'off');
                title(ax, sprintf('%s — No frame (max %g)', m, maxF));
                continue;
            end

            img = fridge_read_frame(m, fEff, S.hdrs, S.files);

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
            title(ax, sprintf('%s — Frame %d/%d%s', m, fEff, maxF, note));
        end
    end

    function [fEff, status] = effectiveFrame(modality, fReq, targetTime)
        if nargin < 3
            targetTime = [];
        end

        maxF = S.maxFrames(modality);
        if isnan(maxF) || maxF < 1
            fEff   = NaN;
            status = 'missing';
            return;
        end

        % Time-driven: pick the nearest timestamp for this modality and hold
        % the first/last frame outside its bounds so all panes stay aligned.
        if strcmp(S.sliderMode,'time') && hasFridgeTimes() && ...
                isKey(S.fridgeTimesMap, modality)
            tVec = S.fridgeTimesMap(modality);
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
        S.currentHsi   = struct('sensor','', 'time', NaT);

        lblFile.Text   = 'Filename: -';
        lblFrames.Text = 'Frames: -';
        lblMem.Text    = '';
        lblPixel.Text  = 'Pixel: -';
        lblValue.Text  = 'Value: -';
        lblTime.Text   = 'Time: -';

        for i = 1:numel(modalities)
            ax = axMap(modalities{i});
            cla(ax);
            axis(ax,'off');
            title(ax, modalities{i});
        end

        btnPrev.Enable     = 'off';
        btnNext.Enable     = 'off';
        btnSave.Enable     = 'off';
        goField.Enable     = 'off';
        btnGo.Enable       = 'off';
        goField.Value      = 1;
        goField.Limits     = [1 Inf];

        frameSlider.Enable = 'off';
        frameSlider.Limits = [1 2];
        frameSlider.Value  = 1;

        cla(cerbAxLWIR); title(cerbAxLWIR,'CERB LWIR');
        cla(cerbAxVNIR); title(cerbAxVNIR,'CERB VNIR');
        cla(mxAx);       title(mxAx,'MX20 SW');

        S.cerb = struct('LWIR',[],'VNIR',[]);
        S.mx20 = struct('hdr',[],'ctx',[]);
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
