classdef RMBV_Export
%RMBV_Export Off-screen montage snapshot and video exporter for RMBV.
%
% Usage:
%   RMBV_Export.exportSnapshot(S, opts);
%   RMBV_Export.exportVideo(S, opts);
%   frame = RMBV_Export.renderMontageFrame(S, timeOrPlan);
%
% The exporter tiles normalized FRIDGE and HSI panes into an RGB montage,
% overlays timestamps/labels, and writes frames directly via VideoWriter
% (preferring MPEG-4, falling back to Motion JPEG AVI). Frames are rendered
% off-screen using arrays, not UI captures, so export works even when the
% viewer UI is minimized.
%
% Key options (see individual methods for full defaults):
%   opts.targetSize - [H W] output resolution (default [1080 1920])
%   opts.fps        - Video frame rate (default 15)
%   opts.frameStep  - Export every N FRIDGE frames (default 1)
%   opts.timeStep   - Export every N seconds when timestamps exist
%   opts.outputPath - Destination file for snapshot or video

methods(Static)
    function exportSnapshot(S, opts)
        if nargin < 2
            opts = struct();
        end
        opts = applyDefaults(opts, struct('time', [], 'targetSize', [1080 1920], ...
                                          'outputPath', '', 'includeLabels', true));

        layoutSpec = computeLayoutSpec(S, opts.targetSize);
        plan = buildSinglePlan(S, opts.time);
        cache = struct();
        [frameRGB, ~] = renderMontageFrameInternal(S, plan, layoutSpec, cache, opts.includeLabels);

        if isempty(opts.outputPath)
            return;
        end

        imwrite(frameRGB, opts.outputPath);
    end

    function exportVideo(S, opts)
        if nargin < 2
            opts = struct();
        end
        opts = applyDefaults(opts, struct('targetSize', [1080 1920], ...
                                          'fps', 15, ...
                                          'frameStep', 1, ...
                                          'timeStep', [], ...
                                          'times', [], ...
                                          'timeBase', 'master', ...
                                          'outputPath', '', ...
                                          'includeLabels', true, ...
                                          'parentFigure', []));

        if isempty(opts.outputPath)
            error('RMBV_Export:MissingOutput', 'opts.outputPath is required for exportVideo.');
        end

        % Ensure any lingering dialogs are always cleaned up when we exit.
        dialogGuard = onCleanup(@()forceCloseExportDialogs());

        prepDlg = openProgress(opts.parentFigure, 'Preparing export plan...', 1, true);
        prepDlgClose = prepDlg; % capture current handle for cleanup
        prepCloser = onCleanup(@()closeProgress(prepDlgClose));

        layoutSpec = computeLayoutSpec(S, opts.targetSize);
        plans = buildPlanList(S, opts);
        nFrames = numel(plans);
        if nFrames < 1
            error('RMBV_Export:NoFrames', 'No frames were selected for export.');
        end

        closeProgress(prepDlg); %#ok<NASGU>

        writer = makeVideoWriter(opts.outputPath, opts.fps);
        open(writer);

        dlg = openProgress(opts.parentFigure, 'Exporting montage video...', nFrames, false);
        dlgClose = dlg; % capture handle value before it can change
        dlgCloser = onCleanup(@()closeProgress(dlgClose));
        cache = struct();
        cancelHit = false;
        started = tic;
        try
            for ii = 1:nFrames
                [frameRGB, cache] = renderMontageFrameInternal(S, plans(ii), layoutSpec, cache, opts.includeLabels);
                writeVideo(writer, frameRGB);
                [cancelHit, dlg] = updateProgress(dlg, ii, nFrames, started);
                if cancelHit
                    break;
                end
            end
        catch err
            close(writer);
            closeProgress(dlgClose);
            rethrow(err);
        end

        close(writer);
        closeProgress(dlgClose);
        forceCloseExportDialogs();
        if cancelHit
            warning('RMBV_Export:Canceled', 'Export canceled by user. Video may be incomplete.');
        end
    end

    function [nFrames, times, meta] = estimateFrameCount(S, opts)
        if nargin < 2
            opts = struct();
        end
        opts = applyDefaults(opts, struct('frameStep', 1, 'timeStep', [], 'times', [], ...
                                          'fps', 15, 'timeBase', 'master'));
        [times, meta] = deriveTimes(S, opts);
        nFrames = numel(times);
    end

    function [frameRGB, cacheOut] = renderMontageFrame(S, timeOrPlan, layoutSpec, cacheIn)
        if nargin < 3 || isempty(layoutSpec)
            layoutSpec = computeLayoutSpec(S, [1080 1920]);
        end
        if nargin < 4
            cacheIn = struct();
        end

        if isstruct(timeOrPlan) && isfield(timeOrPlan, 'frameMap')
            plan = timeOrPlan;
        else
            plan = buildSinglePlan(S, timeOrPlan);
        end
        [frameRGB, cacheOut] = renderMontageFrameInternal(S, plan, layoutSpec, cacheIn, true);
    end

    end
end

%--------------------------------------------------------------------------
function idx = selectHsiScanIndex(varargin)
    % selectHsiScanIndex(currentTimeSec, totalDurationSec, nItems)
    % Extra trailing arguments are ignored to tolerate legacy call sites.
    currentTimeSec    = getArg(varargin, 1, 0);
    totalDurationSec  = getArg(varargin, 2, 0);
    nItems            = getArg(varargin, 3, NaN);

    if isempty(nItems) || nItems < 1
        idx = NaN;
        return;
    end
    if nItems == 1
        idx = 1;
        return;
    end

    if isempty(totalDurationSec) || ~isfinite(totalDurationSec)
        totalDurationSec = 0;
    end
    if isempty(currentTimeSec) || ~isfinite(currentTimeSec)
        currentTimeSec = 0;
    end

    if totalDurationSec <= 0
        idx = 1;
        return;
    end

    tinyEps = eps;
    u = min(max(currentTimeSec, 0), max(totalDurationSec - tinyEps, 0));
    delta = totalDurationSec / nItems;
    if delta <= 0
        idx = 1;
        return;
    end
    idx = floor(u / delta) + 1;
    idx = min(max(idx, 1), nItems);
end

%--------------------------------------------------------------------------
function [timelineSeconds, totalDurationSec] = computeExportTimelineSeconds(times, opts, varargin)
    timelineSeconds = [];
    totalDurationSec = 0;
    if isempty(times)
        return;
    end

    if isdatetimeVector(times) && ~any(isnat(times))
        timelineSeconds = seconds(times - times(1));
    elseif isnumeric(times) && isvector(times)
        times = times(:);
        timelineSeconds = times - times(1);
    else
        fps = getfieldOr(opts, 'fps', 1);
        if isempty(fps) || ~isfinite(fps) || fps <= 0
            fps = 1;
        end
        timelineSeconds = (0:(numel(times)-1))' ./ fps;
    end

    if ~isempty(timelineSeconds)
        totalDurationSec = max(timelineSeconds(:));
    end

end

%--------------------------------------------------------------------------
function hsiScanSchedule = buildHsiScanSchedule(S)
    hsiScanSchedule = containers.Map('KeyType','char','ValueType','any');
    if ~isfield(S, 'hsiGroupsMap') || isempty(S.hsiGroupsMap)
        return;
    end

    keys = S.hsiGroupsMap.keys;
    for kk = 1:numel(keys)
        paneKey = keyify(keys{kk});
        data = getOr(S.hsiGroupsMap, paneKey, []);
        if isempty(data) || ~isfield(data, 'groups') || isempty(data.groups)
            continue;
        end

        entries = struct('groupIdx', {}, 'scanIdx', {}, 'item', {}, 'time', {});
        groups = data.groups;
        for gg = 1:numel(groups)
            grp = groups(gg);
            if ~isfield(grp, 'items') || isempty(grp.items)
                continue;
            end
            for jj = 1:numel(grp.items)
                itm = grp.items(jj);
                tVal = getfieldOr(grp, 'time', []);
                entries(end+1) = struct('groupIdx', gg, 'scanIdx', jj, 'item', itm, 'time', tVal); %#ok<AGROW>
            end
        end

        if ~isempty(entries)
            entries = sortHsiEntries(entries);
            hsiScanSchedule(paneKey) = struct('entries', entries);
        end
    end
end

%--------------------------------------------------------------------------
function entriesOut = sortHsiEntries(entriesIn)
    entriesOut = entriesIn;
    if isempty(entriesIn)
        return;
    end

    scanIds = arrayfun(@(e)getfieldOr(e.item, 'scanId', NaN), entriesIn);
    paths = string(arrayfun(@(e)getfieldOr(e.item, 'path', ''), entriesIn, 'UniformOutput', false));
    orig = (1:numel(entriesIn))';
    tbl = table(scanIds(:), paths(:), orig, 'VariableNames', {'scanId','path','orig'});
    tbl = sortrows(tbl, {'scanId','path','orig'});
    ord = tbl.orig;
    entriesOut = entriesIn(ord);
end

%==========================================================================
function plan = buildSinglePlan(S, t)
    plan = struct('time', [], ...
                  'frameMap', containers.Map('KeyType','char','ValueType','double'), ...
                  'holdMap', containers.Map('KeyType','char','ValueType','char'), ...
                  'hsiMap', containers.Map('KeyType','char','ValueType','any'), ...
                  'hsiIndexMap', containers.Map('KeyType','char','ValueType','double'));
    if nargin < 2
        t = [];
    end
    if isempty(t)
        t = currentTime(S);
    end

    plan.time = t;
    activePanels = getActivePanels(S);
    for ii = 1:numel(activePanels)
        pane = activePanels(ii);
        if pane.isHsi
            [evt, idx] = pickHsiEvent(S, t, pane.sensor, pane.modality);
            plan.hsiMap(pane.key) = evt;
            plan.hsiIndexMap(pane.key) = idx;
            continue;
        end
        [idx, holdState] = effectiveFrameForExport(S, pane.key, t);
        plan.frameMap(pane.key) = idx;
        plan.holdMap(pane.key)  = holdState;
    end
end

%--------------------------------------------------------------------------
function plans = buildPlanList(S, opts)
    times = opts.times;
    if isempty(times)
        [times, ~] = deriveTimes(S, opts);
    end
    if isempty(times)
        plans = struct([]);
        return;
    end
    plans(numel(times)) = struct('time', [], 'frameMap', [], 'holdMap', [], 'hsiMap', [], 'hsiIndexMap', []);
    lastHsiEvt = containers.Map('KeyType','char','ValueType','any');
    lastHsiIdx = containers.Map('KeyType','char','ValueType','double');
    [timelineSeconds, totalDurationSec] = computeExportTimelineSeconds(times, opts);
    hsiScanSchedule = buildHsiScanSchedule(S);
    for ii = 1:numel(times)
        plan = buildSinglePlan(S, times(ii));
        if isa(plan.hsiMap, 'containers.Map')
            keys = plan.hsiMap.keys;
            for kCell = keys
                k = kCell{1};
                evtVal = plan.hsiMap(k);
                if isempty(evtVal) && isKey(lastHsiEvt, k)
                    plan.hsiMap(k) = lastHsiEvt(k);
                    if isa(plan.hsiIndexMap, 'containers.Map') && isKey(lastHsiIdx, k)
                        plan.hsiIndexMap(k) = lastHsiIdx(k);
                    end
                elseif ~isempty(evtVal)
                    lastHsiEvt(k) = evtVal;
                    if isa(plan.hsiIndexMap, 'containers.Map') && isKey(plan.hsiIndexMap, k)
                        lastHsiIdx(k) = plan.hsiIndexMap(k);
                    end
                end
            end
        end
        plan = applyHsiExportCycling(S, plan, hsiScanSchedule, timelineSeconds, totalDurationSec, ii);
        plans(ii) = plan;
    end
end

%--------------------------------------------------------------------------
function planOut = applyHsiExportCycling(varargin)
    % applyHsiExportCycling(S, planIn, hsiScanSchedule, timelineSeconds, totalDurationSec, frameIdx, ...)
    % Accepts optional trailing arguments to tolerate legacy callers.
    S                 = getArg(varargin, 1, struct());
    planIn            = getArg(varargin, 2, struct());
    hsiScanSchedule   = getArg(varargin, 3, containers.Map('KeyType','char','ValueType','any'));
    timelineSeconds   = getArg(varargin, 4, []);
    totalDurationSec  = getArg(varargin, 5, []);
    frameIdx          = getArg(varargin, 6, NaN);

    planOut = planIn;
    if ~isstruct(planIn) || ~isfield(planIn, 'hsiMap') || ...
            ~isa(planIn.hsiMap, 'containers.Map') || planIn.hsiMap.Count < 1 || ...
            ~isstruct(S) || ~isfield(S, 'hsiGroupsMap') || isempty(S.hsiGroupsMap)
        return;
    end

    if nargin < 3 || isempty(hsiScanSchedule)
        hsiScanSchedule = containers.Map('KeyType','char','ValueType','any');
    end

    paneKeys = planIn.hsiMap.keys;
    for kCell = paneKeys
        paneKey = keyify(kCell{1});
        schedule = getOr(hsiScanSchedule, paneKey, struct('entries', struct([])));
        entries = getfieldOr(schedule, 'entries', struct([]));
        if isempty(entries)
            planOut.hsiMap(paneKey) = [];
            continue;
        end

        nItems = numel(entries);
        if nItems < 1
            planOut.hsiMap(paneKey) = [];
            continue;
        end

        currentTimeSec = 0;
        if nargin >= 4 && ~isempty(timelineSeconds) && frameIdx >= 1 && frameIdx <= numel(timelineSeconds)
            currentTimeSec = timelineSeconds(frameIdx);
        end
        if nargin < 5 || isempty(totalDurationSec) || ~isfinite(totalDurationSec)
            totalDurationSec = 0;
        end

        scanIdx = selectHsiScanIndex(currentTimeSec, totalDurationSec, nItems);
        entry = entries(scanIdx);
        item = entry.item;
        evtTime = getfieldOr(entry, 'time', []);
        evt = struct('sensor', item.sensor, 'modality', item.modality, 'path', item.path, ...
            'time', evtTime, 'groupIdx', entry.groupIdx, 'scanIdx', entry.scanIdx, ...
            'scanCount', nItems);
        if isfield(item, 'label')
            evt.scanLabel = item.label;
        end
        planOut.hsiMap(paneKey) = evt;
        if isa(planOut.hsiIndexMap, 'containers.Map')
            planOut.hsiIndexMap(paneKey) = entry.scanIdx;
        end
    end
end

%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
function [times, meta] = deriveTimes(S, opts)
    meta = struct('timeBase', getfieldOr(opts, 'timeBase'), ...
                  'masterModality', '', 'startTime', NaT, 'endTime', NaT, ...
                  'cadenceSeconds', NaN);
    times = opts.times;
    if ~isempty(times)
        if isdatetimeVector(times)
            meta.startTime = times(1);
            meta.endTime   = times(end);
        end
        return;
    end

    if ~isfield(opts, 'timeBase') || isempty(opts.timeBase)
        opts.timeBase = 'master';
    end

    switch lower(string(opts.timeBase))
        case "fixed"
            [times, meta] = deriveFixedTimes(S, opts, meta);
        otherwise
            [times, meta] = deriveMasterTimes(S, opts, meta);
    end
end

%--------------------------------------------------------------------------
function [times, meta] = deriveMasterTimes(S, opts, meta)
    if nargin < 3 || isempty(meta)
        meta = struct('timeBase','master','masterModality','', ...
                     'startTime', NaT, 'endTime', NaT, 'cadenceSeconds', NaN);
    end

    stats = collectTimelineStats(S);
    meta.masterModality = stats.densestModality;
    meta.startTime = stats.tMin;
    meta.endTime   = stats.tMax;

    dtExport = stats.dt;
    dtMinClamp = seconds(1/60); % avoid runaway >60 fps exports
    if isnat(meta.startTime) || isnat(meta.endTime) || isnan(seconds(dtExport))
        times = fallbackTimeline(S, opts);
    else
        spanSeconds = seconds(meta.endTime - meta.startTime);
        dtSeconds = seconds(dtExport);
        if ~isfinite(spanSeconds) || spanSeconds <= 0 || ~isfinite(dtSeconds)
            times = fallbackTimeline(S, opts);
        else
            if dtSeconds < seconds(dtMinClamp)
                dtSeconds = seconds(dtMinClamp);
            end

            maxFramesSafe = 1e6; % guard against enormous allocations
            nEstimate = floor(spanSeconds / dtSeconds) + 1;
            if nEstimate > maxFramesSafe
                dtSeconds = spanSeconds / maxFramesSafe;
                dtSeconds = max(dtSeconds, seconds(dtMinClamp));
                nEstimate = floor(spanSeconds / dtSeconds) + 1;
            end

            times = meta.startTime + seconds((0:(nEstimate-1))' .* dtSeconds);
            if ~isempty(times) && times(end) < meta.endTime
                times(end+1,1) = meta.endTime; %#ok<AGROW>
            end
            dtExport = seconds(dtSeconds);
        end
    end

    if isempty(times)
        return;
    end
    if isempty(opts.frameStep) || opts.frameStep < 1
        opts.frameStep = 1;
    end
    idx = 1:opts.frameStep:numel(times);
    times = times(idx);
    meta.startTime = times(1);
    meta.endTime   = times(end);
    meta.cadenceSeconds = seconds(dtExport);
end

%--------------------------------------------------------------------------
function [times, meta] = deriveFixedTimes(S, opts, meta)
    if nargin < 3 || isempty(meta)
        meta = struct('timeBase','fixed','masterModality','', ...
                      'startTime', NaT, 'endTime', NaT, 'cadenceSeconds', NaN);
    end
    times = [];
    dtSeconds = [];
    if ~isempty(opts.timeStep)
        dtSeconds = opts.timeStep;
    elseif isfield(opts, 'fps') && ~isempty(opts.fps) && opts.fps > 0
        dtSeconds = 1 / opts.fps;
    end

    [tStart, tEnd] = fridgeTimeBounds(S);
    if isnat(tStart) || isnat(tEnd) || isempty(dtSeconds) || dtSeconds <= 0
        times = fallbackTimeline(S, opts);
    else
        times = tStart:seconds(dtSeconds):tEnd;
        if ~isempty(times) && times(end) < tEnd
            times(end+1) = tEnd; %#ok<AGROW>
        end
    end

    if isempty(times)
        return;
    end
    if isempty(opts.frameStep) || opts.frameStep < 1
        opts.frameStep = 1;
    end
    idx = 1:opts.frameStep:numel(times);
    times = times(idx);
    meta.startTime = times(1);
    meta.endTime   = times(end);
    meta.cadenceSeconds = dtSeconds;
end

%--------------------------------------------------------------------------
function [masterMod, masterTimes] = selectMasterTimebase(S)
    stats = collectTimelineStats(S);
    masterMod = stats.densestModality;
    masterTimes = stats.densestTimes;
    if ~isempty(masterTimes)
        masterTimes = unique(masterTimes);
        masterTimes = sort(masterTimes);
    end
end

%--------------------------------------------------------------------------
function stats = collectTimelineStats(S)
    stats = struct('tMin', NaT, 'tMax', NaT, 'dt', seconds(NaN), ...
                   'densestModality','', 'densestTimes', datetime.empty(0,1));
    if ~isfield(S, 'fridgeTimesMap') || ~isa(S.fridgeTimesMap, 'containers.Map')
        return;
    end

    keys = S.fridgeTimesMap.keys;
    bestDt = Inf;
    for ii = 1:numel(keys)
        k = keyify(keys{ii});
        if isfield(S, 'exists') && isa(S.exists, 'containers.Map') && ...
                isKey(S.exists, k) && ~S.exists(k)
            continue;
        end

        rawTimes = S.fridgeTimesMap(k);
        maxF = getOr(S.maxFrames, k, numel(rawTimes));
        tVec = sanitizeTimes(rawTimes, maxF);
        if isempty(tVec)
            continue;
        end

        if isnat(stats.tMin) || tVec(1) < stats.tMin
            stats.tMin = tVec(1);
        end
        if isnat(stats.tMax) || tVec(end) > stats.tMax
            stats.tMax = tVec(end);
        end

        dts = seconds(diff(tVec));
        dts = dts(~isnan(dts) & dts > 0);
        if isempty(dts)
            continue;
        end
        medDt = median(dts);
        if medDt < bestDt
            bestDt = medDt;
            stats.densestModality = k;
            stats.densestTimes = tVec;
        end
    end

    hsiTimes = collectHsiTimes(S);
    if ~isempty(hsiTimes)
        if isnat(stats.tMin) || hsiTimes(1) < stats.tMin
            stats.tMin = hsiTimes(1);
        end
        if isnat(stats.tMax) || hsiTimes(end) > stats.tMax
            stats.tMax = hsiTimes(end);
        end
    end

    if isfinite(bestDt)
        stats.dt = seconds(bestDt);
    end
end

%--------------------------------------------------------------------------
function times = fallbackTimeline(S, opts)
    times = datetime.empty(0,1);
    if isdatetimeVector(S.timelineTimes)
        if isempty(opts.frameStep) || opts.frameStep < 1
            opts.frameStep = 1;
        end
        idx = 1:opts.frameStep:numel(S.timelineTimes);
        times = S.timelineTimes(idx);
    elseif isfield(S, 'nFrames') && S.nFrames > 0
        if isempty(opts.frameStep) || opts.frameStep < 1
            opts.frameStep = 1;
        end
        n = numel(1:opts.frameStep:S.nFrames);
        times = NaT(n,1);
    end
end

%--------------------------------------------------------------------------
function [tStart, tEnd] = fridgeTimeBounds(S)
    tStart = NaT;
    tEnd   = NaT;
    if isfield(S, 'fridgeTimesMap') && isa(S.fridgeTimesMap, 'containers.Map')
        keys = S.fridgeTimesMap.keys;
        for ii = 1:numel(keys)
            tVec = sanitizeTimes(S.fridgeTimesMap(keys{ii}));
            if isempty(tVec)
                continue;
            end
            if isnat(tStart) || tVec(1) < tStart
                tStart = tVec(1);
            end
            if isnat(tEnd) || tVec(end) > tEnd
                tEnd = tVec(end);
            end
        end
    end

    hsiTimes = collectHsiTimes(S);
    if ~isempty(hsiTimes)
        if isnat(tStart) || hsiTimes(1) < tStart
            tStart = hsiTimes(1);
        end
        if isnat(tEnd) || hsiTimes(end) > tEnd
            tEnd = hsiTimes(end);
        end
    end
end

%--------------------------------------------------------------------------
function tVec = collectHsiTimes(S)
    tVec = datetime.empty(0,1);
    if ~isfield(S, 'hsiEvents') || isempty(S.hsiEvents)
        return;
    end
    times = arrayfun(@(e) effectiveHsiTime(S, e), S.hsiEvents);
    times = times(~isnat(times));
    if isempty(times)
        return;
    end
    tVec = sanitizeTimes(sort(times(:)));
end

%--------------------------------------------------------------------------
function tVec = sanitizeTimes(tVec, maxCount)
    if nargin < 2
        maxCount = [];
    end
    if ~isdatetimeVector(tVec)
        tVec = datetime.empty(0,1);
        return;
    end
    tVec = tVec(:);
    tVec = tVec(~isnat(tVec));
    if isempty(tVec)
        return;
    end

    y = year(tVec);
    validMask = isfinite(y) & y >= 0 & y <= 9999;
    tVec = tVec(validMask);
    if isempty(tVec)
        return;
    end

    tVec = sort(tVec);
    if ~isempty(maxCount)
        tVec = tVec(1:min(numel(tVec), maxCount));
    end
end

%--------------------------------------------------------------------------
function layoutSpec = computeLayoutSpec(S, targetSize)
    if nargin < 2 || isempty(targetSize)
        targetSize = [1080 1920];
    end
    activePanels = getActivePanels(S);
    n = numel(activePanels);
    maxCols = 3;
    cols = min(maxCols, max(1, n));
    rows = ceil(max(1, n) / cols);

    layoutSpec = struct();
    layoutSpec.targetSize = targetSize;
    layoutSpec.panels     = activePanels;
    layoutSpec.rows       = rows;
    layoutSpec.cols       = cols;
end

%--------------------------------------------------------------------------
function [frameRGB, cacheOut] = renderMontageFrameInternal(S, plan, layoutSpec, cacheIn, includeLabels)
    if nargin < 5
        includeLabels = true;
    end
    cacheOut = cacheIn;

    H = layoutSpec.targetSize(1);
    W = layoutSpec.targetSize(2);
    tileH = floor(H / layoutSpec.rows);
    tileW = floor(W / layoutSpec.cols);

    frameRGB = uint8(255 * ones(H, W, 3));

    for idx = 1:numel(layoutSpec.panels)
        pane = layoutSpec.panels(idx);
        key = keyify(pane.key);
        r = ceil(idx / layoutSpec.cols);
        c = mod(idx-1, layoutSpec.cols) + 1;
        yStart = (r-1)*tileH + 1;
        yEnd   = min(r*tileH, H);
        xStart = (c-1)*tileW + 1;
        xEnd   = min(c*tileW, W);

        [tile, cacheOut] = renderTile(S, plan, pane, cacheOut);
        tile = letterboxToSize(tile, yEnd - yStart + 1, xEnd - xStart + 1, 0);
        frameRGB(yStart:yEnd, xStart:xEnd, :) = tile;

        if includeLabels
            lbl = buildTileLabel(plan, pane);
            frameRGB = overlayText(frameRGB, [xStart+5, yStart+5], lbl);
            holdLbl = holdLabel(plan, key);
            if ~isempty(holdLbl)
                frameRGB = overlayText(frameRGB, [xStart+5, yEnd-25], holdLbl);
            end
        end
    end

    if includeLabels
        tsLbl = formatTimestamp(plan.time);
        frameRGB = overlayText(frameRGB, [10, H-30], tsLbl);
    end
end

%--------------------------------------------------------------------------
function [tileRGB, cacheOut] = renderTile(S, plan, pane, cacheIn)
    cacheOut = cacheIn;
    cacheKey = makeCacheFieldName(['pane_' pane.key]);
    cached = struct('idx', NaN, 'img', []);
    if isfield(cacheOut, cacheKey)
        cached = cacheOut.(cacheKey);
    end

    if pane.isHsi
        evt = getOr(plan.hsiMap, pane.key, []);
        evtKey = hsiCacheKey(evt);
        cachedEvtKey = '';
        if isfield(cached, 'evtKey')
            cachedEvtKey = cached.evtKey;
        end
        if isempty(evt) && ~isempty(cached) && isfield(cached, 'img') && ~isempty(cached.img)
            tileRGB = cached.img;
            return;
        end
        if ~isempty(cached) && strcmp(cachedEvtKey, evtKey)
            tileRGB = cached.img;
            return;
        end
        tileRGB = renderHsiTile(evt);
        cacheOut.(cacheKey) = struct('evtKey', evtKey, 'img', tileRGB);
        return;
    end

    key = pane.key;
    fIdx = NaN;
    if isfield(plan, 'frameMap') && isa(plan.frameMap, 'containers.Map') && isKey(plan.frameMap, key)
        fIdx = plan.frameMap(key);
    end
    if ~isnan(cached.idx) && cached.idx == fIdx && ~isempty(cached.img)
        tileRGB = cached.img;
        return;
    end

    tileRGB = renderFridgeTile(S, key, fIdx);
    cacheOut.(cacheKey) = struct('idx', fIdx, 'img', tileRGB);
end

%--------------------------------------------------------------------------
function tileRGB = renderFridgeTile(S, modality, fIdx)
    if isnan(fIdx) || fIdx < 1
        tileRGB = uint8(255 * ones(100, 160, 3));
        return;
    end

    try
        img = fridge_read_frame(modality, fIdx, S.hdrs, S.files);
    catch
        tileRGB = uint8(255 * ones(100, 160, 3));
        return;
    end

    tileRGB = normalizePane(img, modality);
end

%--------------------------------------------------------------------------
function tileRGB = renderHsiTile(evt)
    if isempty(evt) || ~isstruct(evt) || ~isfield(evt, 'sensor') || isempty(evt.sensor)
        tileRGB = uint8(255 * ones(100, 160, 3));
        return;
    end

    try
        switch upper(evt.sensor)
            case 'CERB'
                img = loadCerberusContext(evt.path);
                img = rot90(img, -1);
            case 'MX20'
                img = loadCerberusContext(evt.path);
                img = rot90(img, -1);
            case 'FAST'
                img = loadCerberusContext(evt.path);
                img = rot90(img, -1);
            otherwise
                img = [];
        end
    catch
        img = [];
    end

    tileRGB = normalizePane(img, evtModality(evt));
end

%--------------------------------------------------------------------------
function lbl = paneTitleFromEvent(evt, pane)
    lbl = 'HSI — none';
    if nargin < 2
        pane = struct('sensor','', 'modality','');
    end
    if isempty(evt)
        if isfield(pane,'sensor') && ~isempty(pane.sensor)
            lbl = pane.sensor;
            if isfield(pane,'modality') && ~isempty(pane.modality)
                lbl = sprintf('%s %s', lbl, pane.modality);
            end
        end
        return;
    end

    lbl = evt.sensor;
    if isfield(evt,'modality') && ~isempty(evt.modality)
        lbl = sprintf('%s %s', lbl, evt.modality);
    elseif isfield(pane,'modality') && ~isempty(pane.modality)
        lbl = sprintf('%s %s', lbl, pane.modality);
    end
    if isfield(evt,'scanIdx') && ~isempty(evt.scanIdx) && isfinite(evt.scanIdx)
        if isfield(evt,'scanCount') && ~isempty(evt.scanCount) && isfinite(evt.scanCount)
            lbl = sprintf('%s — Scan %d/%d', lbl, evt.scanIdx, evt.scanCount);
        else
            lbl = sprintf('%s — Scan %d', lbl, evt.scanIdx);
        end
    end
    if isfield(evt,'scanLabel') && ~isempty(evt.scanLabel)
        lbl = sprintf('%s — %s', lbl, evt.scanLabel);
    end
    if isfield(evt,'time') && isdatetime(evt.time) && ~isnat(evt.time)
        lbl = sprintf('%s — %s', lbl, datestr(evt.time,'yyyy-mm-dd HH:MM:SS.FFF'));
    end
end

%--------------------------------------------------------------------------
function modOut = evtModality(evt)
    modOut = 'HSI';
    if isstruct(evt) && isfield(evt,'modality') && ~isempty(evt.modality)
        modOut = keyify(evt.modality);
    end
end

%--------------------------------------------------------------------------
function lbl = buildTileLabel(plan, pane)
    key = pane.key;
    if pane.isHsi
        evt = getOr(plan.hsiMap, key, []);
        lbl = paneTitleFromEvent(evt, pane);
        return;
    end

    lbl = key;
    if isfield(plan, 'frameMap') && isa(plan.frameMap, 'containers.Map') && isKey(plan.frameMap, key)
        fIdx = plan.frameMap(key);
        if ~isnan(fIdx)
            lbl = sprintf('%s — Frame %d', key, fIdx);
        end
    end
end

%--------------------------------------------------------------------------
function lbl = holdLabel(plan, key)
    lbl = '';
    if ~isfield(plan, 'holdMap') || ~isa(plan.holdMap, 'containers.Map') || ~isKey(plan.holdMap, key)
        return;
    end
    state = plan.holdMap(key);
    switch string(state)
        case "pre"
            lbl = 'HOLD first frame';
        case "post"
            lbl = 'HOLD last frame';
        otherwise
            lbl = '';
    end
end

%--------------------------------------------------------------------------
function rgb = normalizePane(img, modality)
    if nargin < 2
        modality = '';
    end
    if isempty(img)
        rgb = uint8(255 * ones(100, 160, 3));
        return;
    end

    if strcmpi(modality, 'VIS-COLOR') && ndims(img)==3 && size(img,3)==3
        rgb = normalizeColor(img);
        return;
    end

    img8 = toUint8(img);
    rgb  = repmat(img8, [1 1 3]);
end

%--------------------------------------------------------------------------
function rgb = normalizeColor(img)
    if isa(img, 'uint8')
        rgb = img;
        return;
    end
    img = double(img);
    mx = max(img(:));
    if mx > 0
        img = img ./ mx;
    end
    img = min(max(img, 0), 1);
    rgb = uint8(round(img * 255));
end

%--------------------------------------------------------------------------
function [idx, holdState] = effectiveFrameForExport(S, modality, targetTime)
    modality = keyify(modality);
    holdState = 'none';
    maxF = getOr(S.maxFrames, modality, NaN);
    if isnan(maxF) || maxF < 1
        idx = NaN;
        return;
    end

    if isdatetimeVector(getOr(S.fridgeTimesMap, modality, datetime.empty(0,1)))
        tVec = S.fridgeTimesMap(modality);
        tVec = tVec(~isnat(tVec));
        if isempty(targetTime) || ~isdatetime(targetTime)
            targetTime = currentTime(S);
        end
        if isdatetime(targetTime)
            if targetTime <= tVec(1)
                idx = 1;
                holdState = 'pre';
            elseif targetTime >= tVec(end)
                idx = numel(tVec);
                holdState = 'post';
            else
                [~, idxSel] = min(abs(tVec - targetTime));
                idx = idxSel;
            end
            idx = min(max(1, idx), maxF);
            return;
        end
    end

    if isfield(S, 'behavior') && strcmp(S.behavior, 'loop')
        idx = mod(max(0, round(currentFrame(S))-1), maxF) + 1;
    else
        idx = min(max(1, round(currentFrame(S))), maxF);
    end
end

%--------------------------------------------------------------------------
function [evt, idx] = pickHsiEvent(S, targetTime, sensor, modality)
    evt = [];
    idx = NaN;
    if ~isfield(S, 'hsiEvents') || isempty(S.hsiEvents)
        return;
    end

    mask = true(numel(S.hsiEvents), 1);
    if nargin >= 3 && ~isempty(sensor)
        maskSensor = reshape(strcmpi({S.hsiEvents.sensor}, sensor), [], 1);
        mask = mask & maskSensor;
    end
    if nargin >= 4 && ~isempty(modality)
        maskMod = reshape(strcmpi({S.hsiEvents.modality}, modality), [], 1);
        mask = mask & maskMod;
    end
    evtList = S.hsiEvents(mask);
    if isempty(evtList)
        return;
    end

    effTimesAll = arrayfun(@(e) effectiveHsiTime(S, e), evtList);
    validMask = ~isnat(effTimesAll);
    if ~any(validMask)
        return;
    end
    effTimes = effTimesAll(validMask);
    evtValid = evtList(validMask);
    if nargin < 2 || isempty(targetTime) || ~isdatetime(targetTime)
        targetTime = currentTime(S);
    end
    if isempty(targetTime) || ~isdatetime(targetTime) || (isdatetime(targetTime) && any(isnat(targetTime)))
        targetTime = effTimes(1);
    end
    if isempty(targetTime) || ~isdatetime(targetTime)
        [~, idxLocal] = min(effTimes - min(effTimes));
    else
        diffs = abs(effTimes - targetTime);
        if all(isnan(diffs))
            [~, idxLocal] = min(effTimes - min(effTimes));
        else
            [~, idxLocal] = min(diffs);
        end
    end
    idxLocal = min(max(1, idxLocal), numel(evtValid));
    evt = evtValid(idxLocal);

    idxCandidates = find(mask);
    idxValid = idxCandidates(validMask);
    if numel(idxValid) >= idxLocal
        idx = idxValid(idxLocal);
    end
end

%--------------------------------------------------------------------------
function tEff = effectiveHsiTime(S, evt)
    tEff = evt.time;
    if ~isfield(evt,'path') || isempty(evt.path)
        return;
    end
    key = evt.path;
    if isfield(S, 'hsiPreciseCache') && isa(S.hsiPreciseCache, 'containers.Map') && isKey(S.hsiPreciseCache, key)
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
end

%--------------------------------------------------------------------------
function img = safeResize(img, targetHW)
    if isempty(img)
        img = uint8(255 * ones(targetHW(1), targetHW(2), 3));
        return;
    end
    if exist('imresize','file') == 2
        img = imresize(img, targetHW);
    else
        img = basicResize(img, targetHW);
    end
end

%--------------------------------------------------------------------------
function imgOut = basicResize(imgIn, targetHW)
    imgOut = uint8(255 * ones(targetHW(1), targetHW(2), size(imgIn,3)));
    sz = size(imgIn);
    scaleY = targetHW(1) / sz(1);
    scaleX = targetHW(2) / sz(2);
    for y = 1:targetHW(1)
        srcY = max(1, min(sz(1), round(y / scaleY)));
        for x = 1:targetHW(2)
            srcX = max(1, min(sz(2), round(x / scaleX)));
            imgOut(y,x,:) = imgIn(srcY, srcX, :);
        end
    end
end

%--------------------------------------------------------------------------
function txt = formatTimestamp(t)
    if isdatetime(t) && ~isnat(t)
        txt = datestr(t, 'yyyy-mm-dd HH:MM:SS.FFF');
    else
        txt = 'Time: (unknown)';
    end
end

%--------------------------------------------------------------------------
function img = overlayText(img, position, txt)
    if ~(ischar(txt) || (isstring(txt) && isscalar(txt)))
        return;
    end
    if exist('insertText','file') == 2
        img = insertText(img, position, char(txt), 'FontSize', 14, ...
            'TextColor', 'white', 'BoxColor', 'black', 'BoxOpacity', 0.6);
        return;
    end

    hFig = figure('Visible','off','Position',[0 0 size(img,2) size(img,1)]);
    ax = axes('Parent', hFig);
    imshow(img, 'Parent', ax);
    text(ax, position(1), position(2), char(txt), 'Color','w', 'FontSize', 12, ...
        'BackgroundColor','k','Margin',1);
    frame = getframe(ax);
    img = frame.cdata;
    close(hFig);
end

%--------------------------------------------------------------------------
function imgOut = letterboxToSize(imgIn, targetH, targetW, padValue)
    if nargin < 4 || isempty(padValue)
        padValue = 0;
    end
    if isempty(imgIn)
        imgIn = uint8(255 * ones(1,1,3));
    end
    sz = size(imgIn);
    if numel(sz) < 3
        sz(3) = 1;
    end
    scale = min(targetH / sz(1), targetW / sz(2));
    scale = max(scale, eps);
    newH = max(1, floor(sz(1) * scale));
    newW = max(1, floor(sz(2) * scale));
    resized = safeResize(imgIn, [newH newW]);

    imgOut = repmat(uint8(padValue), targetH, targetW, 3);
    yOffset = floor((targetH - newH) / 2);
    xOffset = floor((targetW - newW) / 2);
    yRange = (1:newH) + yOffset;
    xRange = (1:newW) + xOffset;
    imgOut(yRange, xRange, :) = resized;
end

%--------------------------------------------------------------------------
function writer = makeVideoWriter(pathOut, fps)
    [~,~,ext] = fileparts(pathOut);
    try
        if strcmpi(ext, '.mp4')
            writer = VideoWriter(pathOut, 'MPEG-4');
        else
            writer = VideoWriter(pathOut);
        end
    catch
        writer = VideoWriter(pathOut, 'Motion JPEG AVI');
    end
    writer.FrameRate = fps;
end

%--------------------------------------------------------------------------
function dlg = openProgress(parentFig, msg, ~, indeterminate)
    if nargin < 4
        indeterminate = false;
    end
    dlg = [];
    if exist('uiprogressdlg','file') == 2
        try
            dlgHandle = uiprogressdlg(parentFig, 'Title','Export', 'Message', msg, ...
                'Cancelable', true, 'Value', 0, 'Indeterminate', logical(indeterminate));
            try
                dlgHandle.Tag = 'RMBV_EXPORT_PROGRESS';
                if isprop(dlgHandle, 'CloseRequestFcn')
                    dlgHandle.CloseRequestFcn = @(src,evt)delete(src);
                end
            catch
            end
            dlg = struct('type','uiprogress', 'h', dlgHandle);
            dlg.h.UserData = struct('Canceled', false);
            progressRegistry('add', dlg);
            return;
        catch
            dlg = [];
        end
    end

    % Fallback to waitbar so users still see progress
    try
        if indeterminate
            h = waitbar(0, msg, 'Name','Export', 'CreateCancelBtn', ...
                'setappdata(gcbf,''Canceling'',true)');
        else
            h = waitbar(0, msg, 'Name','Export', 'CreateCancelBtn', ...
                'setappdata(gcbf,''Canceling'',true)');
        end
        setappdata(h, 'Canceling', false);
        try
            set(h, 'Tag', 'RMBV_EXPORT_PROGRESS', 'CloseRequestFcn', @(src,evt)delete(src));
        catch
        end
        dlg = struct('type','waitbar', 'h', h);
        progressRegistry('add', dlg);
    catch
        dlg = [];
    end
end

%--------------------------------------------------------------------------
function [cancelHit, dlg] = updateProgress(dlg, ii, nFrames, startedTic)
    cancelHit = false;
    if isempty(dlg)
        return;
    end

    val = ii / nFrames;
    etaMsg = '';
    if nargin >= 4 && ~isempty(startedTic)
        elapsed = toc(startedTic);
        rate = ii / max(elapsed, eps);
        etaSeconds = (nFrames - ii) / max(rate, eps);
        etaRounded = max(0, ceil(etaSeconds));
        etaHours = floor(etaRounded / 3600);
        etaMinutes = floor(mod(etaRounded, 3600) / 60);
        etaSecondsRem = mod(etaRounded, 60);
        etaMsg = sprintf(' (ETA %02d:%02d:%02d)', etaHours, etaMinutes, etaSecondsRem);
    end
    msg = sprintf('Writing frame %d of %d%s', ii, nFrames, etaMsg);

    switch dlg.type
        case 'uiprogress'
            try
                dlg.h.Value = val;
                if isprop(dlg.h, 'Message')
                    dlg.h.Message = msg;
                end
                pause(0.001);
                if isfield(dlg.h.UserData, 'Canceled') && dlg.h.UserData.Canceled
                    cancelHit = true;
                end
                if isprop(dlg.h, 'CancelRequested') && dlg.h.CancelRequested
                    cancelHit = true;
                end
            catch
                dlg = [];
            end
        case 'waitbar'
            if ~ishandle(dlg.h)
                dlg = [];
                return;
            end
            try
                waitbar(val, dlg.h, msg);
                cancelHit = getappdata(dlg.h, 'Canceling');
                drawnow limitrate;
            catch
                dlg = [];
            end
    end
end

%--------------------------------------------------------------------------
function closeProgress(dlg)
    if isempty(dlg)
        return;
    end
    reallyCloseProgressHandle(dlg);
    progressRegistry('remove', dlg);
    progressRegistry('closeall');
end

%--------------------------------------------------------------------------
function reallyCloseProgressHandle(dlg)
    if isempty(dlg)
        return;
    end
    try
        switch dlg.type
            case 'uiprogress'
                try
                    close(dlg.h);
                catch
                    delete(dlg.h);
                end
            case 'waitbar'
                if ishandle(dlg.h)
                    try
                        close(dlg.h);
                    catch
                        delete(dlg.h);
                    end
                end
        end
    catch
    end
end

%--------------------------------------------------------------------------
function t = currentTime(S)
    t = [];
    if isfield(S, 'timelineTimes') && isdatetimeVector(S.timelineTimes) && ...
            isfield(S, 'frame') && S.frame >= 1 && S.frame <= numel(S.timelineTimes)
        t = S.timelineTimes(S.frame);
    end
end

%--------------------------------------------------------------------------
function f = currentFrame(S)
    f = 1;
    if isfield(S, 'frame') && ~isempty(S.frame)
        f = S.frame;
    end
end

%--------------------------------------------------------------------------
function tf = isdatetimeVector(v)
    tf = isdatetime(v) && ~isempty(v);
end

%--------------------------------------------------------------------------
function panels = getActivePanels(S)
    panels = struct('key', {}, 'isHsi', {}, 'sensor', {}, 'modality', {});
    if isfield(S, 'exists') && ~isempty(S.exists)
        keys = S.exists.keys;
        for ii = 1:numel(keys)
            k = keyify(keys{ii});
            if getOr(S.exists, k, false)
                panels(end+1) = struct('key', k, 'isHsi', false, 'sensor','', 'modality',''); %#ok<AGROW>
            end
        end
    end

    hsiPanels = getActiveHsiPanels(S);
    if ~isempty(hsiPanels)
        panels = [panels hsiPanels]; %#ok<AGROW>
    end
end

%--------------------------------------------------------------------------
function hsiPanels = getActiveHsiPanels(S)
    hsiPanels = struct('key', {}, 'isHsi', {}, 'sensor', {}, 'modality', {});
    if isfield(S, 'enableHSI') && ~S.enableHSI
        return;
    end

    if isfield(S, 'cerb')
        mods = {'LWIR','VNIR'};
        for ii = 1:numel(mods)
            mod = mods{ii};
            if ~isempty(getfieldOr(S.cerb, mod))
                hsiPanels(end+1) = struct('key', sprintf('HSI-CERB-%s', mod), ...
                    'isHsi', true, 'sensor', 'CERB', 'modality', mod); %#ok<AGROW>
            end
        end
    end

    if isfield(S, 'mx20') && isfield(S.mx20, 'ctx') && ~isempty(S.mx20.ctx)
        hsiPanels(end+1) = struct('key','HSI-MX20-SWIR', 'isHsi', true, ...
            'sensor','MX20', 'modality','SWIR'); %#ok<AGROW>
    end

    if isfield(S, 'fast') && ~isempty(S.fast)
        mods = fieldnames(S.fast);
        for ii = 1:numel(mods)
            mod = keyify(mods{ii});
            if ~isempty(S.fast.(mods{ii}))
                hsiPanels(end+1) = struct('key', sprintf('HSI-FAST-%s', mod), ...
                    'isHsi', true, 'sensor','FAST', 'modality', mod); %#ok<AGROW>
            end
        end
    end

    if isempty(hsiPanels) && isfield(S, 'hsiEvents') && ~isempty(S.hsiEvents)
        sensors = unique({S.hsiEvents.sensor});
        for ii = 1:numel(sensors)
            sensor = sensors{ii};
            modSet = unique({S.hsiEvents(strcmp({S.hsiEvents.sensor}, sensor)).modality});
            if isempty(modSet)
                modSet = {''};
            end
            for jj = 1:numel(modSet)
                mod = modSet{jj};
                hsiPanels(end+1) = struct('key', hsiKey(sensor, mod), 'isHsi', true, ...
                    'sensor', sensor, 'modality', mod); %#ok<AGROW>
            end
        end
    end
end

%--------------------------------------------------------------------------
function val = getfieldOr(s, name, defaultVal)
    if nargin < 3
        defaultVal = [];
    end
    if isstruct(s) && isfield(s, name)
        val = s.(name);
    else
        val = defaultVal;
    end
end

%--------------------------------------------------------------------------
function progressRegistry(action, dlg)
    persistent registry;
    if isempty(registry)
        registry = {};
    end

    if nargin < 1 || isempty(action)
        action = 'list';
    end

    registry = pruneInvalidProgress(registry);

    switch lower(string(action))
        case "add"
            if nargin >= 2 && ~isempty(dlg)
                registry{end+1} = dlg; %#ok<AGROW>
            end
        case "remove"
            if nargin >= 2 && ~isempty(dlg)
                keep = true(size(registry));
                for ii = 1:numel(registry)
                    keep(ii) = ~sameProgressHandle(registry{ii}, dlg);
                end
                registry = registry(keep);
            end
        case "closeall"
            for ii = 1:numel(registry)
                try
                    reallyCloseProgressHandle(registry{ii});
                catch
                end
            end
            registry = {};
        case "list"
            % no-op; allows inspection for debugging
    end
end

%--------------------------------------------------------------------------
function cleanList = pruneInvalidProgress(listIn)
    if isempty(listIn)
        cleanList = listIn;
        return;
    end
    keep = true(size(listIn));
    for ii = 1:numel(listIn)
        dlg = listIn{ii};
        if isempty(dlg)
            keep(ii) = false;
            continue;
        end
        if ~isfield(dlg, 'h')
            keep(ii) = false;
            continue;
        end
        try
            switch dlg.type
                case 'uiprogress'
                    keep(ii) = isvalid(dlg.h);
                case 'waitbar'
                    keep(ii) = ishghandle(dlg.h);
            end
        catch
            keep(ii) = false;
        end
    end
    cleanList = listIn(keep);
end

%--------------------------------------------------------------------------
function tf = sameProgressHandle(a, b)
    tf = false;
    if isempty(a) || isempty(b)
        return;
    end
    if ~isfield(a, 'type') || ~isfield(b, 'type') || ~strcmp(a.type, b.type)
        return;
    end
    if ~isfield(a, 'h') || ~isfield(b, 'h')
        return;
    end
    try
        switch a.type
            case 'uiprogress'
                tf = isequal(a.h, b.h);
            case 'waitbar'
                tf = isequal(double(a.h), double(b.h));
        end
    catch
        tf = false;
    end
end

%--------------------------------------------------------------------------
function forceCloseExportDialogs()
    progressRegistry('closeall');
    try
        tagged = findall(groot, 'Tag','RMBV_EXPORT_PROGRESS');
        for ii = 1:numel(tagged)
            try
                delete(tagged(ii));
            catch
            end
        end

        % Close any known waitbar or progress figures tagged/named for export
        lingering = findall(groot, 'Type','figure');
        for ii = 1:numel(lingering)
            fig = lingering(ii);
            try
                hasName = isprop(fig, 'Name') && strcmp(get(fig,'Name'), 'Export');
                hasTitle = isprop(fig, 'Title') && strcmp(get(fig,'Title'), 'Export');
                isWaitbarTag = isprop(fig, 'Tag') && strcmp(get(fig,'Tag'), 'TMWWaitbar');
                hasExportTag = isprop(fig, 'Tag') && strcmp(get(fig,'Tag'), 'RMBV_EXPORT_PROGRESS');
                if hasName || hasTitle || isWaitbarTag || hasExportTag
                    try
                        close(fig);
                    catch
                        delete(fig);
                    end
                end
            catch
            end
        end

        % Close any uiprogress dialogs that may not present as figures
        try
            prgs = findall(groot, '-isa', 'matlab.ui.dialog.ProgressDialog');
            for ii = 1:numel(prgs)
                dlg = prgs(ii);
                try
                    close(dlg);
                catch
                    delete(dlg);
                end
            end
        catch
        end

        drawnow;
    catch
    end
end

%--------------------------------------------------------------------------
function fname = makeCacheFieldName(key)
    fname = regexprep(key, '[^A-Za-z0-9_]', '_');
    if isempty(fname)
        fname = 'pane';
    end
    if ~isletter(fname(1)) && fname(1) ~= '_'
        fname = ['pane_' fname];
    end
end

%--------------------------------------------------------------------------
function keyOut = keyify(keyIn)
    if isstring(keyIn)
        if numel(keyIn) >= 1
            keyOut = char(keyIn(1));
        else
            keyOut = '';
        end
    elseif iscellstr(keyIn) || (iscell(keyIn) && numel(keyIn)==1 && ischar(keyIn{1}))
        keyOut = keyIn{1};
    else
        keyOut = keyIn;
    end
end

%--------------------------------------------------------------------------
function key = hsiKey(sensor, modality)
    if nargin < 2
        modality = '';
    end
    if isempty(modality)
        key = sprintf('HSI-%s', upper(sensor));
    else
        key = sprintf('HSI-%s-%s', upper(sensor), upper(modality));
    end
end

%--------------------------------------------------------------------------
function v = getOr(mapObj, k, defaultVal)
    if nargin < 3
        defaultVal = [];
    end
    k = keyify(k);
    try
        if isa(mapObj, 'containers.Map') && isKey(mapObj, k)
            v = mapObj(k);
        else
            v = defaultVal;
        end
    catch
        v = defaultVal;
    end
end

%--------------------------------------------------------------------------
function evtKey = hsiCacheKey(evt)
    evtKey = 'none';
    if isempty(evt) || ~isstruct(evt)
        return;
    end
    parts = {};
    if isfield(evt,'sensor'), parts{end+1} = evt.sensor; end
    if isfield(evt,'path'), parts{end+1} = evt.path; end
    if isfield(evt,'modality'), parts{end+1} = evt.modality; end
    if isfield(evt,'groupIdx'), parts{end+1} = sprintf('g%d', evt.groupIdx); end
    if isfield(evt,'scanIdx'), parts{end+1} = sprintf('s%d', evt.scanIdx); end
    evtKey = strjoin(parts, '|');
end

%--------------------------------------------------------------------------
function val = getArg(args, idx, defaultVal)
    if nargin < 3
        defaultVal = [];
    end
    if numel(args) >= idx
        val = args{idx};
    else
        val = defaultVal;
    end
end

%--------------------------------------------------------------------------
function optsOut = applyDefaults(optsIn, defaults)
    optsOut = defaults;
    fields = fieldnames(optsIn);
    for ii = 1:numel(fields)
        optsOut.(fields{ii}) = optsIn.(fields{ii});
    end
end

%--------------------------------------------------------------------------
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
