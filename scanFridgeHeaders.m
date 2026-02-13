function [fridgeInstancesByDay, dateListOut] = scanFridgeHeaders( ...
    dataRoot, dateList, fridgePattern, defaultDurationSec)
% scanFridgeHeaders
% Recursively scans dataRoot for FRIDGE header files matching fridgePattern,
% parses timestamps & wavelengths using parseFridgeHeader(), and groups
% them into "captures" by captureKey (e.g. AARO_core_7).
%
% For each capture:
%   - Collect all (start,end,wavelength) intervals for that capture
%   - For each wavelength, compute duration = maxEnd - minStart
%   - Pick the wavelength with the LONGEST duration
%   - Create ONE FRIDGE instance for the capture using that wavelength's
%     min start & max end
%
% Then bucket those instances into days for the timeline based on startTime.
%
% Inputs:
%   dataRoot          - root directory to scan
%   dateList          - Nx1 datetime array of days of interest
%   fridgePattern     - filename pattern for FRIDGE headers (e.g. '*FRIDGE*hdr*')
%   defaultDurationSec- duration used when tEnd == NaT (rare, just fallback)
%
% Output:
%   fridgeInstancesByDay - cell array, one cell per dateList entry, each cell
%                          is an array of structs with fields:
%                          .startTime, .endTime, .wavelength, .path
%   dateListOut          - date list actually used for bucketing (derived when
%                          input dateList is empty)

% Lightweight header parse cache avoids repeatedly reparsing unchanged .hdr files.
persistent hdrParseCache;
if isempty(hdrParseCache)
    hdrParseCache = containers.Map('KeyType','char','ValueType','any');
end
perfEnabled = any(strcmp(getenv('RMBV_PERF'), {'1','true','TRUE','on','ON'}));

    if nargin < 2 || isempty(dateList)
        dateList = datetime.empty(0,1);
    end

    dateListOut = dateList;
    nDaysOut    = numel(dateListOut);

    emptyStruct = struct( ...
        'startTime', datetime.empty(0,1), ...
        'endTime',   datetime.empty(0,1), ...
        'wavelength', {{}}, ...
        'path',      {{}} );

    fridgeInstancesByDay = cell(nDaysOut,1);
    for k = 1:nDaysOut
        fridgeInstancesByDay{k} = emptyStruct;
    end

    if nargin < 4
        error('scanFridgeHeaders:NotEnoughInputs', ...
            'Usage: scanFridgeHeaders(dataRoot, dateList, fridgePattern, defaultDurationSec)');
    end

    if isempty(dataRoot) || ~isfolder(dataRoot)
        warning('scanFridgeHeaders:InvalidRoot', ...
            'Data root "%s" is not a valid folder.', dataRoot);
        return;
    end

    % ---------- FIND ALL FRIDGE HEADERS ----------
    searchPattern = fullfile(dataRoot, '**', fridgePattern);
    d = dir(searchPattern);

    if isempty(d)
        fprintf('scanFridgeHeaders: No files found with pattern "%s" under "%s".\n', ...
            fridgePattern, dataRoot);
        return;
    end

    allStartTimes   = datetime.empty(0,1);
    allEndTimes     = datetime.empty(0,1);
    allWaves        = {};
    allPaths        = {};
    allCaptureKeys  = {};
    nSkipped        = 0;

    for i = 1:numel(d)
        fpath = fullfile(d(i).folder, d(i).name);

        cacheKey = sprintf('%s|%d|%.0f', fpath, d(i).bytes, d(i).datenum);
        if isKey(hdrParseCache, cacheKey)
            parsed = hdrParseCache(cacheKey);
            tStart = parsed.tStart;
            tEnd = parsed.tEnd;
            waveLabel = parsed.waveLabel;
            captureKey = parsed.captureKey;
        else
            tParse = tic;
            [tStart, tEnd, waveLabel, captureKey] = parseFridgeHeader(fpath);
            hdrParseCache(cacheKey) = struct('tStart', tStart, 'tEnd', tEnd, ...
                'waveLabel', waveLabel, 'captureKey', captureKey);
            if perfEnabled
                fprintf('[perf] parseFridgeHeader %s: %.1f ms\n', d(i).name, toc(tParse)*1000);
            end
        end

        if isnat(tStart)
            % No usable time info
            nSkipped = nSkipped + 1;
            continue;
        end

        if isnat(tEnd)
            tEnd = tStart + seconds(defaultDurationSec);
        end

        allStartTimes(end+1,1)  = tStart;              %#ok<AGROW>
        allEndTimes(end+1,1)    = tEnd;                %#ok<AGROW>
        allWaves{end+1,1}       = waveLabel;           %#ok<AGROW>
        allPaths{end+1,1}       = fpath;               %#ok<AGROW>
        allCaptureKeys{end+1,1} = char(captureKey);    %#ok<AGROW>
    end

    fprintf('\nscanFridgeHeaders: processed %d headers, skipped %d (no time).\n', ...
        numel(allStartTimes), nSkipped);

    if isempty(allStartTimes)
        return;
    end

    % ---------- GROUP BY CAPTURE KEY ----------
    captureKeysUnique = unique(allCaptureKeys);
    nCaptures = numel(captureKeysUnique);

    % We'll build one "instance" per capture
    captureInstances(nCaptures,1) = struct( ...
        'startTime', datetime, ...
        'endTime',   datetime, ...
        'wavelength', '', ...
        'path',      '' );

    for ci = 1:nCaptures
        key = captureKeysUnique{ci};

        idxKey = strcmp(allCaptureKeys, key);

        keyStarts = allStartTimes(idxKey);
        keyEnds   = allEndTimes(idxKey);
        keyWaves  = allWaves(idxKey);
        keyPaths  = allPaths(idxKey);

        % Per-wavelength duration within this capture
        wavesUnique = unique(keyWaves);
        bestWave      = '';
        bestDur       = duration(0,0,0);
        bestStart     = NaT;
        bestEnd       = NaT;
        bestPath      = '';

        for w = 1:numel(wavesUnique)
            wlab = wavesUnique{w};
            idxW = strcmp(keyWaves, wlab);

            wStarts = keyStarts(idxW);
            wEnds   = keyEnds(idxW);

            % Combined interval for this wavelength in this capture
            wStartMin = min(wStarts);
            wEndMax   = max(wEnds);
            wDur      = wEndMax - wStartMin;

            if wDur > bestDur
                bestDur   = wDur;
                bestWave  = wlab;
                bestStart = wStartMin;
                bestEnd   = wEndMax;
                % Just grab the first path associated with this wavelength
                wPaths = keyPaths(idxW);
                bestPath = wPaths{1};
            end
        end

        % If something went wrong, skip this capture
        if isnat(bestStart) || isnat(bestEnd)
            continue;
        end

        captureInstances(ci).startTime  = bestStart;
        captureInstances(ci).endTime    = bestEnd;
        captureInstances(ci).wavelength = bestWave;
        captureInstances(ci).path       = bestPath;
    end

    % Remove any uninitialized captures (just in case)
    validIdx = ~arrayfun(@(s) isnat(s.startTime), captureInstances);
    captureInstances = captureInstances(validIdx);

    % Derive date list if caller did not provide one
    if isempty(dateListOut)
        dateListOut = unique(dateshift([captureInstances.startTime],'start','day'));
    end

    nDays = numel(dateListOut);
    fridgeInstancesByDay = cell(nDays,1);
    for k = 1:nDays
        fridgeInstancesByDay{k} = emptyStruct;
    end

    % ---------- BUCKET CAPTURE INSTANCES BY DAY ----------
    for di = 1:nDays
        dayStart = dateshift(dateListOut(di), 'start', 'day');
        dayEnd   = dateshift(dateListOut(di), 'end',   'day');

        % Decide which captures belong to this day (by start time)
        instIdx = arrayfun(@(s) s.startTime >= dayStart && s.startTime <= dayEnd, ...
                           captureInstances);

        instToday = captureInstances(instIdx);

        if isempty(instToday)
            fridgeInstancesByDay{di} = emptyStruct;
        else
            fridgeInstancesByDay{di} = instToday;
        end
    end
end
