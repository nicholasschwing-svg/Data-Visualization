function [fastTimesByDay, fastMetaByDay, dateListOut, modalities] = scanFastFiles( ...
    dataRoot, dateList)
% scanFastFiles
%   Recursively scans dataRoot for FAST HSI headers (.hdr) and buckets them
%   by day and modality. Modalities are inferred from filenames using
%   parseFastFilenameTime.
%
% Inputs:
%   dataRoot  - root directory to scan
%   dateList  - Nx1 datetime array of days of interest (can be empty)
%
% Outputs:
%   fastTimesByDay - containers.Map keyed by modality. Each value is a
%                    cell array (one cell per dateList entry) containing
%                    datetime vectors.
%   fastMetaByDay  - containers.Map keyed by modality. Each value is a
%                    cell array (one cell per dateList entry) containing
%                    structs with fields .time and .paths.
%   dateListOut    - date list actually used (derived when input is empty)
%   modalities     - cellstr of modalities discovered (upper-case)

    if nargin < 2 || isempty(dateList)
        dateList = datetime.empty(0,1);
    end

    dateListOut = dateList;
    fastTimesByDay = containers.Map('KeyType','char','ValueType','any');
    fastMetaByDay  = containers.Map('KeyType','char','ValueType','any');
    modalities     = {};

    if isempty(dataRoot) || ~isfolder(dataRoot)
        warning('scanFastFiles:InvalidRoot', ...
            'Data root "%s" is not a valid folder.', dataRoot);
        return;
    end

    % Prefer headers to avoid double-counting .hsic/.hdr pairs.
    searchPattern = fullfile(dataRoot, '**', '*cal_hsi.hdr');
    d = dir(searchPattern);

    if isempty(d)
        fprintf('scanFastFiles: No FAST files found under %s.\n', dataRoot);
        return;
    end

    allTimes      = datetime.empty(0,1);
    allPaths      = {};
    allModalities = {};
    nSkipped      = 0;

    for i = 1:numel(d)
        fname = d(i).name;
        fpath = fullfile(d(i).folder, fname);

        % Require FAST marker somewhere in the path to avoid cross-sensor
        % contamination.
        if ~contains(lower(fpath), 'fast')
            nSkipped = nSkipped + 1;
            continue;
        end

        [dt, mod, pointId] = parseFastFilenameTime(fname); %#ok<ASGLU>
        if isnat(dt) || isempty(mod)
            nSkipped = nSkipped + 1;
            continue;
        end

        allTimes(end+1,1)      = dt; %#ok<AGROW>
        allPaths{end+1,1}      = fpath; %#ok<AGROW>
        allModalities{end+1,1} = mod; %#ok<AGROW>
    end

    fprintf('scanFastFiles: processed %d files, skipped %d (no timestamp/modality match).\n', ...
        numel(allTimes), nSkipped);

    if isempty(allTimes)
        return;
    end

    modalities = unique(allModalities);

    if isempty(dateListOut)
        dateListOut = unique(dateshift(allTimes,'start','day'));
    end

    nDays = numel(dateListOut);

    % Initialize per-modality buckets
    for m = 1:numel(modalities)
        key = modalities{m};

        timesCell = cell(nDays,1);
        metaCell  = cell(nDays,1);
        for k = 1:nDays
            timesCell{k} = datetime.empty(0,1);
            metaCell{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
        end

        fastTimesByDay(key) = timesCell;
        fastMetaByDay(key)  = metaCell;
    end

    % Bucket events per day and modality
    for di = 1:nDays
        dayStart = dateshift(dateListOut(di), 'start', 'day');
        dayEnd   = dateshift(dateListOut(di), 'end',   'day');

        idxDay = allTimes >= dayStart & allTimes <= dayEnd;

        timesDay      = allTimes(idxDay);
        pathsDay      = allPaths(idxDay);
        modsDay       = allModalities(idxDay);

        for m = 1:numel(modalities)
            key = modalities{m};
            inMod = strcmp(modsDay, key);

            timesCell = fastTimesByDay(key);
            metaCell  = fastMetaByDay(key);

            timesCell{di} = timesDay(inMod);
            metaCell{di}  = struct('time',  timesDay(inMod), ...
                                   'paths', {pathsDay(inMod)});

            fastTimesByDay(key) = timesCell;
            fastMetaByDay(key)  = metaCell;
        end
    end
end
