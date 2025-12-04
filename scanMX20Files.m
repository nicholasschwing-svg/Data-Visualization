function [mxTimesByDay, mxMetaByDay, dateListOut] = scanMX20Files( ...
    dataRoot, dateList, fnameTimePattern)
% scanMX20Files
% Recursively scans dataRoot for MX20 HSI headers (.hdr), parses timestamps
% from their filenames using fnameTimePattern, and buckets them into per-day
% cell arrays. When the caller provides an empty dateList, the set of dates
% is derived from the data that was found.
%
% Inputs:
%   dataRoot         - HSI root (or top-level root containing HSI)
%   dateList         - Nx1 datetime array of days of interest (can be empty)
%   fnameTimePattern - regexp to extract timestamp from filename
%
% Outputs:
%   mxTimesByDay     - cell array, one cell per dateList entry with datetime vector
%   mxMetaByDay      - cell array, one cell per dateList entry with struct:
%                      .time  (datetime vector)
%                      .paths (cell array of full file paths)
%   dateListOut      - date list actually used for bucketing (derived when
%                      input dateList is empty)

    if nargin < 2 || isempty(dateList)
        dateList = datetime.empty(0,1);
    end

    dateListOut  = dateList;
    nDaysOut     = numel(dateListOut);
    mxTimesByDay = cell(nDaysOut,1);
    mxMetaByDay  = cell(nDaysOut,1);
    for k = 1:nDaysOut
        mxTimesByDay{k} = datetime.empty(0,1);
        mxMetaByDay{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
    end

    if isempty(dataRoot) || ~isfolder(dataRoot)
        warning('scanMX20Files:InvalidRoot', ...
            'Data root "%s" is not a valid folder.', dataRoot);
        return;
    end

    % For MX20 we use header files (*.hdr) as the canonical handle. Search
    % recursively in case the folder layout differs between datasets.
    searchPattern = fullfile(dataRoot, '**', '*.hdr');
    d = dir(searchPattern);

    if isempty(d)
        return;
    end

    allTimes = datetime.empty(0,1);
    allPaths = {};
    nSkipped = 0;

    for i = 1:numel(d)
        fname = d(i).name;
        fpath = fullfile(d(i).folder, fname);

        dt = parseCerbFilenameTime(fname, fnameTimePattern);  % reuse helper
        if isnat(dt)
            nSkipped = nSkipped + 1;
            continue;
        end

        allTimes(end+1,1) = dt; %#ok<AGROW>
        allPaths{end+1,1} = fpath; %#ok<AGROW>
    end

    fprintf('scanMX20Files: processed %d files, skipped %d (no timestamp match).\n', ...
        numel(allTimes), nSkipped);

    if isempty(allTimes)
        return;
    end

    if isempty(dateListOut)
        dateListOut = unique(dateshift(allTimes,'start','day'));
    end

    nDays = numel(dateListOut);
    mxTimesByDay = cell(nDays, 1);
    mxMetaByDay  = cell(nDays, 1);
    for k = 1:nDays
        mxTimesByDay{k} = datetime.empty(0,1);
        mxMetaByDay{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
    end

    for di = 1:nDays
        dayStart = dateshift(dateListOut(di), 'start', 'day');
        dayEnd   = dateshift(dateListOut(di), 'end',   'day');

        idx = allTimes >= dayStart & allTimes <= dayEnd;
        mxTimesByDay{di} = allTimes(idx);
        mxMetaByDay{di}  = struct('time',  allTimes(idx), ...
                                  'paths', {allPaths(idx)});
    end
end
