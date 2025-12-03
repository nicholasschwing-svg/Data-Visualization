function [cerbTimesByDay, cerbMetaByDay] = scanCerberusFiles( ...
    dataRoot, dateList, cerbPattern, fnameTimePattern)
% scanCerberusFiles
% Recursively scans dataRoot for CERBERUS files matching cerbPattern,
% parses timestamps from their filenames using fnameTimePattern,
% and buckets them into per-day cell arrays.
%
% Inputs:
%   dataRoot         - root directory to scan
%   dateList         - Nx1 datetime array of days of interest
%   cerbPattern      - filename pattern for CERBERUS files (e.g. '*cal_hsi*')
%   fnameTimePattern - regexp to extract timestamp from filename
%
% Outputs:
%   cerbTimesByDay   - cell array, one cell per dateList entry with datetime vector
%   cerbMetaByDay    - cell array, one cell per dateList entry with struct:
%                      .time  (datetime vector)
%                      .paths (cell array of full file paths)

    nDays = numel(dateList);
    cerbTimesByDay = cell(nDays, 1);
    cerbMetaByDay  = cell(nDays, 1);

    % Initialize outputs
    for k = 1:nDays
        cerbTimesByDay{k} = datetime.empty(0,1);
        cerbMetaByDay{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
    end

    if nargin < 4
        error('scanCerberusFiles:NotEnoughInputs', ...
            'Usage: scanCerberusFiles(dataRoot, dateList, cerbPattern, fnameTimePattern)');
    end

    if isempty(dataRoot) || ~isfolder(dataRoot)
        warning('scanCerberusFiles:InvalidRoot', ...
            'Data root "%s" is not a valid folder.', dataRoot);
        return;
    end

    % Recursively find CERBERUS data files
    %
    % NOTE: '**' wildcard requires R2016b or newer. If your MATLAB
    % is older, you can replace this with a genpath-based recursive dir.
    searchPattern = fullfile(dataRoot, '**', cerbPattern);
    d = dir(searchPattern);

    if isempty(d)
        fprintf('scanCerberusFiles: No files found with pattern "%s" under "%s".\n', ...
            cerbPattern, dataRoot);
        return;
    end

    allTimes = datetime.empty(0,1);
    allPaths = {};
    nSkipped = 0;

    for i = 1:numel(d)
        fname = d(i).name;
        fpath = fullfile(d(i).folder, fname);

        % Parse time from filename
        dt = parseCerbFilenameTime(fname, fnameTimePattern);
        if isnat(dt)
            nSkipped = nSkipped + 1;
            continue;
        end

        allTimes(end+1,1) = dt; %#ok<AGROW>
        allPaths{end+1,1} = fpath; %#ok<AGROW>
    end

    fprintf('\nscanCerberusFiles: processed %d files, skipped %d (no timestamp match).\n', ...
        numel(allTimes), nSkipped);

    if isempty(allTimes)
        % Nothing parsed successfully
        return;
    end

    % Bucket into days specified by dateList
    for di = 1:nDays
        dayStart = dateshift(dateList(di), 'start', 'day');
        dayEnd   = dateshift(dateList(di), 'end',   'day');

        idx = allTimes >= dayStart & allTimes <= dayEnd;

        cerbTimesByDay{di} = allTimes(idx);

        % Store metadata: one struct per day with times + paths
        cerbMetaByDay{di} = struct( ...
            'time',  allTimes(idx), ...
            'paths', {allPaths(idx)} );
    end
end
