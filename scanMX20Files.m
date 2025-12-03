function [mxTimesByDay, mxMetaByDay] = scanMX20Files( ...
    dataRoot, dateList, fnameTimePattern)
% scanMX20Files
% Scans ROOT/HSI/MX20/<mm-dd> for MX20 HSI files,
% parses timestamps from filenames using fnameTimePattern,
% and buckets them into per-day cell arrays.
%
% Inputs:
%   dataRoot         - top-level root (contains FRIDGE and HSI)
%   dateList         - Nx1 datetime array of days of interest
%   fnameTimePattern - regexp to extract timestamp from filename
%
% Outputs:
%   mxTimesByDay     - cell array, one cell per dateList entry with datetime vector
%   mxMetaByDay      - cell array, one cell per dateList entry with struct:
%                      .time  (datetime vector)
%                      .paths (cell array of full file paths)

    nDays = numel(dateList);
    mxTimesByDay = cell(nDays, 1);
    mxMetaByDay  = cell(nDays, 1);

    for k = 1:nDays
        mxTimesByDay{k} = datetime.empty(0,1);
        mxMetaByDay{k}  = struct('time', datetime.empty(0,1), 'paths', {{}}); %#ok<CCAT>
    end

    if isempty(dataRoot) || ~isfolder(dataRoot)
        warning('scanMX20Files:InvalidRoot', ...
            'Data root "%s" is not a valid folder.', dataRoot);
        return;
    end

    for di = 1:nDays
        dayStr = datestr(dateList(di), 'mm-dd');  % '11-18', etc.

        mxDayRoot = fullfile(dataRoot, 'HSI', 'MX20', dayStr);
        if ~isfolder(mxDayRoot)
            % No MX20 for this day
            continue;
        end

        % For MX20 we use the header files (*.hdr) as the canonical handle
        searchPattern = fullfile(mxDayRoot, '**', '*.hdr');
        d = dir(searchPattern);

        if isempty(d)
            continue;
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

        fprintf('scanMX20Files: day %s: %d files, skipped %d (no timestamp match).\n', ...
            datestr(dateList(di),'mm/dd'), numel(allTimes), nSkipped);

        mxTimesByDay{di} = allTimes;
        mxMetaByDay{di}  = struct('time',  allTimes, ...
                                  'paths', {allPaths});
    end
end
