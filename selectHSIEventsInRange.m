function [cerbSel, mxSel, fastSel] = selectHSIEventsInRange( ...
    xMin, xMax, yMin, yMax, ...
    cerbUD, mxUD, fastUDMap, cerbY, mxY, fastY, ...
    cerbEnabled, mxEnabled, fastEnabledMap)
% selectHSIEventsInRange
%   Given a horizontal time range [xMin,xMax] in hours-of-day and a
%   vertical range [yMin,yMax] in axis units, choose which HSI sensor
%   (CERBERUS, MX20, FAST, or multiple) should be considered and return
%   the earliest event in that range for each.
%
% Inputs:
%   cerbUD.meta.paths, cerbUD.times - as stored in CERB scatter UserData
%   mxUD.meta.paths,   mxUD.times   - as stored in MX20 scatter UserData
%   fastUDMap                      - containers.Map of modality->UserData
%   fastY                          - y-position for FAST markers
%   fastEnabledMap                 - containers.Map of modality->logical
%
% Outputs:
%   cerbSel, mxSel, fastSel - structs with fields:
%       .has  (logical)
%       .time (datetime or [])
%       .path (char or '')
%       .modality (FAST only)

    % Default outputs
    cerbSel = struct('has', false, 'time', [], 'path', '', ...
                     'timesInRange', datetime.empty(0,1), ...
                     'pathsInRange', {{}});
    mxSel   = struct('has', false, 'time', [], 'path', '', ...
                     'timesInRange', datetime.empty(0,1), ...
                     'pathsInRange', {{}});
    fastSel = struct('has', false, 'time', [], 'path', '', 'modality', '', ...
                     'timesInRange', datetime.empty(0,1), ...
                     'pathsInRange', {{}}, 'modalitiesInRange', {{}});

    % Decide which rows we target based on vertical coverage of the
    % rectangle. Any box that overlaps either HSI row should be considered
    % valid; if it overlaps neither row, default to allowing both.
    hsiThresh = 0.10;
    overlapsCerb = (yMax >= (cerbY - hsiThresh)) && (yMin <= (cerbY + hsiThresh));
    overlapsMx   = (yMax >= (mxY   - hsiThresh)) && (yMin <= (mxY   + hsiThresh));
    overlapsFast = (yMax >= (fastY - hsiThresh)) && (yMin <= (fastY + hsiThresh));

    selectCerb = overlapsCerb || (~overlapsCerb && ~overlapsMx && ~overlapsFast);
    selectMx   = overlapsMx   || (~overlapsCerb && ~overlapsMx && ~overlapsFast);
    selectFast = overlapsFast || (~overlapsCerb && ~overlapsMx && ~overlapsFast);

    % Apply checkbox enables
    if ~cerbEnabled
        selectCerb = false;
    end
    if ~mxEnabled
        selectMx = false;
    end
    if nargin < 11 || isempty(fastEnabledMap)
        selectFast = false;
    end
    % Any modality must be enabled to allow FAST selection
    if selectFast && isa(fastEnabledMap, 'containers.Map')
        enabledVals = fastEnabledMap.values;
        if ~any(cellfun(@(v) logical(v), enabledVals))
            selectFast = false;
        end
    end

    % ----- CERBERUS -----
    if selectCerb && ~isempty(cerbUD) && isfield(cerbUD,'times') && ~isempty(cerbUD.times)
        timesToday = cerbUD.times;
        metaToday  = cerbUD.meta;

        h = hour(timesToday) + minute(timesToday)/60 + second(timesToday)/3600;
        inRange = (h >= xMin) & (h <= xMax);

        idxCandidates = find(inRange);
        cerbSel.timesInRange = timesToday(idxCandidates);
        cerbSel.pathsInRange = metaToday.paths(idxCandidates);

        if ~isempty(idxCandidates)
            % Earliest in time
            [~, idxLocal] = min(timesToday(idxCandidates));  % datetime min
            idxCerb = idxCandidates(idxLocal);

            cerbSel.has  = true;
            cerbSel.time = timesToday(idxCerb);
            cerbSel.path = metaToday.paths{idxCerb};
        end
    end

    % ----- MX20 -----
    if selectMx && ~isempty(mxUD) && isfield(mxUD,'times') && ~isempty(mxUD.times)
        mxTimesToday = mxUD.times;
        mxMetaToday  = mxUD.meta;

        hm = hour(mxTimesToday) + minute(mxTimesToday)/60 + second(mxTimesToday)/3600;
        inRangeMx = (hm >= xMin) & (hm <= xMax);

        idxMxCandidates = find(inRangeMx);
        mxSel.timesInRange = mxTimesToday(idxMxCandidates);
        mxSel.pathsInRange = mxMetaToday.paths(idxMxCandidates);

        if ~isempty(idxMxCandidates)
            [~, idxLocal] = min(mxTimesToday(idxMxCandidates));  % earliest
            idxMx = idxMxCandidates(idxLocal);

            mxSel.has  = true;
            mxSel.time = mxTimesToday(idxMx);
            mxSel.path = mxMetaToday.paths{idxMx};
        end
    end

    % ----- FAST (aggregate across modalities) -----
    if selectFast && isa(fastUDMap, 'containers.Map') && ~isempty(fastUDMap)
        fastKeys = fastUDMap.keys;
        allTimes = datetime.empty(0,1);
        allPaths = {};
        allMods  = {};

        for ii = 1:numel(fastKeys)
            key = fastKeys{ii};
            if isKey(fastEnabledMap, key) && ~fastEnabledMap(key)
                continue;
            end
            ud = fastUDMap(key);
            if isempty(ud) || ~isfield(ud,'times') || isempty(ud.times)
                continue;
            end
            tVec = ud.times;
            h = hour(tVec) + minute(tVec)/60 + second(tVec)/3600;
            inRange = (h >= xMin) & (h <= xMax);
            idxCandidates = find(inRange);
            fastSel.timesInRange = [fastSel.timesInRange; tVec(idxCandidates)]; %#ok<AGROW>
            fastSel.pathsInRange = [fastSel.pathsInRange; ud.meta.paths(idxCandidates)]; %#ok<AGROW>
            fastSel.modalitiesInRange = [fastSel.modalitiesInRange; repmat({key}, numel(idxCandidates), 1)]; %#ok<AGROW>

            allTimes = [allTimes; tVec(idxCandidates)]; %#ok<AGROW>
            allPaths = [allPaths; ud.meta.paths(idxCandidates)]; %#ok<AGROW>
            allMods  = [allMods; repmat({key}, numel(idxCandidates), 1)]; %#ok<AGROW>
        end

        if ~isempty(allTimes)
            [~, idxLocal] = min(allTimes);
            fastSel.has      = true;
            fastSel.time     = allTimes(idxLocal);
            fastSel.path     = allPaths{idxLocal};
            fastSel.modality = allMods{idxLocal};
        end
    end
end
