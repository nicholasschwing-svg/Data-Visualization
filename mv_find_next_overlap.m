function [overlapTs, info] = mv_find_next_overlap(playheadTs, direction, sensorTimesMap, toleranceMs, minSensors, snapMode, masterSensor)
%MV_FIND_NEXT_OVERLAP Find overlap timestamp by scanning candidate times.
    if nargin < 4 || isempty(toleranceMs), toleranceMs = 250; end
    if nargin < 5 || isempty(minSensors), minSensors = 2; end
    if nargin < 6 || isempty(snapMode), snapMode = 'ANY'; end
    if nargin < 7, masterSensor = ''; end

    overlapTs = [];
    info = struct('found', false, 'degraded', false, 'reason', '');
    if isempty(sensorTimesMap) || sensorTimesMap.Count == 0 || isempty(playheadTs)
        info.reason = 'no-data';
        return;
    end

    keys = sensorTimesMap.keys;
    cands = [];
    if strcmp(snapMode,'MASTER') && ~isempty(masterSensor) && isKey(sensorTimesMap, masterSensor)
        cands = sensorTimesMap(masterSensor);
    else
        for i = 1:numel(keys)
            tVec = sensorTimesMap(keys{i});
            cands = [cands; tVec(:)]; %#ok<AGROW>
        end
    end
    if isempty(cands)
        info.reason = 'no-candidates';
        return;
    end
    cands = unique(sort(cands));

    if isdatetime(cands)
        ahead = cands > playheadTs;
        behind = cands < playheadTs;
    else
        ahead = cands > playheadTs;
        behind = cands < playheadTs;
    end
    if direction >= 0
        scan = cands(ahead);
    else
        scan = flipud(cands(behind));
    end

    cap = min(10000, numel(scan));
    bestCount = 0; bestTs = [];
    for i = 1:cap
        t = scan(i);
        count = 0;
        for k = 1:numel(keys)
            res = mv_resolve_panel_sample(sensorTimesMap(keys{k}), t, struct('toleranceMs', toleranceMs));
            if strcmp(res.status, 'OK')
                count = count + 1;
            end
        end
        if count >= minSensors
            overlapTs = t;
            info.found = true;
            return;
        end
        if count > bestCount
            bestCount = count;
            bestTs = t;
        end
    end

    if strcmp(snapMode,'ALL') && ~isempty(bestTs)
        overlapTs = bestTs;
        info.found = true;
        info.degraded = true;
        info.reason = 'degraded-any';
    else
        info.reason = 'not-found';
    end
end
