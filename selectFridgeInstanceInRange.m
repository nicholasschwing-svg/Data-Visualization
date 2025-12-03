function fridgeSel = selectFridgeInstanceInRange(xMin, xMax, insts, anchorTime)
% selectFridgeInstanceInRange
%   Choose a FRIDGE instance whose time span overlaps the given range
%   [xMin,xMax] in hours-of-day. If an anchorTime (datetime) is provided,
%   prefer the instance whose interval is closest to that anchor (zero
%   distance when the anchor falls inside the interval). Otherwise, fall
%   back to the earliest overlapping instance.
%
% insts is the per-day struct array used in the timeline:
%   .startTime, .endTime, .path, .wavelength

    fridgeSel = struct('has', false, 'instance', []);

    if nargin < 4
        anchorTime = [];
    end

    if isempty(insts) || isempty(insts(1).startTime)
        return;
    end

    nInst = numel(insts);
    xs    = zeros(nInst,1);
    xe    = zeros(nInst,1);

    for kk = 1:nInst
        t1 = insts(kk).startTime;
        t2 = insts(kk).endTime;
        xs(kk) = hour(t1) + minute(t1)/60 + second(t1)/3600;
        xe(kk) = hour(t2) + minute(t2)/60 + second(t2)/3600;
    end

    % Overlap condition: [xs, xe] intersects [xMin, xMax]
    overlaps = (xe >= xMin) & (xs <= xMax);

    idxOverlap = find(overlaps);
    if isempty(idxOverlap)
        return;
    end

    if ~isempty(anchorTime) && isdatetime(anchorTime) && ~isnat(anchorTime)
        % Measure distance from anchor to each overlapping interval.
        % If the anchor is inside the interval, distance is zero.
        dist = zeros(numel(idxOverlap),1);
        for ii = 1:numel(idxOverlap)
            k = idxOverlap(ii);
            tStart = insts(k).startTime;
            tEnd   = insts(k).endTime;
            if anchorTime >= tStart && anchorTime <= tEnd
                dist(ii) = 0;
            else
                dist(ii) = min(abs(anchorTime - tStart), abs(anchorTime - tEnd));
            end
        end
        [~, idxLocal] = min(dist);
    else
        % Earliest by startTime
        startTimes = [insts(idxOverlap).startTime];
        [~, idxLocal] = min(startTimes);
    end

    idxFridge = idxOverlap(idxLocal);

    fridgeSel.has      = true;
    fridgeSel.instance = insts(idxFridge);
end
