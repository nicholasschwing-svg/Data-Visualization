function fridgeSel = selectFridgeInstanceInRange(xMin, xMax, insts)
% selectFridgeInstanceInRange
%   Choose the earliest FRIDGE instance whose time span overlaps the
%   given range [xMin,xMax] in hours-of-day.
%
% insts is the per-day struct array used in the timeline:
%   .startTime, .endTime, .path, .wavelength

    fridgeSel = struct('has', false, 'instance', []);

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

    % Earliest by startTime
    startTimes = [insts(idxOverlap).startTime];
    [~, idxLocal] = min(startTimes);
    idxFridge = idxOverlap(idxLocal);

    fridgeSel.has      = true;
    fridgeSel.instance = insts(idxFridge);
end
