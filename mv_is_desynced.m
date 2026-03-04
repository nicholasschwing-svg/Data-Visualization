function tf = mv_is_desynced(playheadTs, panelTimes, toleranceMs)
%MV_IS_DESYNCED True when any panel differs from playhead by tolerance.
    if nargin < 3 || isempty(toleranceMs)
        toleranceMs = 250;
    end
    tf = false;
    if isempty(playheadTs) || isempty(panelTimes)
        return;
    end
    for i = 1:numel(panelTimes)
        t = panelTimes{i};
        if isempty(t)
            continue;
        end
        if isdatetime(t)
            d = abs(milliseconds(t - playheadTs));
        else
            d = abs(t - playheadTs);
        end
        if d > toleranceMs
            tf = true;
            return;
        end
    end
end
