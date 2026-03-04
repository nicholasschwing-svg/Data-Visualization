function [sampleTs, idx, deltaMs] = mv_get_nearest_sample(times, targetTs)
%MV_GET_NEAREST_SAMPLE Find nearest timestamp with binary search.
    sampleTs = [];
    idx = NaN;
    deltaMs = NaN;
    if isempty(times) || isempty(targetTs)
        return;
    end

    times = times(:);
    if isdatetime(times)
        times = sort(times(~isnat(times)));
        if isempty(times) || ~isdatetime(targetTs) || isnat(targetTs)
            return;
        end
        nums = posixtime(times) * 1000;
        tNum = posixtime(targetTs) * 1000;
    else
        nums = sort(times(~isnan(times)));
        tNum = targetTs;
    end

    n = numel(nums);
    if n == 0
        return;
    end
    if tNum <= nums(1)
        idx = 1;
    elseif tNum >= nums(end)
        idx = n;
    else
        hi = find(nums >= tNum, 1, 'first');
        lo = hi - 1;
        if abs(nums(hi) - tNum) < abs(nums(lo) - tNum)
            idx = hi;
        else
            idx = lo;
        end
    end

    if isdatetime(times)
        sampleTs = times(idx);
    else
        sampleTs = nums(idx);
    end
    deltaMs = nums(idx) - tNum;
end
