function out = mv_resolve_panel_sample(sampleTimes, targetTs, opts)
%MV_RESOLVE_PANEL_SAMPLE Resolve panel sample near target playhead.
% out.status: OK | NO_SAMPLE_NEARBY | HOLDING_LAST
    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts,'toleranceMs') || isempty(opts.toleranceMs)
        opts.toleranceMs = 250;
    end
    if ~isfield(opts,'allowHold') || isempty(opts.allowHold)
        opts.allowHold = false;
    end

    [sampleTs, idx, deltaMs] = mv_get_nearest_sample(sampleTimes, targetTs);
    out = struct('sample', sampleTs, 'idx', idx, 'deltaMs', deltaMs, 'status', 'NO_SAMPLE_NEARBY');
    if isempty(sampleTs) || isnan(idx)
        if opts.allowHold
            out.status = 'HOLDING_LAST';
        end
        return;
    end
    if abs(deltaMs) <= opts.toleranceMs
        out.status = 'OK';
    elseif opts.allowHold
        out.status = 'HOLDING_LAST';
    else
        out.status = 'NO_SAMPLE_NEARBY';
    end
end
