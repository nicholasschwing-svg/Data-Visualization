function [dt, dtMap] = fridge_derive_times_from_hdrs(hdrsMap, existsMap)
% FRIDGE_DERIVE_TIMES_FROM_HDRS
%   Build per-modality datetime vectors from ENVI "band names" when they
%   are strings like "yyyy-MM-dd HH:mm:ss.SSS".
%
% hdrsMap, existsMap: containers.Map keyed by modality.
%
% Outputs:
%   dt     - the first non-empty modality time vector (for backward compat)
%   dtMap  - containers.Map(modality -> datetime vector or [])

    dt    = [];
    dtMap = containers.Map('KeyType','char','ValueType','any');

    % Preserve the historical modality priority so "dt" still matches the
    % first populated band-name list.
    tryMods = {'LWIR','MWIR','SWIR','MONO','VIS-COLOR'};
    for ii = 1:numel(tryMods)
        m = tryMods{ii};
        if isstring(m) && isscalar(m)
            m = char(m);
        end
        dtMap(m) = datetime.empty(0,1);

        if ~existsMap(m)
            continue;
        end

        hdr = hdrsMap(m);

        if isfield(hdr,'band_names')
            rawNames = hdr.band_names;
        elseif isfield(hdr,'bandNames')
            rawNames = hdr.bandNames;
        else
            continue;
        end

        if ischar(rawNames)
            parts = regexp(rawNames, ',', 'split');
            parts = strtrim(parts);
        else
            parts = rawNames;
        end

        parts = parts(~cellfun(@isempty,parts));
        if isempty(parts)
            continue;
        end

        dtCandidate = datetime.empty(0,1);
        try
            dtCandidate = datetime(parts, ...
                'InputFormat','yyyy-MM-dd HH:mm:ss.SSSSSS');
        catch
            try
                dtCandidate = datetime(parts, ...
                    'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            catch
                dtCandidate = datetime.empty(0,1);
            end
        end

        if ~isempty(dtCandidate)
            dtMap(m) = dtCandidate(:);
            if isempty(dt)
                dt = dtMap(m);
            end
        end
    end
end
