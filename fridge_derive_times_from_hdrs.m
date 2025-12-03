function dt = fridge_derive_times_from_hdrs(hdrsMap, existsMap)
% FRIDGE_DERIVE_TIMES_FROM_HDRS
%   Build a datetime vector from ENVI "band names" when they are
%   strings like "yyyy-MM-dd HH:mm:ss.SSS".
%
% hdrsMap, existsMap: containers.Map keyed by modality.

    dt = [];

    tryMods = {'LWIR','MWIR','SWIR','MONO','VIS-COLOR'};
    for ii = 1:numel(tryMods)
        m = tryMods{ii};
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

        try
            dtCandidate = datetime(parts, ...
                'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        catch
            continue;
        end

        if ~isempty(dtCandidate)
            dt = dtCandidate(:);
            return;
        end
    end
end
