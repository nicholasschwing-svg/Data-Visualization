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

    function kOut = keyify(kIn)
        if isstring(kIn) && isscalar(kIn)
            kOut = char(kIn);
        else
            kOut = kIn;
        end
    end

    function mapOut = normalizeMapKeys(mapIn)
        if isempty(mapIn) || ~isa(mapIn,'containers.Map')
            mapOut = mapIn;
            return;
        end
        mapOut = containers.Map('KeyType','char','ValueType','any');
        keysIn = mapIn.keys;
        for jj = 1:numel(keysIn)
            k = keyify(keysIn{jj});
            mapOut(k) = mapIn(keysIn{jj});
        end
    end

    function v = getOr(mapObj, k, defaultVal)
        if nargin < 3
            defaultVal = [];
        end
        k = keyify(k);
        if isa(mapObj,'containers.Map') && isKey(mapObj, k)
            v = mapObj(k);
        else
            v = defaultVal;
        end
    end

    hdrsMap   = normalizeMapKeys(hdrsMap);
    existsMap = normalizeMapKeys(existsMap);

    % Preserve the historical modality priority so "dt" still matches the
    % first populated band-name list.
    tryMods = {'LWIR','MWIR','SWIR','MONO','VIS-COLOR'};
    for ii = 1:numel(tryMods)
        m = keyify(tryMods{ii});
        dtMap(m) = datetime.empty(0,1);

        if ~getOr(existsMap, m, false)
            continue;
        end

        hdr = getOr(hdrsMap, m, []);

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
