function [tStart, tEnd, waveLabel, captureKey] = parseFridgeHeader(headerPath)
% parseFridgeHeader
% Extracts:
%   - tStart     : datetime (earliest frame timestamp)
%   - tEnd       : datetime (latest frame timestamp)
%   - waveLabel  : 'LWIR','MWIR','SWIR','MONO','VIS', or 'UNKNOWN'
%   - captureKey : string identifying the capture without wavelength
%                  e.g. AARO_core_7 for AARO_core_7_LWIR.hdr
%
% Assumes FRIDGE headers look like ENVI headers with:
%   band names = {
%   2024-11-18 14:22:18.390,
%   2024-11-18 14:22:18.406,
%   ...
%   }
%
% and that wavelength is encoded in the filename, e.g.:
%   AARO_core_7_LWIR.hdr

    tStart     = NaT;
    tEnd       = NaT;
    waveLabel  = 'UNKNOWN';
    captureKey = "";

    if ~isfile(headerPath)
        return;
    end

    try
        txt = fileread(headerPath);
    catch
        % Could not read file
        return;
    end

    % ---------- PARSE ALL BAND TIMESTAMPS ----------
    %
    % Match patterns like: 2024-11-18 14:22:18.390
    timeMatches = regexp(txt, ...
        '(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3})', ...
        'match');

    if ~isempty(timeMatches)
        try
            dts = datetime(timeMatches, ...
                'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
            tStart = min(dts);
            tEnd   = max(dts);
        catch
            % If parsing fails, leave tStart/tEnd as NaT
        end
    end

    % ---------- DERIVE WAVELENGTH & CAPTURE KEY FROM FILENAME ----------

    [~, baseName, ~] = fileparts(headerPath);   % e.g. 'AARO_core_7_LWIR'

    % Look for pattern: <anything>_(LWIR|MWIR|SWIR|MONO|VIS)
    m = regexp(baseName, ...
        '^(.*)_(LWIR|MWIR|SWIR|MONO|VIS)$', ...
        'tokens', 'once');

    if ~isempty(m)
        captureKey = string(m{1});   % 'AARO_core_7'
        waveLabel  = upper(m{2});    % 'LWIR', etc.
    else
        % Fallback: try to find wavelength token anywhere in the name
        waves = {'LWIR','MWIR','SWIR','MONO','VIS'};
        for w = 1:numel(waves)
            if contains(upper(baseName), waves{w})
                waveLabel = waves{w};
                % Remove "_<wave>" from end if present for captureKey
                captureKey = regexprep(baseName, ['_' waves{w} '$'], '');
                captureKey = string(captureKey);
                break;
            end
        end

        if captureKey == ""
            % As a last resort, just use the basename as capture key
            captureKey = string(baseName);
        end
    end
end
