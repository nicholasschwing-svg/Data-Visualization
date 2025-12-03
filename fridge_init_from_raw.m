function [filesMap, hdrsMap, existsMap, maxFramesMap, nFrames, fridgeTimes, fridgeTimesMap] = ...
    fridge_init_from_raw(pathStr, prefix, modalities)
% FRIDGE_INIT_FROM_RAW
%   Scan all FRIDGE modalities for a given prefix and directory.
%
% Inputs:
%   pathStr   - folder containing the RAW/HDR files
%   prefix    - filename prefix (before "_LWIR.raw", etc.)
%   modalities- cell array of modality names (e.g. {'LWIR','MWIR',...})
%
% Outputs:
%   filesMap     - containers.Map of modality -> raw file path
%   hdrsMap      - containers.Map of modality -> header struct
%   existsMap    - containers.Map of modality -> logical
%   maxFramesMap - containers.Map of modality -> max usable frames
%   nFrames      - overall max frames across modalities (>=1)
%   fridgeTimes  - datetime vector from ENVI "band names" (or [] if none)
%   fridgeTimesMap - containers.Map of modality -> datetime vector (or [])

    filesMap     = containers.Map(modalities, repmat({''},1,numel(modalities)));
    hdrsMap      = containers.Map(modalities, repmat({[]},1,numel(modalities)));
    existsMap    = containers.Map(modalities, num2cell(false(1,numel(modalities))));
    maxFramesMap = containers.Map(modalities, num2cell(nan(1,numel(modalities))));

    frameCounts = nan(1,numel(modalities));

    for i = 1:numel(modalities)
        m = modalities{i};

        rawPath = fullfile(pathStr, sprintf('%s_%s.raw', prefix, m));
        filesMap(m) = rawPath;

        hdrPath = strrep(rawPath, '.raw', '.hdr');
        existsMap(m) = isfile(rawPath) && isfile(hdrPath);

        if ~existsMap(m)
            maxFramesMap(m) = NaN;
            continue;
        end

        hdr = readENVIHeader(hdrPath);
        % Inject band names from the raw header text so per-frame timestamps
        % are available even though readENVIHeader only returns core fields.
        hdr = injectBandNamesFromText(hdrPath, hdr);
        hdrsMap(m) = hdr;

        % VIS: possibly 3 bands per frame, others 1 band per frame
        if strcmp(m,'VIS-COLOR') && isfield(hdr,'bands') && ...
                hdr.bands >= 3 && mod(double(hdr.bands),3) == 0
            hdrFrames = floor(double(hdr.bands)/3);
        else
            hdrFrames = hdr.bands;
        end

        fsFrames = framesFromFile(m, hdr, rawPath);   % your existing helper

        cand = [hdrFrames, fsFrames];
        cand = cand(~isnan(cand) & cand > 0);
        if isempty(cand)
            frameCounts(i) = NaN;
        else
            frameCounts(i) = min(cand);
        end
        maxFramesMap(m) = frameCounts(i);
    end

    valid = frameCounts(~isnan(frameCounts) & frameCounts > 0);
    if isempty(valid)
        nFrames = 1;
    else
        nFrames = max(1, floor(max(valid)));
    end

    [fridgeTimes, fridgeTimesMap] = fridge_derive_times_from_hdrs(hdrsMap, existsMap);

end

function hdrOut = injectBandNamesFromText(hdrPath, hdrIn)
% injectBandNamesFromText
%   readENVIHeader ignores band_names, but the FRIDGE headers carry per-frame
%   timestamps inside that field. Parse them directly from the header text so
%   downstream time alignment can use them.

    hdrOut = hdrIn;

    if isfield(hdrIn, 'band_names') && ~isempty(hdrIn.band_names)
        return;
    end

    try
        txt = fileread(hdrPath);
    catch
        hdrOut.band_names = {};
        return;
    end

    matches = regexp(txt, '(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3})', 'match');
    if isempty(matches)
        hdrOut.band_names = {};
    else
        hdrOut.band_names = matches(:)';
    end

    % Alias sometimes used by other loaders
    hdrOut.bandNames = hdrOut.band_names;
end
