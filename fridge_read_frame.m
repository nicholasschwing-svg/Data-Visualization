function img = fridge_read_frame(modality, frameIdx, hdrsMap, filesMap)
% FRIDGE_READ_FRAME
%   Read a single FRIDGE frame for a given modality and index.
%
% Handles VIS-COLOR (3 bands per frame) and mono-band modalities.

    % Normalize modality keys so callers can pass string scalars without
    % triggering containers.Map key-type errors.
    if isstring(modality) && isscalar(modality)
        modality = char(modality);
    end

    if ~(isa(hdrsMap,'containers.Map') && isa(filesMap,'containers.Map'))
        error('Header/file maps are invalid for modality %s.', modality);
    end

    if ~(isKey(hdrsMap, modality) && isKey(filesMap, modality))
        error('Modality %s not found in provided header/file maps.', modality);
    end

    hdr  = hdrsMap(modality);
    file = filesMap(modality);

    [dtype, bps] = enviDataType(hdr.dataType);

    machine = 'ieee-le';
    if isfield(hdr,'byteOrder') && hdr.byteOrder == 1
        machine = 'ieee-be';
    end

    % VIS: 3 bands/frame; others: 1 band/frame
    if strcmp(modality,'VIS-COLOR') && isfield(hdr,'bands') && ...
            hdr.bands >= 3 && mod(double(hdr.bands),3) == 0

        bandsPerFrame = 3;
        nFramesVis    = floor(double(hdr.bands)/bandsPerFrame);
        if frameIdx > nFramesVis
            error('Requested VIS-COLOR frame %d exceeds available %d.', ...
                  frameIdx, nFramesVis);
        end

        offset = hdr.headerOffset + ...
                 (frameIdx-1) * hdr.samples * hdr.lines * bandsPerFrame * bps;

        fid = fopen(file,'r',machine);
        if fid < 0
            error('Cannot open %s', file);
        end
        cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fseek(fid, offset, 'bof');

        nPix = hdr.samples * hdr.lines * 3;
        dat  = fread(fid, nPix, ['*' dtype]);
        if numel(dat) < nPix
            error('Unexpected EOF when reading VIS-COLOR frame.');
        end
        dat = reshape(dat, [hdr.samples, hdr.lines, 3]);  % x, y, band
        img = permute(dat, [2 1 3]);                      % y, x, band

    else
        offset = hdr.headerOffset + ...
                 (frameIdx-1) * hdr.samples * hdr.lines * bps;

        fid = fopen(file,'r',machine);
        if fid < 0
            error('Cannot open %s', file);
        end
        cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fseek(fid, offset, 'bof');

        img = fread(fid, [hdr.samples, hdr.lines], ['*' dtype]);
        img = img.';   % ENVI BSQ → lines x samples
    end
end
