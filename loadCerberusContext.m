function ctx = loadCerberusContext(hsicFile)
%LOADCERBERUSCONTEXT Load a CERBERUS .hsic cube and compute a context image.
%
%   ctx = loadCerberusContext(hsicFile)
%
%   Uses the ENVI .hdr next to the .hsic (or .hdr) file and handles
%   BSQ/BIL/BIP interleaves correctly. The output ctx is [lines x samples]
%   and is a simple mean over spectral bands.

    if ~isfile(hsicFile)
        error('loadCerberusContext:FileNotFound', ...
              'File not found: %s', hsicFile);
    end

    [path, name, ext] = fileparts(hsicFile);
    if strcmpi(ext, '.hdr')
        hdrFile = hsicFile;
        dataCandidates = {
            fullfile(path, [name '.hsic']), ...
            fullfile(path, [name '.raw']), ...
            fullfile(path, [name '.img']) ...
            };
        dataFile = '';
        for k = 1:numel(dataCandidates)
            if isfile(dataCandidates{k})
                dataFile = dataCandidates{k};
                break;
            end
        end
        if isempty(dataFile)
            dataFile = hsicFile; % fall back to the provided path
        end
    else
        dataFile = hsicFile;
        if ~strcmpi(ext,'.hsic')
            warning('loadCerberusContext:Extension', ...
                    'Expected .hsic file, got "%s". Proceeding anyway.', ext);
        end
        hdrFile = fullfile(path, [name '.hdr']);
    end

    if ~isfile(hdrFile)
        error('loadCerberusContext:HdrMissing', ...
              'Expected header file not found: %s', hdrFile);
    end

    hdr = readENVIHeader(hdrFile);

    requiredFields = {'samples','lines','bands','dataType','headerOffset','interleave'};
    for k = 1:numel(requiredFields)
        fld = requiredFields{k};
        if ~isfield(hdr,fld) || isempty(hdr.(fld))
            error('loadCerberusContext:BadHeader', ...
                  'Header missing field "%s".', fld);
        end
    end

    [dtype, ~] = enviDataType(hdr.dataType);

    machine = 'ieee-le';
    if isfield(hdr,'byteOrder') && isequal(hdr.byteOrder,1)
        machine = 'ieee-be';
    end

    fid = fopen(dataFile,'r',machine);
    if fid < 0
        error('loadCerberusContext:OpenFail', ...
              'Could not open file: %s', dataFile);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    % Skip header bytes
    fseek(fid, hdr.headerOffset, 'bof');

    % Number of samples in the cube
    nTotal = hdr.samples * hdr.lines * hdr.bands;

    % Read all data as a vector
    raw = fread(fid, nTotal, ['*' dtype]);
    if numel(raw) ~= nTotal
        warning('loadCerberusContext:ShortRead', ...
                'Expected %d samples, read %d. Padding with zeros.', ...
                nTotal, numel(raw));
        raw(end+1:nTotal) = 0;
    end

    % Reshape based on interleave
    inter = lower(strtrim(hdr.interleave));
    switch inter
        case 'bsq'
            % [band][line][sample]; read band by band
            tmp  = reshape(raw, [hdr.samples, hdr.lines, hdr.bands]);  % [samples x lines x bands]
            cube = permute(tmp, [2, 1, 3]);                           % [lines x samples x bands]

        case 'bil'
            % [line][band][sample]; records are "band rows" per line
            % Read as [samples x (bands*lines)], then reshape
            tmp2 = reshape(raw, [hdr.samples, hdr.bands*hdr.lines]);       % [samples x (bands*lines)]
            tmp3 = reshape(tmp2, [hdr.samples, hdr.bands, hdr.lines]);     % [samples x bands x lines]
            cube = permute(tmp3, [3, 1, 2]);                               % [lines x samples x bands]

        case 'bip'
            % [line][sample][band]; pixels are full spectra
            % Read as [bands x (samples*lines)], then reshape
            tmp2 = reshape(raw, [hdr.bands, hdr.samples*hdr.lines]);       % [bands x (samples*lines)]
            tmp3 = reshape(tmp2, [hdr.bands, hdr.samples, hdr.lines]);     % [bands x samples x lines]
            cube = permute(tmp3, [3, 2, 1]);                               % [lines x samples x bands]

        otherwise
            error('loadCerberusContext:Interleave', ...
                  'Unsupported interleave type "%s".', hdr.interleave);
    end

    % Simple context image: mean over bands
    ctx = mean(cube, 3);  % double; use imshow(ctx,[]) for display
end
