function img = fridge_read_frame(modality, frameIdx, hdrsMap, filesMap, varargin)
% FRIDGE_READ_FRAME
%   Read a single FRIDGE frame for a given modality and index.
%
% Handles VIS-COLOR (3 bands per frame) and mono-band modalities.
%
% Optional name/value options:
%   'UseCache'      (default true)
%   'MaxCacheBytes' (default 256 MB)
%   'PerfEnabled'   (default false)

    persistent frameCache cacheOrder cacheBytes cacheMaxBytes;

    opts = parseOptions(varargin{:});

    if isempty(frameCache)
        frameCache = containers.Map('KeyType','char','ValueType','any');
        cacheOrder = {};
        cacheBytes = 0;
        cacheMaxBytes = opts.MaxCacheBytes;
    end

    cacheMaxBytes = max(16*1024*1024, opts.MaxCacheBytes);

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

    if ~isfile(file)
        error('Cannot open %s', file);
    end

    cacheKey = '';
    if opts.UseCache
        cacheKey = buildCacheKey(modality, file, frameIdx);
        if isKey(frameCache, cacheKey)
            entry = frameCache(cacheKey);
            img = entry.img;
            touchCacheKey(cacheKey);
            if opts.PerfEnabled
                fprintf('[perf] fridge_read_frame cache hit (%s #%d)\n', modality, frameIdx);
            end
            return;
        end
    end

    tRead = tic;

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

    if opts.PerfEnabled
        fprintf('[perf] fridge_read_frame disk read (%s #%d): %.1f ms\n', ...
            modality, frameIdx, toc(tRead)*1000);
    end

    if opts.UseCache && ~isempty(cacheKey)
        entry = struct('img', img, 'bytes', estimateBytes(img));
        frameCache(cacheKey) = entry;
        cacheOrder{end+1} = cacheKey; %#ok<AGROW>
        cacheBytes = cacheBytes + entry.bytes;
        evictIfNeeded();
    end

    function optsOut = parseOptions(vararginIn)
        optsOut = struct('UseCache', true, ...
                         'MaxCacheBytes', 256*1024*1024, ...
                         'PerfEnabled', false);
        if isempty(vararginIn)
            return;
        end
        for ii = 1:2:numel(vararginIn)
            if ii+1 > numel(vararginIn)
                break;
            end
            key = vararginIn{ii};
            val = vararginIn{ii+1};
            if ~ischar(key) && ~(isstring(key) && isscalar(key))
                continue;
            end
            key = lower(char(key));
            switch key
                case 'usecache'
                    optsOut.UseCache = logical(val);
                case 'maxcachebytes'
                    optsOut.MaxCacheBytes = double(val);
                case 'perfenabled'
                    optsOut.PerfEnabled = logical(val);
            end
        end
    end

    function key = buildCacheKey(modalityIn, fileIn, frameIn)
        info = dir(fileIn);
        bytes = 0;
        stamp = 0;
        if ~isempty(info)
            bytes = info.bytes;
            stamp = posixtime(datetime(info.datenum, 'ConvertFrom','datenum'));
        end
        key = sprintf('%s|%s|%d|%d|%.0f', modalityIn, fileIn, frameIn, bytes, stamp);
    end

    function b = estimateBytes(arr)
        info = whos('arr');
        b = info.bytes;
    end

    function touchCacheKey(key)
        idx = find(strcmp(cacheOrder, key), 1, 'first');
        if ~isempty(idx)
            cacheOrder(idx) = [];
            cacheOrder{end+1} = key;
        end
    end

    function evictIfNeeded()
        while cacheBytes > cacheMaxBytes && ~isempty(cacheOrder)
            victim = cacheOrder{1};
            cacheOrder(1) = [];
            if isKey(frameCache, victim)
                entryVictim = frameCache(victim);
                cacheBytes = max(0, cacheBytes - entryVictim.bytes);
                remove(frameCache, victim);
            end
        end
    end
end
