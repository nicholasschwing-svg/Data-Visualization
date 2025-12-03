function hdr = readENVIHeader(hdrFile)
%READENVIHEADER Read basic ENVI .hdr metadata into a struct.
%
%   hdr = readENVIHeader(hdrFile)
%
%   Fields:
%       samples, lines, bands, dataType, headerOffset, interleave, byteOrder

hdr = struct('samples',[], 'lines',[], 'bands',[], 'dataType',[], ...
             'headerOffset',0, 'interleave','bsq', 'byteOrder',0);

fid = fopen(hdrFile,'r');
if fid < 0
    error('readENVIHeader:OpenFail', 'Cannot open HDR: %s', hdrFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

t = fgetl(fid);
while ischar(t)
    parts = strsplit(t, '=');
    if numel(parts) == 2
        key = strtrim(lower(parts{1}));
        val = strtrim(parts{2});
        switch key
            case 'samples'
                hdr.samples = cleanNumber(val);
            case 'lines'
                hdr.lines = cleanNumber(val);
            case 'bands'
                hdr.bands = cleanNumber(val);
            case 'data type'
                hdr.dataType = cleanNumber(val);
            case 'header offset'
                hdr.headerOffset = cleanNumber(val);
            case 'interleave'
                hdr.interleave = lower(val);
            case 'byte order'
                hdr.byteOrder = cleanNumber(val);
        end
    end
    t = fgetl(fid);
end

if ~strcmpi(hdr.interleave,'bsq')
    warning('readENVIHeader:Interleave', ...
        'Viewer assumes BSQ; detected "%s". Results may be incorrect.', ...
        hdr.interleave);
end
end
