function name = rmbv_default_filename(prefix, ext)
%RMBV_DEFAULT_FILENAME Build a timestamped default filename for exports.
%   name = RMBV_DEFAULT_FILENAME(prefix, ext) returns a filename of the form
%   "<prefix>_YYYYMMDD_HHMMSS.<ext>". The extension can be supplied with or
%   without a leading dot. If prefix or ext are omitted, sensible defaults
%   are used.

    if nargin < 1 || isempty(prefix)
        prefix = 'montage';
    end
    if nargin < 2 || isempty(ext)
        ext = 'png';
    end

    if startsWith(ext, '.')
        extOut = ext;
    else
        extOut = ['.' ext];
    end

    stamp = datestr(datetime('now'), 'yyyymmdd_HHMMSS');
    name = sprintf('%s_%s%s', prefix, stamp, extOut);
end
