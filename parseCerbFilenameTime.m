function dt = parseCerbFilenameTime(fname, fnameTimePattern)
% parseCerbFilenameTime
% Extracts a datetime from a CERBERUS filename using a regexp pattern.
%
% Example filename:
%   2024-11-19_15-11-38_LWIR_Scan_00198_cal_hsi
%
% Example pattern:
%   '(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_(?<hour>\d{2})-(?<min>\d{2})-(?<sec>\d{2})'
%
% Returns:
%   dt - datetime, or NaT if parsing fails.

    tok = regexp(fname, fnameTimePattern, 'names');

    if isempty(tok)
        dt = NaT;
        return;
    end

    try
        yr  = str2double(tok.year);
        mo  = str2double(tok.month);
        dy  = str2double(tok.day);
        hh  = str2double(tok.hour);
        mm  = str2double(tok.min);
        ss  = str2double(tok.sec);

        dt = datetime(yr, mo, dy, hh, mm, ss);
    catch
        dt = NaT;
    end
end
