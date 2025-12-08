function [dt, modality, pointId] = parseFastFilenameTime(fname)
%PARSEFASTFILENAMETIME Extract datetime + modality from FAST HSI filename.
%   [dt, modality, pointId] = parseFastFilenameTime(fname)
%   parses filenames like:
%       2024-11-20_1817_17_point-00_LWIR_cal_hsi.hdr
%       2024-11-20_1817_17_point-00_LWIR_cal_hsi.hsic
%   Returns NaT when the timestamp cannot be parsed.
%
%   modality is returned as upper-case text (e.g., 'LWIR', 'VNIR').
%   pointId captures the optional numeric run identifier after 'point-'.

    % Extract core tokens: date/time, optional point-XX, modality
    tok = regexp(fname, ...
        ['(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_' ...
         '(?<hour>\d{2})(?<min>\d{2})_(?<sec>\d{2})' ...
         '(?:_point-(?<point>\d+))?_' ...
         '(?<modality>[A-Za-z0-9-]+?)_' ...
         'cal_hsi'], ...
        'names');

    if isempty(tok)
        dt       = NaT;
        modality = '';
        pointId  = '';
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

    modality = upper(strrep(tok.modality, '-', ''));
    if isfield(tok, 'point') && ~isempty(tok.point)
        pointId = tok.point;
    else
        pointId = '';
    end
end
