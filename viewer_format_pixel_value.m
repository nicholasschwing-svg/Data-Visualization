function s = viewer_format_pixel_value(v)
% VIEWER_FORMAT_PIXEL_VALUE
%   Format scalar / small-vector pixel values for display.

    v = squeeze(v);
    if isscalar(v)
        s = sprintf('%.6g', v);
    else
        v = v(:).';
        if numel(v) <= 6
            s = sprintf('[%s]', ...
                strjoin(arrayfun(@(x)sprintf('%.6g',x), v, 'UniformOutput', false), ', '));
        else
            s = sprintf('[%s, ...] (len=%d)', ...
                strjoin(arrayfun(@(x)sprintf('%.6g',x), v(1:6), 'UniformOutput', false), ', '), ...
                numel(v));
        end
    end
end
