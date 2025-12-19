function [img, method] = rmbv_capture_frame(figHandle)
%RMBV_CAPTURE_FRAME Capture the current viewer figure into an RGB array.
%   [img, method] = RMBV_CAPTURE_FRAME(figHandle) captures the provided
%   uifigure using exportapp when available (MATLAB R2020a+) and falls back to
%   getframe for older runtimes. The returned method string notes which path
%   succeeded.

    if nargin < 1 || isempty(figHandle) || ~isvalid(figHandle)
        error('A valid uifigure handle is required for capture.');
    end

    method = 'getframe';

    if exist('exportapp','file') == 2 || exist('exportapp','builtin') == 5
        tmpPng = [tempname '.png'];
        try
            exportapp(figHandle, tmpPng);
            img = imread(tmpPng);
            delete(tmpPng);
            method = 'exportapp';
            return;
        catch
            if isfile(tmpPng)
                delete(tmpPng);
            end
            % Fall back to getframe below
        end
    end

    frameStruct = getframe(figHandle);
    img = frameStruct.cdata;
end
