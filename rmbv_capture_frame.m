function [img, method] = rmbv_capture_frame(figHandle, components, forceGetframe)
%RMBV_CAPTURE_FRAME Capture the current viewer into an RGB array.
%   [img, method] = RMBV_CAPTURE_FRAME(figHandle) captures the provided
%   uifigure using exportapp when available (MATLAB R2020a+) and falls back to
%   getframe for older runtimes. If COMPONENTS is provided, the captured image
%   is cropped to the union of the specified UI components so exports can focus
%   on the montage instead of the full application chrome. Set FORCEGETFRAME to
%   true to skip exportapp even when present (useful for performance-oriented
%   captures).

    if nargin < 1 || isempty(figHandle) || ~isvalid(figHandle)
        error('A valid uifigure handle is required for capture.');
    end

    if nargin < 3
        forceGetframe = false;
    end

    method = 'getframe';

    if ~forceGetframe && (exist('exportapp','file') == 2 || exist('exportapp','builtin') == 5)
        tmpPng = [tempname '.png'];
        try
            exportapp(figHandle, tmpPng);
            img = imread(tmpPng);
            delete(tmpPng);
            method = 'exportapp';
        catch
            if isfile(tmpPng)
                delete(tmpPng);
            end
            % Fall back to getframe below
            frameStruct = getframe(figHandle);
            img = frameStruct.cdata;
        end
    else
        frameStruct = getframe(figHandle);
        img = frameStruct.cdata;
    end

    if nargin >= 2 && ~isempty(components)
        img = cropToComponents(img, figHandle, components);
    end
end

function imgOut = cropToComponents(imgIn, figHandle, components)
    imgOut = imgIn;
    if isempty(components)
        return;
    end

    if ~iscell(components)
        components = num2cell(components);
    end

    validComps = components(cellfun(@(c) ~isempty(c) && isvalid(c), components));
    if isempty(validComps)
        return;
    end

    posUnion = [];
    for ii = 1:numel(validComps)
        pos = getpixelposition(validComps{ii}, true);
        if isempty(posUnion)
            posUnion = pos;
        else
            x1 = min(posUnion(1), pos(1));
            y1 = min(posUnion(2), pos(2));
            x2 = max(posUnion(1)+posUnion(3), pos(1)+pos(3));
            y2 = max(posUnion(2)+posUnion(4), pos(2)+pos(4));
            posUnion = [x1, y1, x2 - x1, y2 - y1];
        end
    end

    if isempty(posUnion)
        return;
    end

    pad = 6;  % small padding to avoid clipping borders
    posUnion(1) = posUnion(1) - pad;
    posUnion(2) = posUnion(2) - pad;
    posUnion(3) = posUnion(3) + 2*pad;
    posUnion(4) = posUnion(4) + 2*pad;

    figPos = getpixelposition(figHandle);
    figW = figPos(3);
    figH = figPos(4);

    x1 = max(0, floor(posUnion(1)));
    y1 = max(0, floor(posUnion(2)));
    x2 = min(figW, ceil(posUnion(1) + posUnion(3)));
    y2 = min(figH, ceil(posUnion(2) + posUnion(4)));

    if x2 <= x1 || y2 <= y1
        return;
    end

    % Convert from bottom-left origin (component positions) to image row/col
    % coordinates (top-left origin in cdata).
    rowStart = max(1, size(imgIn,1) - y2 + 1);
    rowEnd   = min(size(imgIn,1), size(imgIn,1) - y1);
    colStart = max(1, x1 + 1);
    colEnd   = min(size(imgIn,2), x2);

    imgOut = imgIn(rowStart:rowEnd, colStart:colEnd, :);
end
