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

function [info, warnMsg] = rmbv_prepare_video_writer(targetPath, fps)
%RMBV_PREPARE_VIDEO_WRITER Open a VideoWriter or GIF fallback for montage export.
%   [info, warnMsg] = RMBV_PREPARE_VIDEO_WRITER(targetPath, fps) attempts to
%   create a VideoWriter for the requested path. If MPEG-4 is unavailable, the
%   helper falls back to Motion JPEG AVI, and finally to an animated GIF.
%   The returned info struct contains:
%     mode       - 'video' or 'gif'
%     path       - output path actually used
%     writer     - VideoWriter handle when mode == 'video'
%     frameDelay - per-frame delay (for GIF)
%     profile    - profile name used for VideoWriter
%   warnMsg summarizes any fallback decision for user messaging.

    if nargin < 2 || isempty(fps) || ~isnumeric(fps)
        fps = 10;
    end
    fps = max(1, fps);

    info = struct('mode','', 'path','', 'writer',[], 'frameDelay',1/fps, ...
                  'profile','');
    warnMsg = '';

    [baseDir, baseName, ext] = fileparts(targetPath);
    if isempty(ext)
        ext = '.mp4';
        targetPath = fullfile(baseDir, [baseName ext]);
    end

    switch lower(ext)
        case '.gif'
            info.mode = 'gif';
            info.path = targetPath;
            return;

        case '.avi'
            try
                vw = VideoWriter(targetPath, 'Motion JPEG AVI');
                vw.FrameRate = fps;
                open(vw);
                info.mode    = 'video';
                info.path    = targetPath;
                info.writer  = vw;
                info.profile = 'Motion JPEG AVI';
                return;
            catch ME
                warnMsg = sprintf('AVI unavailable (%s). Falling back to GIF.', ME.message);
                % Continue to GIF fallback
            end

        otherwise % mp4 or unknown → try MPEG-4
            try
                vw = VideoWriter(targetPath, 'MPEG-4');
                vw.FrameRate = fps;
                open(vw);
                info.mode    = 'video';
                info.path    = targetPath;
                info.writer  = vw;
                info.profile = 'MPEG-4';
                return;
            catch ME
                warnMsg = sprintf('MP4 unavailable (%s). ', ME.message);
                % Attempt AVI fallback next
                aviPath = fullfile(baseDir, [baseName '.avi']);
                try
                    vw = VideoWriter(aviPath, 'Motion JPEG AVI');
                    vw.FrameRate = fps;
                    open(vw);
                    info.mode    = 'video';
                    info.path    = aviPath;
                    info.writer  = vw;
                    info.profile = 'Motion JPEG AVI';
                    warnMsg = [warnMsg sprintf('Using AVI fallback: %s', aviPath)];
                    return;
                catch ME2
                    warnMsg = [warnMsg sprintf('AVI unavailable (%s). Falling back to GIF.', ME2.message)];
                end
            end
    end

    % Final fallback: GIF
    info.mode = 'gif';
    info.path = fullfile(baseDir, [baseName '.gif']);
end
