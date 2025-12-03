function nFrames = framesFromFile(modality, hdr, rawPath)
%FRAMESFROMFILE Infer number of frames from RAW file size and header.
%
%   nFrames = framesFromFile(modality, hdr, rawPath)
%
%   For VIS-COLOR, assumes 3 bands per frame if total frames is divisible
%   by 3.

[~, bps] = enviDataType(hdr.dataType);

info = dir(rawPath);
if isempty(info)
    nFrames = NaN;
    return;
end

payload = max(0, info.bytes - hdr.headerOffset);
perFrameBytes = hdr.samples * hdr.lines * bps;
if perFrameBytes <= 0
    nFrames = NaN;
    return;
end

n = floor(double(payload) / double(perFrameBytes));

if strcmp(modality, 'VIS-COLOR') && mod(n,3) == 0
    n = n / 3;
end

nFrames = n;
end
