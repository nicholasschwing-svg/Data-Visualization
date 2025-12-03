function img8 = toUint8(img)
%TOUINT8 Normalize an arbitrary numeric image to uint8.
%
%   img8 = toUint8(img)
%
%   If img is already uint8, it is returned unchanged. Otherwise the data
%   is linearly scaled to [0,255].

if isa(img,'uint8')
    img8 = img;
    return;
end

img = double(img);
if isempty(img)
    img8 = uint8([]);
    return;
end

lo = min(img(:));
hi = max(img(:));

if hi > lo
    img8 = uint8(255*(img - lo)/(hi - lo));
else
    img8 = uint8(img);
end
end
