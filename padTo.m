function out = padTo(in, targetHW)
%PADTO Pad image to target [H W] with white (255) background.
%
%   out = padTo(in, [H W])

[h,w] = size(in);
H = targetHW(1);
W = targetHW(2);

out = uint8(255*ones(H,W));
out(1:h,1:w) = uint8(in);
end
