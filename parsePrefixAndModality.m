function [prefix, modality] = parsePrefixAndModality(fname)
%PARSEPREFIXANDMODALITY Extract prefix and modality from RAW filename.
%   [prefix, modality] = parsePrefixAndModality(fname)
%
%   Expected patterns like:
%       PREFIX_LWIR.raw
%       PREFIX_MWIR.raw
%       PREFIX_SWIR.raw
%       PREFIX_MONO.raw
%       PREFIX_VIS-COLOR.raw

tok = regexp(fname, ...
    '^(.*)_(LWIR|MWIR|SWIR|MONO|VIS-COLOR)\.raw$', ...
    'tokens','once');

if isempty(tok)
    prefix   = '';
    modality = '';
else
    prefix   = tok{1};
    modality = tok{2};
end
end
