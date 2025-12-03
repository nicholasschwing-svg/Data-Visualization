function n = cleanNumber(val)
%CLEANNUMBER Parse ENVI-style numeric string value to double.
%
%   ENVI headers sometimes wrap values in { } or use commas.

val = regexprep(val, '[\{\}]', '');
val = regexprep(val, ',', ' ');
val = regexp(val, '^[^;]+', 'match', 'once');
n   = str2double(strtrim(val));
end
