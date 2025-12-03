function [dtype, bps] = enviDataType(code)
%ENVIDATATYPE Map ENVI data type code to MATLAB class and bytes-per-sample.
%
%   [dtype, bps] = enviDataType(code)
%
%   See ENVI documentation for data type codes.

switch code
    case 1
        dtype = 'uint8';   bps = 1;
    case 2
        dtype = 'int16';   bps = 2;
    case 3
        dtype = 'int32';   bps = 4;
    case 4
        dtype = 'single';  bps = 4;
    case 5
        dtype = 'double';  bps = 8;
    case 12
        dtype = 'uint16';  bps = 2;
    case 13
        dtype = 'uint32';  bps = 4;
    case 14
        dtype = 'int64';   bps = 8;
    case 15
        dtype = 'uint64';  bps = 8;
    otherwise
        error('enviDataType:Unsupported', ...
              'Unsupported ENVI data type: %d', code);
end
end
