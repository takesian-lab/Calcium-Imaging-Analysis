function A = read_prepended_1d_array(fp, datatype)
% READ_PREPENDED_ARRAY -- read array saved by LV with prepended size.
% Usage: A = read_prepended_array(fp)
%

if nargin < 2, datatype = 'uint32'; end

sz = fread(fp, 1, 'uint32');
A = fread(fp, sz, datatype);
