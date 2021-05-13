function A = read_prepended_2d_array(fp, datatype)
% READ_PREPENDED_ARRAY -- read array saved by LV with prepended size.
% Usage: A = read_prepended_array(fp)
%
if nargin < 2, datatype = 'double'; end
sz = fread(fp, [1 2], 'uint32');
if any(sz) == 0, 
   A = [];
   return;
end

A = fread(fp, prod(sz), datatype);
A = reshape(A, fliplr(sz));
