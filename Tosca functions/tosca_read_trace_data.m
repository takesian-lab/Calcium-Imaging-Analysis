function TR = tosca_read_trace_data(FN)
% TOSCA_READ_TRACE_DATA -- read diagnostic trace data
%
% Usage: TR = tosca_read_trace_data(fn)
% 
% Inputs:
%   FN : filename, e.g. 'XYZ999-Session1-Run12.trace.txt'
%
% Output:
%   TR : structure array containing trace data (by trial).
%

if ~contains(FN, '.trace')
   FN = strrep(FN, '.txt', '.trace.txt');
end

fp = fopen(FN, 'rt');
if fp < 0
   error('Could not find file: %s', FN);
end

c = textscan(fp, '%f %d %s %f %s', 'Delimiter', '\t', 'HeaderLines', 1);
fclose(fp);

t = c{1};
e = c{2};
m = c{3};
d = c{4};
s = c{5};

TR = struct( ...
   'Time', [], ...
   'Event', [], ...
   'Message', {}, ...
   'Data', [], ...
   'Source', {} ...
);

itrial = [find(e == 1); length(e)+1];

for kt = 1:length(itrial)-1,
   idx = itrial(kt) : (itrial(kt+1)-1);
   TR(kt).Time = t(idx);
   TR(kt).Event = e(idx);
   TR(kt).Message = m{idx};
   TR(kt).Data = d(idx);
   TR(kt).Source = s{idx};
end
