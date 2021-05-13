function S = tosca_read_ai(Params, Data, Trial)
% TOSCA_READ_AI -- read analog input data for a single Tosca trial.
% Usage: S = tosca_read_ai(Params, Trial)
% 
% Inputs:
%   Params, Data : structure returned by TOSCA_READ_RUN
%   Trial        : trial number
%
% Output:
%   S : structure containing analog input data.
%

if isempty(Data),
   fileNum = Trial;
elseif isfield(Data{Trial}, 'N'),
   fileNum = Data{Trial}.N;
else
   fileNum = Trial;
end

if isnan(fileNum) || Params.AI.active == 0,
   % Immediate error, no data for this trial
   S = [];
   return;
end

% Construct trial data filename
[folder, fn] = fileparts(Params.Info.Filename);
fn = fullfile(folder, sprintf('%s.aidi%02d.bin', fn, fileNum));
if ~exist(fn, 'file'),
   error('File does not exist: %s', fn);
end

data = [];

fp = fopen(fn, 'rb', 'ieee-be');

iorder = read_prepended_1d_array(fp, 'int32');

while ~feof(fp),
   di = read_prepended_2d_array(fp); 
   y = read_prepended_2d_array(fp);
   data = [data; y];
end

fclose(fp);

S.names = Params.AI.AI;
S.t = (0:size(data, 1)-1) / Params.AI.Fs;
S.data = data(:, iorder+1);

itr = find(strcmpi(S.names, 'trial'));
if ~isempty(itr),
   istart = find(diff(S.data(:,itr)) > 3.3, 1, 'first');
   if ~isempty(istart),
      S.t = S.t - S.t(istart);
   end
end
