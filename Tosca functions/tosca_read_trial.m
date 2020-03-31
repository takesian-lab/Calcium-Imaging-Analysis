function S = tosca_read_trial(Params, Data, Trial)
% TOSCA_READ_TRIAL -- read detailed data for a single Tosca trial.
% Usage: S = tosca_read_trial(Params, Data)
% 
% Inputs:
%   Params, Data : structures returned by TOSCA_READ_RUN
%   Trial        : trial number
%
% Output:
%   S : structure containing trial data.
%

% Construct trial data filename
[folder, fn] = fileparts(Params.Info.Filename);
fn = fullfile(folder, sprintf('%s-Trial%02d.di.txt', fn, Trial));
if ~exist(fn, 'file'),
   error('File does not exist in current path: %s', fn);
end

fp = fopen(fn, 'rt');

% Parse header row
s = fgetl(fp);
c = textscan(s, '%s', 'delimiter', '\t');
names = c{1};
for k = 1:length(names),
   names{k} = strrep(names{k}, ' ', '_');
   names{k} = strrep(names{k}, '(', '');
   names{k} = strrep(names{k}, ')', '');
end

c = textscan(fp, '%f', 'delimiter', '\t');
nc = length(names);
nr = length(c{1})/nc;

val = reshape(c{1}, nc, nr);
val = val';

S = Data{Trial};
S.DigitalNames = names;

for k = 1:length(names),
   S.(names{k}) = val(:, k)';
end

fclose(fp);
