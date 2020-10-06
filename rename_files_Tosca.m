%% Rename files to be numbered

% Get all files in the current folder
folder =('\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\CnL070520F1\Tosca_CnL070520F1\Session 19'); 
cd (folder);
files = dir([folder '/*.txt']);
% 
% cd folder
% Loop through each file
for id = 1:length(files)
%addpath (folder);
f = files(id).name;
% = sprintf('%06d',id);
f_new = strrep(f,'CnL070520M3','CnL070520F1'); %strrep(f,'a','b') - within file "f" replace a with b;
movefile(f,f_new);
end
