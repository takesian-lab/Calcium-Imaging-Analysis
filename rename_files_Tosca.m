%% Rename files to be numbered

% Get all files in the current folder
folder =('Z:\Nick\Microglia Project\COVID_timeline\LD031920F3\Tosca_LD031920F3\Session 4'); 
cd (folder);
files = dir([folder '/*.txt']);
% 
% cd folder
% Loop through each file
for id = 1:length(files)
%addpath (folder);
f = files(id).name;
% = sprintf('%06d',id);
f_new = strrep(f,'LD-031920-F3','LD031920F3'); %strrep(f,'a','b') - within file "f" replace a with b;
movefile(f,f_new);
end
