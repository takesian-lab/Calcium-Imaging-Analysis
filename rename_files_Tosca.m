%% Rename files to be numbered

% Get all files in the current folder
folder =('\\apollo\research\ENT\Takesian Lab\Carolyn\Behavior\SSRI_mice\Cn0012621M4\Tosca_Cn0012621M4\Session 15'); 
cd (folder);
files = dir([folder '/*.txt']);
% 
% cd folder
% Loop through each file
for id = 1:length(files)
%addpath (folder);
f = files(id).name;
% = sprintf('%06d',id);
f_new = strrep(f,'Cn0-012621-M4','Cn0012621M4'); %strrep(f,'a','b') - within file "f" replace a with b;
movefile(f,f_new);
end
