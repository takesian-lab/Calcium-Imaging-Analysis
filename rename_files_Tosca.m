%% Rename files to be numbered

% Get all files in the current folder
folder =('D:\2P analysis\2P local data\Vivek\Tosca Data\YC062419M3_VK_9_28\Session 1'); 
cd (folder);
files = dir([folder '/*.txt']);
% 
% cd folder
% Loop through each file
for id = 1:length(files)
%addpath (folder);
f = files(id).name;
% = sprintf('%06d',id);
f_new = strrep(f,'_VK_9_28',''); %strrep(f,'a','b') - within file "f" replace b with a;
movefile(f,f_new);
end
