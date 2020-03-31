%% Rename files to be numbered

% Get all files in the current folder
folder =('C:\2P analysis\2P local data\Carolyn\VxAB121018M5\2019-04-01\BOT_VxAB121018M5_2019-04-01-005'); 
% addpath (folder); 
files = dir([folder '/*.tif']);
% 
% cd folder
% Loop through each file
A=length(files);
for id = 1:A
%addpath (folder);
f = files(id).name;
% = sprintf('%06d',id);
f_new = strrep(f,'.ome',''); %strrep(f,'a','b') - within file "f" replace b with a;
movefile(f,f_new);
end
