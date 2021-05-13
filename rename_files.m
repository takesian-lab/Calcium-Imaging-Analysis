%% Rename files to be numbered

% Get all files in the current folder
folder =('Z:\Carolyn\2P Imaging data\5HT sensor\Cn0012621F1\2021-04-22\BOT_Cn0012621F1_RF_5HTsensor-001'); 
cd(folder)
files = dir([folder '/*.tif']);
% 
% cd folder
% Loop through each file
A=length(files);
for id = 41130:A
%addpath (folder);
f = files(id).name;
% = sprintf('%06d',id);
f_new = strrep(f,'CnO012621F1','Cn0012621F1'); %strrep(f,'a','b') - within file "f" replace a with b;
movefile(f,f_new);
end
