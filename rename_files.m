%% Rename files to be numbered

% Get all files in the current folder
folder =('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\NxDC030220F2\Tosca_NxDC030220F2\Session 1'); 
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
f_new = strrep(f,'NxDC030920F2','NxDC030220F2'); %strrep(f,'a','b') - within file "f" replace a with b;
movefile(f,f_new);
end
