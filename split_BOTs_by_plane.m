%% split_BOTs_by_plane
% Split BOTs into separate folders based on imaging plane for continuous
% data collected with the piezo

%% CD to folder in question

[path,folder_name,~] = fileparts(pwd);
filedir = dir;
filenames = {filedir(:).name};
tiffs_str = string(filenames(contains(filenames,'.ome.tif')))';
tiffs = char(tiffs_str);

%The filename format should be as follows:
%BOT_FILENAME-###_Cycle#####_Ch2_######.ome.tif
% ### is the block number
% ##### is the cycle number (used for t-series)
% ###### is the tif/frame number

%Filepath access
planeA = 13; %number of chars to start of frame number
planeB = 8;  %number of chars to end of frame number
cycleA = 23;
cycleB = 19;

filelist_lengths = strlength(deblank(tiffs_str)); %Get length of each filename
L = unique(filelist_lengths); %How many unique lengths are there
if length(L) > 1; error('Check filenames'); end

planes = double(string(tiffs(:,L-planeA:L-planeB)));
cycles = double(string(tiffs(:,L-cycleA:L-cycleB)));
unique_planes = unique(planes);

for p = 1:length(unique_planes)
    mkdir([path '/' folder_name], ['plane' num2str(p-1)]); %Make new folder
    
    %Move files
    current_cycle_files = tiffs_str(planes == unique_planes(p));
    for c = 1:length(current_cycle_files)
        copyfile(current_cycle_files(c), [path '/'  folder_name '/' 'plane' num2str(p-1)])
    end
end

disp('Done processing')