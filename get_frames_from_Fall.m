function [Frame_set, Frame_set_range] = get_frames_from_Fall(ops, block_name, displayTable) 
% This function gets frame numbers for each recording block from the Suite2p Fall.mat file
% Returns Frame_set and Frame_set_range for the block specified by block_name
% Optionally displays a complete table of block names and frame numbers for the current Fall.mat
%
% Argument(s): 
%   ops (struct from Fall.mat)
%   block_name (optional string) - if no block_name is included, function will
%   display table  without returning anything
%   displayTable (optional 0 or 1) - option to display table
% 
% Returns:
%   Frame_set (double) - full list of frames in current set
%   Frame_set_range (double) - first and last frame in each set
% 
% Notes:
%
%
% TODO: Allow the code accommodate instances where different blocks in
% Fall.mat have the same block_name (i.e. use full filepath to distinguish them)
% Search 'TODO'

%% Define function options

%If no block_name is specified, simply display table
returnFrameSet = 1;
if nargin < 2
    returnFrameSet = 0;
    displayTable = 1;
end

%Default is to display block information in table format
if nargin < 3
    displayTable = 1;
end

%% Magic numbers

redChannel = 1; %Number of the red channel, i.e. Ch1 or Ch2

%Filepath access
frameA = 13; %number of chars to start of frame number
frameB = 8;  %number of chars to end of frame number
blockA = 32; %number of chars to start of block number
blockB = 30; %number of chars to end of block number
chanA  = 15; %number of chars to channel number

%% Extract frame and block numbers from the filenames in ops

% Get filelist from Fall.mat:
filelist_char = ops.filelist;
filelist_str = string(ops.filelist);

%The filelist format should be as follows:
%D:/FILEPATH/BOT_FILENAME-###/BOT_FILENAME-###_Cycle#####_Ch2_######.ome.tif
% ### is the block number
% ##### is the cycle number (used for t-series)
% ###### is the tif/frame number

frames = nan(size(filelist_str)); %Prep variables
blocks = nan(size(filelist_str)); %Prep variables
paths = strings(size(filelist_str)); %Prep variables
chans = nan(size(filelist_str)); %Prep variables

%Accommodate for different filename lengths by processing separately
filelist_lengths = strlength(deblank(filelist_str)); %Get length of each filename
unique_lengths = unique(filelist_lengths); %How many unique lengths are there

for i = 1:length(unique_lengths)
    L = unique_lengths(i);
    currentFiles = filelist_lengths == L;
    frames(currentFiles,:) = double(string(filelist_char(currentFiles,L-frameA:L-frameB)));
    blocks(currentFiles,:) = double(string(filelist_char(currentFiles,L-blockA:L-blockB)));
    paths(currentFiles,:) = string(filelist_char(currentFiles,1:L-blockB));
    chans(currentFiles,:) = double(string(filelist_char(currentFiles,L-chanA)));
end

%If data has both red and green channels, delete channel 1
if sum(chans == redChannel) > 0
    frames(chans == redChannel,:) = [];
    blocks(chans == redChannel,:) = [];
    paths(chans == redChannel,:) = [];
end

%Get first and last frame of each block by finding unique filepaths:
uniquePaths = unique(paths);
blockNumbers = nan(length(uniquePaths),1); %Prep variable
blockFrames = nan(length(uniquePaths),2); %Prep variable

for i = 1:length(uniquePaths)
    currentPath = uniquePaths(i);
    matchingPaths = strcmp(paths, currentPath);
    %If there was an error with some frames because of incorrect filenaming
    %replace block number with Nan to allow the code to finish
    if sum(~isnan(blocks(matchingPaths,:))) == 0
        blockNumbers(i,1) = nan;
        warning('Some files have been assigned Nan as a block number')
    else
        blockNumbers(i,1) = unique(blocks(matchingPaths,:));
    end
    blockFrames(i,1) = find(matchingPaths, 1, 'first');
    blockFrames(i,2) = find(matchingPaths, 1, 'last');
end

%Remove blocks that are only one frame long (this could be from a
%SingleImage file in the same folder as the other data)
isOneFrame = (blockFrames(:,2) - blockFrames(:,1)) == 0;
blockNumbers(isOneFrame,:) = [];
blockFrames(isOneFrame,:) = [];
uniquePaths(isOneFrame,:) = [];

%Get block names by splitting filepaths based on both / and \
if size(uniquePaths,1) == 1
    blockNames_temp1 = split(uniquePaths, '/');
    blockNames_temp2 = split(blockNames_temp1(end), '\');
    blockNames = blockNames_temp2(end);
else
    blockNames_temp1 = split(uniquePaths, '/');
    blockNames_temp2 = split(blockNames_temp1(:,end), '\');
    blockNames = blockNames_temp2(:,end);
end

%Display table of block numbers and frames to user
if displayTable
    format long
    disp(table(blockNames, blockNumbers, blockFrames))
end

if returnFrameSet
    %Convert blockFrames into Frame_set based on block_name
    matching_blocks = strcmp(blockNames, block_name);
    char_block_name = char(block_name);
    
    if isequal(char_block_name(1:10),'Zcorrected')
        matching_blocks = 1; %Accommodate Z-corrected data
    elseif sum(matching_blocks) == 0
        error('block_name is not contained in dataset')
    elseif sum(matching_blocks) > 1
        error('Multiple blocks with the same block_name');
    end
    currentFrames = blockFrames(matching_blocks,:);
    Frame_set = currentFrames(1,1):currentFrames(1,2);
    Frame_set_range = blockFrames(matching_blocks,:); %Use for troubleshooting
    disp(['Current frame set is:   ' num2str(Frame_set_range)])
end


