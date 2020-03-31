function [Frame_set] = get_frames_from_Fall(ops, Imaging_Block, displayTable)
%Get frame numbers for each recording block from the suite2p Fall.mat file
%Handles only one Imaging_Block at a time
%displayTable is 0 or 1 (default) depending on if you want to view a table of the
%complete block numbers and frames for the current Fall.mat

if nargin < 3
    displayTable = 1;
end

%Get filelist from Fall.mat:
filelist_char = ops.filelist;
filelist_str = string(ops.filelist);

%Extract frame and block numbers from the filenames:

%The filelist format should be as follows:
%D:/FILEPATH/BOT_MOUSEID_DETAILS-###/BOT_MOUSEID_DETAILS-###_Cycle#####_Ch2_######.ome.tif
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
    frames(currentFiles,:) = double(string(filelist_char(currentFiles,L-13:L-8)));
    blocks(currentFiles,:) = double(string(filelist_char(currentFiles,L-32:L-30)));
    paths(currentFiles,:) = string(filelist_char(currentFiles,1:L-34));
    chans(currentFiles,:) = double(string(filelist_char(currentFiles,L-15)));
end

%If data has both red and green channels, delete channel 1
if sum(chans == 1) > 0
    frames(chans == 1,:) = [];
    blocks(chans == 1,:) = [];
    paths(chans == 1,:) = [];
end

%Get first and last frame of each block:
%Find unique filepaths to avoid combining blocks with the same block number
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

%Convert blockFrames into Frame_set based on Imaging_Block

if sum(blockNumbers == Imaging_Block) == 0
    error('Imaging_Block is not contained in dataset')
elseif sum(blockNumbers == Imaging_Block) > 1
    warning('Multiple blocks with the same number');
end
currentFrames = blockFrames(blockNumbers == Imaging_Block,:);
Frame_set = currentFrames(1,1):currentFrames(1,2);

%Display table of block numbers and frames to user
if displayTable
    format long
    display(table(blockNumbers, blockFrames))
end
