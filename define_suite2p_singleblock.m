function [block] = define_suite2p_singleblock(block, user_ops, displayTable)
% DOCUMENTATION IN PROGRESS
% 
% This function accesses the Fall.mat suite2p data and stores it in block
% 
% Argument(s): 
%   block (struct)
%   ops (Optional struct) = user-specified ops to check against Fall.ops
%   displayTable (Optional 0 or 1 to display table of blocks & frame numbers
%   from the function get_frames_from_Fall)
% 
% Returns:
%   block (struct)
% 
% Notes:
%
% Variables needed from block.setup:
% -suite2p_path
% -block_name
%
% Uses the function:
% -get_frames_from_Fall
%
% TODO:
% Search 'TODO'

%% Options

%Option to check fall.ops against user-specified ops.mat file
if nargin > 1
    checkOps = user_ops.checkOps;
elseif nargin < 2
    checkOps = 0;
end

%Default is to display block information in table format in get_frames_from_Fall
if nargin < 3
    displayTable = 1;  % Option (0 or 1) to print table of blocks and frame numbers to command line
end

%% Skip this function if Suite2p data is not available

if ismissing(block.setup.suite2p_path)
    disp('Skipping Suite2p data...');
    return
end

disp('Pulling out Suite2p data...');

%% Go to Suite2p folder and load plane0/Fall.mat
%  For multi-plane recordings, there will be one 'plane0' folder per plane

setup = block.setup;
cd(setup.suite2p_path)
[~,currentFolder,~] = fileparts(pwd);

switch currentFolder
    case 'plane0'  %If the user enters a filepath ending in 'plane0', assume there is only one plane
        nPlanes = 1;
    case 'suite2p' %If the user enters a filepath ending in 'suite2p', check how many planes there are
        currentDirectory = dir;
        names = {};
        for n = 1:size(currentDirectory,1)
            names{n,1} = currentDirectory(n).name;
        end
        planeFolders = names(contains(names, 'plane'));
        nPlanes = length(planeFolders);
        if nPlanes > 1
            disp(['Found ' num2str(nPlanes) ' suite2p planes'])
            block.MultiplaneData = true;
        end

        %Regardless of nPlanes, cd to plane0 and use as the default folder for getting ops info
        cd(strcat(setup.suite2p_path, '/plane0'))  
    otherwise
        disp('Check suite2p path. Must end in suite2p or plane0')
        return
end

Fall = load('Fall.mat'); %Must load like this because iscell is a matlab function and might lead to unexpected errors.

%% Get Frame_set using get_frames_from_Fall

[Frame_set,~,Block_position] = get_frames_from_Fall(Fall.ops, setup.block_name, displayTable);
setup.Frame_set = Frame_set;

%Check that Frame_set matches timestamp from Bruker function
if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    warning('Frame_set could not be checked against timestamp')
elseif length(Frame_set) ~= length(block.timestamp)
    error('Frame_set does not match timestamp')
end

%% Check ops

block.ops = get_abridged_ops(Fall.ops);

%Check that neucoeff matches block.ops.neucoeff
if ~isequal(block.setup.constant.neucoeff, block.ops.neucoeff)
    warning(['User neucoeff (' num2str(block.setup.constant.neucoeff) ') does not match Suite2p neucoeff (' num2str(block.ops.neucoeff) ')'])
end

if checkOps
    unmatchingOps = [];
    userOps = [];
    blockOps = [];
    fields = fieldnames(user_ops);
    for f = 1:numel(fields)
        currentField = fields{f};
        if strcmp(currentField, 'checkOps')
            continue
        elseif isfield(block.ops, currentField)
            if user_ops.(currentField) ~= block.ops.(currentField)
                unmatchingOps = [unmatchingOps; string(currentField)];
                userOps = [userOps; user_ops.(currentField)];
                blockOps = [blockOps; block.ops.(currentField)];
            end
        end
    end
    
    if ~isempty(unmatchingOps)
        warning('Some ops do not match the user file.')
        format long
        disp(table(unmatchingOps, userOps, blockOps))
    else
        disp('Ops match the user file.')
    end
end

%% Pull out data from Fall: SINGLE-PLANE DATA
% Fall is too big to keep in its entirety (a couple GB), so just keep the data we'll need

if nPlanes == 1
    %Images for visualization
    block.img.meanImg = Fall.ops.meanImg;
    block.img.refImg = Fall.ops.refImg;
    block.img.max_proj = Fall.ops.max_proj;
    block.img.meanImgE = Fall.ops.meanImgE;
    block.img.Vcorr = Fall.ops.Vcorr;

    %Cell and neuropil data
    %Only keep data from 'is cells' within the block's frame set
    block.iscell = Fall.iscell;
    keep_ind = find(block.iscell(:,1));
    block.cell_number = keep_ind-1; %Subtract 1 for the python to matlab correction of cell label
    block.stat = Fall.stat(1,keep_ind);
    block.F = Fall.F(keep_ind,Frame_set);
    block.Fneu = Fall.Fneu(keep_ind,Frame_set);
    block.spks = Fall.spks(keep_ind,Frame_set);
    redcell = Fall.redcell;
    block.redcell = redcell(keep_ind);
    if any(block.redcell)
        disp('Found red cells')
    end

    %Update zcorr frame set
    if isfield (Fall.ops, 'zcorr')
        disp('Found zcorr')
        block.zcorr = Fall.ops.zcorr(:,Frame_set); %Dimensions are z-stack position vs. frame
    end
end


%% Pull out data from Fall: MULTI-PLANE DATA

if nPlanes > 1   

for n = 1:nPlanes
    currentPlane = strcat('plane', num2str(n - 1));
    disp(['Processing ' currentPlane])
    cd(strcat(setup.suite2p_path, '/', currentPlane))
    Fall = load('Fall.mat');

    %Images for visualization
    block.img.(currentPlane).meanImg = Fall.ops.meanImg;
    block.img.(currentPlane).refImg = Fall.ops.refImg;
    block.img.(currentPlane).max_proj = Fall.ops.max_proj;
    block.img.(currentPlane).meanImgE = Fall.ops.meanImgE;
    block.img.(currentPlane).Vcorr = Fall.ops.Vcorr;

    %Cell and neuropil data
    %Only keep data from 'is cells'
    %Find frame_set for each plane (this is different than the Frame_set above)
    frames_per_folder = Fall.ops.frames_per_folder;
    if Block_position == 1 %or single-block data
        Plane_set = 1:frames_per_folder(Block_position);
    else
        Plane_set = (sum(frames_per_folder(1:Block_position-1)) + 1):sum(frames_per_folder(1:Block_position));
    end
    
    block.iscell.(currentPlane) = Fall.iscell;
    keep_ind = find(block.iscell.(currentPlane)(:,1));
    block.cell_number.(currentPlane) = keep_ind-1; %Subtract 1 for the python to matlab correction of cell label
    block.stat.(currentPlane) = Fall.stat(1,keep_ind);
    block.F.(currentPlane) = Fall.F(keep_ind,Plane_set);
    block.Fneu.(currentPlane) = Fall.Fneu(keep_ind,Plane_set);
    block.spks.(currentPlane) = Fall.spks(keep_ind,Plane_set);
    redcell = Fall.redcell;
    block.redcell.(currentPlane) = redcell(keep_ind);
    if n == 1 && any(block.redcell.(currentPlane)); disp('Found red cells'); end
    
    %Update zcorr frame set
    if isfield (Fall.ops, 'zcorr')
        if n == 1
            block.zcorr = [];
            disp('Found zcorr')
        end
        block.zcorr.(currentPlane) = Fall.ops.zcorr(:,Plane_set); %Dimensions are z-stack position vs. frame
    end
        
end

end
    
%% Save
block.setup = setup;
end