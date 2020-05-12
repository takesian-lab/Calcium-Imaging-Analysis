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

%% Go to Suite2p folder and load Fall.mat

setup = block.setup;
cd(setup.suite2p_path)
Fall = load('Fall.mat'); %Must load like this because iscell is a matlab function and might lead to unexpected errors.


%% Get Frame_set using get_frames_from_Fall

[Frame_set,~] = get_frames_from_Fall(Fall.ops, setup.block_name, displayTable);
setup.Frame_set = Frame_set;

%Check that Frame_set matches timestamp from Bruker function
if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    warning('Frame_set could not be checked against timestamp')
elseif length(Frame_set) ~= length(block.timestamp)
    error('Frame_set does not match timestamp')
end

%% Check ops

block.ops = get_abridged_ops(Fall.ops);

if checkOps
    unmatchingOps = [];
    fields = fieldnames(user_ops);
    for f = 1:numel(fields)
        currentField = fields{f};
        if strcmp(currentField, 'checkOps')
            continue
        elseif isfield(Fall.ops, currentField)
            if user_ops.(currentField) ~= block.ops.(currentField)
                temp = [];
                temp{1,1} = currentField;
                temp{1,2} = user_ops.(currentField);
                temp{1,3} = block.ops.(currentField);
                unmatchingOps = [unmatchingOps; temp];
            end
        end
    end
    
    if ~isempty(unmatchingOps)
        warning('Some ops do not match the user file.')
        table(unmatchingOps)
    else
        disp('Ops match the user file.')
    end
end

%% Pull out data from Fall
% Fall is too big to keep in its entirety (a couple GB), so just keep the data we'll need

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

if isfield(Fall, 'redcell')
    redcell = Fall.redcell; %Not all runs will have red cells
    block.redcell = redcell(keep_ind);
else
    block.redcell = nan;
end

block.setup = setup;
end