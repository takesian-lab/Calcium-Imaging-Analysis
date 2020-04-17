function [block] = define_suite2p_singleblock(block)

if ismissing(block.setup.suite2p_path)
    disp('Skipping Suite2p data...');
    return
end

disp('Pulling out Suite2p data...');

%Needed from setup:
%suite2p_path
%imaging_set

%% Go to Suite2p folder and load Fall.mat

setup = block.setup;
cd(setup.suite2p_path)
Fall = load('Fall.mat'); %Must load like this because iscell is a matlab function and might lead to unexpected errors.


%% Get Frame_set using get_frames_from_Fall

Imaging_Block = setup.imaging_set;
showTable = 1;
Frame_set = get_frames_from_Fall(Fall.ops,Imaging_Block,showTable);
setup.Frame_set = Frame_set;

%Check that Frame_set matches timestamp from Bruker function
if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    warning('Frame_set could not be checked against timestamp')
elseif length(Frame_set) ~= length(block.timestamp)
    error('Frame_set does not match timestamp')
end

%% Pull out data from Fall

% Only keep data that was flagged 'is cell' to keep the size down

%block.ops = ops; %ops is too big...

% This does not accommodate for the Python to Matlab off by one error yet.

block.img.meanImg = Fall.ops.meanImg;
block.img.refImg = Fall.ops.refImg;
block.img.max_proj = Fall.ops.max_proj;
block.img.meanImgE = Fall.ops.meanImgE;
block.img.Vcorr = Fall.ops.Vcorr;
%block.img.sdmov = Fall.ops.sdmov; %Files saved with older versions of suite2p dont have this

block.iscell = Fall.iscell;
keep_ind = find(block.iscell(:,1)); %Only keep data from 'is cells'
block.cell_number = keep_ind-1;
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