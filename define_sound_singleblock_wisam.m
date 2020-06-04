function [block] = define_sound_singleblock_wisam(block)
% [block]=define_sound_singleblock_wisam(block)
%
% DOCUMENTATION IN PROGRESS
% 
% This function pulls out the Bruker-derived timestamps from BOTs and Voltage Recordings
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
% 
% TODO: Magic numbers
% Search 'TODO'

% Check if both BOT and VR files are both 
if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    disp('Skipping Bruker data...');
    return
end

disp('Pulling out Bruker data...');

%Needed from setup:
%block_path
%VR_path
%voltage_recording
%stim_protocol

%Dependent on Tosca data:
%block.New_sound_times
%block.start_time
%block.loco_data

% TODO: Do we want to copy the setup 
% It might be more clear to call 'block.setup' instead
setup = block.setup;

%% Go to block and VR folders and pull out BOT and voltage recording files

% Start with the Voltage recording files - this is what we use to
% Synchronize BOT data and Tosca data

% Now lets get the relevant data from Voltage recording. This is when Tosca
% sends a 5V burst to the Bruker system. The first time this occurs is t=0 on Tosca

% TODO: Do we need to cd into the directory?
if setup.stim_protocol == 4 || setup.voltage_recording == 0
    cd(setup.VR_path) %Widefield
else
    cd(setup.block_path) %2p
end

% Grab file names from the directory
filedir = dir;
filenames = {filedir(:).name};
VR_files = filenames(contains(filenames,'VoltageRecording'));
VR_filename = VR_files{contains(VR_files,'.csv')};

% Align Tosca data with Bruker times
% TODO: Why is this needed? What is happening?
if ~ismissing(block.setup.Tosca_path) % Skip if Tosca info is missing
    
    % Load variables previously obtained during define_behavior_singleblock
    % These are the stimulus onset times
    New_sound_times = block.New_sound_times;
    start_time = block.start_time;
    % TODO: What determnines these dimensions?
    loco_data = block.loco_data;
    
    % Align with VR csv
    display(['Loading ' VR_filename])
    M = csvread(VR_filename, 1,0);
    % TODO: Check if this timestamp alignment is working correctly
    % TODO: Magic numbers
    % Find the first time that tosca sends a signal to VoltageRecording (in seconds)
    start = M(find( M(:,2) > 3, 1 ) )./1000;
    Sound_Time(1,:) = (New_sound_times-start_time)+start;
    loco_times = block.loco_times;
    locTime2 = start + loco_times;
    
    %%%%%%% START HERE
    %%%%%%% START HERE
    %%%%%%% START HERE
    %%%%%%% START HERE
    %%%%%%% START HERE
    
    % Moved this part from define_loco:
    % TODO: Documentation on this script
    [loco_data,active_time] = locomotor_activity(loco_data,VR_filename);
    block.locomotion_data = loco_data; % TRANSFORMED LOCO DATA
    block.active_time = active_time;

    %Record aligned times in block 
    block.Sound_Time = Sound_Time;
    block.locTime = locTime2;
end

% Let's find the time stamp for each frame. This requires to pull out
% the BOT data and correct for the short (>5ms delay) between Voltage
% recording and BOT

cd(setup.block_path)
filedir = dir;
filenames = {filedir(:).name};
BOT_files = filenames(contains(filenames,'botData'));
BOT_filename = BOT_files{contains(BOT_files,'.csv')};

display(['Loading ' BOT_filename])
frame_data = csvread(BOT_filename, 1,0);
timestamp = frame_data(:,1)-frame_data(1,1);% this is where we make that small correction
block.timestamp = timestamp;

%Record filenames
setup.VR_filename = VR_filename;
setup.BOT_filename = BOT_filename;
block.setup = setup;

end



