function [block] = define_sound_singleblock(block)

display('Pulling out Bruker data...');

%Needed from setup:
%block_path
%VR_path
%voltage_recording
%stim_protocol

setup = block.setup;
New_sound_times = block.New_sound_times;
start_time = block.start_time;

%% Go to block and VR folders and pull out BOT and voltage recording files

%start with the Voltage recording files - this is what we use to
%synchronize BOT data and Tosca data

%now lets get the relevant data from Voltage recording. This is when Tosca
%sends a 5V burst to the Bruker system. The first time this occurs is t=0 on Tosca

if setup.stim_protocol == 4 || setup.voltage_recording == 0
    cd(setup.VR_path) %Widefield
else
    cd(setup.block_path) %2p
end

filedir = dir;
filenames = {filedir(:).name};
VR_files = filenames(contains(filenames,'VoltageRecording'));
VR_filename = VR_files{contains(VR_files,'.csv')};

display(['Loading ' VR_filename])
M = csvread(VR_filename, 1,0);
start = M(find( M(:,2) > 3, 1 ) )./1000;% find the first time that tosca sends a signal to VoltageRecording (in seconds)
Sound_Time(1,:)=(New_sound_times-start_time)+start;
loco_times=block.loco_times;
locTime2 = start+loco_times;

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
block.timestamp=timestamp;
block.Sound_Time=Sound_Time;
block.locTime=locTime2;

%Record filenames
setup.VR_filename = VR_filename;
setup.BOT_filename = BOT_filename;
block.setup = setup;

end



