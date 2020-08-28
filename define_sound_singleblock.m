function [block] = define_sound_singleblock(block,constant)
% DOCUMENTATION IN PROGRESS
% 
% This function obtains the Bruker timestamp from the BOT and VR csv files,
% aligns the timestamp with the Tosca data (if available), and stores all
% of this in block
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes: This function uses the Voltage recording file. OVer time, we have
% added addition inputs into our Bruker system. As of 6/25/20, they should
% be the following
% 1)timestamp
% 2)State trigger
% 3) Rep trigger
% 4) trial trigger (we use this to align the timestamps of each trial
% 5) frame trigger from camera (for widefield control only)
% 6) speaker trigger
%   - older data sets will only have columns 1-4-
%
% Variables needed from block.setup:
% -block_path
% -VR_path
% -voltage_recording
% -stim_protocol
%
% Uses the function:
% -locomotor_activity
%
% %If Tosca data is not available these variables won't be created:
% -block.New_sound_times
% -block.start_time
% -block.loco_data
%
% TODO: Remove magic numbers 
% Search 'TODO'

%% Skip this function if Bruker data is not available

if ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    disp('Skipping Bruker data...');
    return
end

disp('Pulling out Bruker data...');

setup = block.setup;

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

% Align Tosca data with Bruker times

if ~ismissing(block.setup.Tosca_path) %Skip if Tosca info is missing
    
    %Load variables previously obtained during define_behavior_singleblock
    New_sound_times = block.New_sound_times;
%     start_time = block.start_time;
    loco_data = block.loco_data;
    
    %Align with VR csv
    display(['Loading ' VR_filename])
    M = csvread(VR_filename, 1,0); %see notes at top for description of M
    
    % find the start of each trial, and align times to it
    Bruker_trial_triggers = M(find( M(:,4) > 3 ) )./1000;
    % each trigger is active for 19-20ms, find the first instance of it 
    diffTrigger = diff(Bruker_trial_triggers);
    Bruker_trial = find(diffTrigger>1)+1;
    % I have the trial start times for all except the first trial, find
    % that trial and add to the list
    Bruker_trial_time = Bruker_trial_triggers(Bruker_trial);
    trial_one = Bruker_trial_triggers(1,1);
    Bruker_trial_time = [trial_one;Bruker_trial_time];
    
    %if there were errors in the run, remove from Bruker trial time
    errors = block.errors;
    Bruker_trial_time(errors,:)=[];
    
    % find the time (normalized from trial start) of sound time
    locomotion_trace =[];
    for i = 1:length(Bruker_trial_time)
        % put sound times on Bruker timescale
        Sound_Time(1,i) = (Bruker_trial_time(i))+New_sound_times(1,i);
        % put locomotor times on Bruker timescale
        Loc_BrukerTime{i} = block.loc_Trial_times{i}(:)+Bruker_trial_time(i);
        % concat loco trials to remake loco trace
        locomotion_trace =[locomotion_trace;Loc_BrukerTime{i}(:)];
      
    end
    
    % set the sound "window" in which to look for locomotor activity. We
    % will use the baseline to end of constant.locowindow. These two
    % numbers are defined at the top of compile_blocks_from_info.m
    base = constant.baseline_length;
    stopwin = constant.locowindow;
    activity = block.loco_data(:,3); % velocity for whole block
    for i = 1:length(Sound_Time)
        sound = Sound_Time(i);
        window = sound + stopwin;
        loc_trial = Loc_BrukerTime{i}(:);
        [c closest_loc_window] = min(abs(loc_trial(:,1)-window));
     act = abs(block.activity_trial{i}(1:closest_loc_window));
     actThrsh = find(act>constant.locoThresh);
     if sum(actThrsh)>1
     active_trials(i) = 1;
     else active_trials(i) = 0;
     end
    
        
    end
   
    loco_times = block.loco_times;
    locTime2 = loco_times; %TODO:needs to be corrected
    
    %Moved this part from define_loco: TODO
    [loco_data,active_time] = locomotor_activity(loco_data,VR_filename,setup);
    block.locomotion_data = loco_data; %TRANSFORMED LOCO DATA
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
block.active_trials = active_trials;
block.locomotion_trace = locomotion_trace;

%Record filenames
setup.VR_filename = VR_filename;
setup.BOT_filename = BOT_filename;
block.setup = setup;

end



