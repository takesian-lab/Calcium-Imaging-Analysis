function [block] = define_sound_singleblock(block)
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
constant = block.setup.constant;

%% Go to block and VR folders and pull out BOT and voltage recording files

%start with the Voltage recording files - this is what we use to
%synchronize BOT data and Tosca data

%now lets get the relevant data from Voltage recording. This is when Tosca
%sends a 5V burst to the Bruker system. The first time this occurs is t=0 on Tosca

if ~ismissing(setup.VR_name)
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
    activity = block.loco_data.speed(:); % velocity for whole block
    for i = 1:length(Sound_Time)
        sound = Sound_Time(i);
        window = sound + stopwin;
        loc_trial = Loc_BrukerTime{i}(:);
        [c closest_loc_window] = min(abs(loc_trial(:,1)-window));
         [c closest_sound] = min(abs(loc_trial(:,1)-sound));
     act = abs(block.loc_Trial_activity{i}(closest_sound:closest_loc_window));
     actThrsh = find(act>constant.locoThresh);
     if sum(actThrsh)>1
        active_trials(i) = 1;
     else
         active_trials(i) = 0;
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
BOT_files_ind = contains(BOT_files,'.csv');
if sum(BOT_files_ind) == 1 %Regular BOT
    BOT_filename = BOT_files{BOT_files_ind};
    display(['Loading ' BOT_filename])
    frame_data = csvread(BOT_filename, 1,0);
    timestamp = frame_data(:,1)-frame_data(1,1);% this is where we make that small correction
    % For z-corrected multiplane data we actually don't want to make this
    % correction so we will fix that below in the XML section
    singleBOT = true;
else %T-series or multiplane BOT
    BOT_filename = nan;
    singleBOT = false;
    
    %This method does not work for multiplane data. Use XML timestamps instead (see below)
    %I'll leave this code here in case it comes in handy for something else
    
%     disp(['Loading ' num2str(sum(BOT_files_ind)) ' BOT files'])
%     frame_data = [];
%     for b = 1:length(BOT_files)
%         if BOT_files_ind(b) == 1
%             temp_frame_data = csvread(BOT_files{b}, 1,0);
%             frame_data = [frame_data; temp_frame_data];
%         end
%     end   
%     timestamp = frame_data(:,1)-frame_data(1,1);
end



%% Read XML file
% We read the XML just to record info that might come in handly, like laser
% power, PMT values, etc. For multiplane data, the XML timestamps are required

disp('Reading XML file')
XML_files = filenames(contains(filenames,'.xml'));  
X = xml2struct(XML_files{1});

if ~ismissing(setup.VR_name) %widefield
    XML = struct;
    XML.activeMode          = X.PVScan.PVStateShard.PVStateValue{1,1}.Attributes.value;
    XML.camera_binFactor    = X.PVScan.PVStateShard.PVStateValue{1,3}.Attributes.value;
    XML.camera_exposureMode = X.PVScan.PVStateShard.PVStateValue{1,4}.Attributes.value;
    XML.camera_exposureTime = X.PVScan.PVStateShard.PVStateValue{1,5}.Attributes.value;
    XML.dwellTime           = X.PVScan.PVStateShard.PVStateValue{1,9}.Attributes.value;
    XML.framePeriod         = X.PVScan.PVStateShard.PVStateValue{1,10}.Attributes.value;
    XML.opticalZoom         = X.PVScan.PVStateShard.PVStateValue{1,21}.Attributes.value;
    
    framerate = ceil(1/str2double(XML.framePeriod));
    if framerate ~= 20
        warning(['Check frame rate. Detected rate is: ' num2str(1/str2double(XML.framePeriod))])
    end
    
else %2p
    XML = struct;
    XML.filename         = XML_files{1};
    XML.activeMode       = X.PVScan.PVStateShard.PVStateValue{1,1}.Attributes.value;
    XML.dwellTime        = X.PVScan.PVStateShard.PVStateValue{1,5}.Attributes.value;
    XML.framePeriod      = X.PVScan.PVStateShard.PVStateValue{1,6}.Attributes.value;
    XML.laserPower       = X.PVScan.PVStateShard.PVStateValue{1,8}.IndexedValue{1,1}.Attributes.value;
    XML.laserWavelength  = X.PVScan.PVStateShard.PVStateValue{1,9}.IndexedValue.Attributes.value;
    XML.PMT1_Gain        = X.PVScan.PVStateShard.PVStateValue{1,19}.IndexedValue{1,1}.Attributes.value;
    XML.PMT2_Gain        = X.PVScan.PVStateShard.PVStateValue{1,19}.IndexedValue{1,2}.Attributes.value;
    XML.linesPerFrame    = X.PVScan.PVStateShard.PVStateValue{1,10}.Attributes.value;
    XML.micronsPerPixelX = X.PVScan.PVStateShard.PVStateValue{1,12}.IndexedValue{1,1}.Attributes.value;
    XML.micronsPerPixelY = X.PVScan.PVStateShard.PVStateValue{1,12}.IndexedValue{1,2}.Attributes.value;
    XML.opticalZoom      = X.PVScan.PVStateShard.PVStateValue{1,17}.Attributes.value;
    XML.x                = X.PVScan.PVStateShard.PVStateValue{1,20}.SubindexedValues{1,1}.SubindexedValue.Attributes.value;
    XML.y                = X.PVScan.PVStateShard.PVStateValue{1,20}.SubindexedValues{1,2}.SubindexedValue.Attributes.value;
    XML.z                = X.PVScan.PVStateShard.PVStateValue{1,20}.SubindexedValues{1,3}.SubindexedValue{1,1}.Attributes.value;
    XML.piezo            = X.PVScan.PVStateShard.PVStateValue{1,20}.SubindexedValues{1,3}.SubindexedValue{1,2}.Attributes.value;

    framerate = ceil(1/str2double(XML.framePeriod));
    if framerate ~= 15 && framerate ~= 30
        warning(['Check frame rate. Detected rate is: ' num2str(1/str2double(XML.framePeriod))])
    end
    
    %Record timestamp
    if length(X.PVScan.Sequence) == 1 %Single plane data
        nFrames = size(X.PVScan.Sequence.Frame,2); %Should equal number of frames in timestamp
        XML.nPlanes = 1;
        XML.voltageRecording_absoluteTime = str2double(X.PVScan.Sequence.VoltageRecording.Attributes.absoluteTime);
        XML.voltageRecording_relativeTime = str2double(X.PVScan.Sequence.VoltageRecording.Attributes.relativeTime);
        XML.absoluteTime = nan(nFrames,1);
        XML.relativeTime = nan(nFrames,1);
        for f = 1:nFrames
            XML.absoluteTime(f) = str2double(X.PVScan.Sequence.Frame{1,f}.Attributes.absoluteTime);
            XML.relativeTime(f) = str2double(X.PVScan.Sequence.Frame{1,f}.Attributes.relativeTime);
        end
        
        %XML.relativeTime is virtually equivalent to timestamp:
        isequal(nFrames, length(timestamp)); %FYI
        isequal(round(XML.relativeTime, 3), round(timestamp,3)); %FYI
     
    elseif length(X.PVScan.Sequence) > 1 %Multiplane data
        XML.nPlanes = size(X.PVScan.Sequence{1,1}.Frame,2);
        framerate = framerate/XML.nPlanes;
        if strcmp(X.PVScan.Sequence{1,1}.Attributes.bidirectionalZ,'False')
            XML.bidirectionalZ = false;
        else
            XML.bidirectionalZ = true;
        end
        %When VR trigger mode is 'Start with next scan (PFI0)' the VR is aligned with frame #1
        XML.voltageRecording_absoluteTime = str2double(X.PVScan.Sequence{1,1}.VoltageRecording.Attributes.absoluteTime);
        XML.voltageRecording_relativeTime = str2double(X.PVScan.Sequence{1,1}.VoltageRecording.Attributes.relativeTime);
        XML.voltageRecording_triggerMode = X.PVScan.Sequence{1,1}.VoltageRecording.Attributes.triggerMode;
        XML.absoluteTime = [];
        XML.relativeTime = [];
        XML.cycle = [];
        XML.index = [];
        XML.position = []; %Piezo position (values are sometimes inf so we can't get all position data from this) 
        
        count = 1;
        for s = 1:size(X.PVScan.Sequence,2)
            for f = 1:size(X.PVScan.Sequence{1,s}.Frame,2)
                XML.cycle(count,1) = str2double(X.PVScan.Sequence{1,s}.Attributes.cycle);
                if size(X.PVScan.Sequence{1,s}.Frame,2) == 1
                    XML.absoluteTime(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame.Attributes.absoluteTime);
                    XML.relativeTime(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame.Attributes.relativeTime);
                    XML.index(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame.Attributes.index);
                    XML.position(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame.PVStateShard.PVStateValue{1,2}.SubindexedValues{1,3}.SubindexedValue{1,2}.Attributes.value);
                else
                    XML.absoluteTime(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame{1,f}.Attributes.absoluteTime);
                    XML.relativeTime(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame{1,f}.Attributes.relativeTime);
                    XML.index(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame{1,f}.Attributes.index);
                    try
                        XML.position(count,1) = str2double(X.PVScan.Sequence{1,s}.Frame{1,f}.PVStateShard.PVStateValue{1,2}.SubindexedValues{1,3}.SubindexedValue{1,2}.Attributes.value);
                    catch
                        XML.position(count,1) = nan;
                    end
                end
                count = count + 1;
            end
        end
        
        %BOT framerate does not account for the retrace frame that PV adds between cycles
        %XML.relativeTime and absoluteTime show that each cycle takes an extra frame to process
        %This is the same for both one-directional and bidirectional data
        %BOT_fs = unique(round(diff(timestamp),5)); %FYI
        %XML_fs = round(diff(XML.relativeTime),5); %FYI
        %XML_fs_betweenCycles = unique(XML_fs(XML.nPlanes:XML.nPlanes:end)); %FYI
        
        %Replace timestamp with XML.relativeTime ONLY if the BOT data hasn't already been corrected using zcorrect_multiplane (in which
        %case there will only be a single BOT in the folder)
        
        if singleBOT %Z-corrected multiplane data (now single plane)
            timestamp = frame_data; %Use frame_data that hasn't been corrected since the times were originally generated from the XML
        else %Multiplane data
            timestamp.combined = XML.relativeTime - XML.relativeTime(1); %Align timestamp to 0 and save
            for n = 1:XML.nPlanes
                plane = n - 1;
                planeName = ['plane' num2str(plane)];
                if  XML.bidirectionalZ
                    nFrames = length(XML.index);
                    nReps = ceil(nFrames/(XML.nPlanes*2));
                    index = repmat([1:XML.nPlanes, fliplr(1:XML.nPlanes)],1,nReps);
                    index = index(1:nFrames)';
                    planeInd = find(index == n);
                else
                    planeInd = find(XML.index == n);  
                end
                timestamp.(planeName) = timestamp.combined(planeInd);
            end
        end
    end
end
  
disp(['Framerate: ' num2str(framerate)])

%% Record XML data, block data, and filenames

setup.XML = XML;
setup.framerate = framerate;
setup.VR_filename = VR_filename;
setup.BOT_filename = BOT_filename;
block.setup = setup;
block.timestamp = timestamp;
block.active_trials = active_trials;
block.locomotion_trace = locomotion_trace;

end
