function [block] = define_behavior_singleblock_wisam(block)
% [block]=define_behavior_singleblock_wisam(block)
%
% DOCUMENTATION IN PROGRESS
% 
% This function compiles the behavior for a single block
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
% 
% TODO: Remove magic numbers 
% TODO: Memory Allocation Needed
% TODO: Ask about the string manipulations 
% Search 'TODO'

% Checking 
if ismissing(block.setup.Tosca_path)
    disp('Skipping Tosca data...');
    return
end

disp('Pulling out Tosca data...');

%Needed from setup:
%Tosca_path
%Tosca_session
%Tosca_run
%mousename
%stim_protocol

%% Go to Tosca folder and pull out files related to setup.Tosca_run

% TODO: We are copying the setup here?
% Would it be more clear to access this via 'block.setup'?
setup = block.setup;
% TODO: do we need to cd into the directory?
cd(setup.Tosca_path)
% TODO: What is happening here? 
% Looking for a string pattern in the tosca file names?
allfiles=dir('*Run*');

% Counters
countblock=1;
counttrials=1;

%%%%% Loop over all tosca files
% length(allfiles) = numBlocks*numTrials + 3 (extra files)
% 1) 'mouseID-Session#-Run#.trace.txt'
% 2) 'mouseID-Session#-Run#.loco.txt'
% 3) 'mouseID-Session#-Run#.txt'
% TODO: Clarify what these files are for
for n = 1:length(allfiles)
    % TODO: Magic Number
    % What's happening?
    % Are we missing edge case(s) here? What about 37>=numel(allfiles(n).name)>=32
    if numel(allfiles(n).name)<=32 % Get all run summaries %27 for other 
        % behaveblock{} is a cell array of filename strings (from the tosca folder)
        % length(behaveblock) = numBlocks
        behaveblock{countblock}=allfiles(n).name;
        countblock=countblock+1;
    elseif numel(allfiles(n).name)>=37 % Get all trials
        % trials{} is a cell array of filename strings (from the tosca folder)
        trials{counttrials}=allfiles(n).name;
        counttrials=counttrials+1;
    end
end

% The above method of indexing files does not work if Run # >9
% This is because MATLAB alphabetizes: Run 1, Run 10, Run 11, Run 2, etc.
% Block name should be format: MouseID-Session##-Run#.txt
% Find Run #, which could be an arbitrary number of digits:
%
% TODO: What is the goal here?
TrialString = '.txt';
% Crop off the file extension
% N characters to remove from end of block name
N = length(TrialString); 
% Pre-allocate nan array 
% length(runs) = numBlocks
runs = nan(size(behaveblock));
% Loop over the blocks
for f = 1:length(behaveblock)
    % TODO: What are these variables?
    % an empty string to start concatenating to
    tempRun = '';
    % This is a flag for the while loop below
    isNumber = 1;
    % counter
    A = 0;

    while isNumber == 1 % Check each character at the end of the 

        possibleNum = behaveblock{f}(end-N-A);
        if ~isnan(str2double(possibleNum))
            tempRun = strcat(possibleNum, tempRun); % combine digits in the number
            A = A+1;
        else
            if A == 0
                error('No run number found.')
            end
            isNumber = 0;
        end
    end
    runs(f) = str2double(tempRun);
end

% Reorder behaveblock based on Run #
[~,sortedIndex] = sort(runs);
behaveblock = behaveblock(sortedIndex);

%% Read data from the run
% TODO: What is happening here?

% TODO: What are these variables?
Var1=[]; 
Var2=[];
% TODO: We using a copying of 'block.setup' here?
% b is for block
b = setup.Tosca_run;

% Reads summary data for Tosca run
% TODO: Is this a hack? 
[Data,Params] = tosca_read_run(behaveblock{b}); % Load block meta-data%%%HACK!!!!!
% TODO: Do we want to enforce the fileformats here?
inblock = trials(contains(trials,['Run' num2str(b) '-'])); % added hyphen to eliminate double digit spurious entries...

% TODO: Why do we need to do this?
if length(inblock)>length(Data)
    inblock=inblock(1:end-1);
end

% TODO: What does this comment mean? 
% Hypothesis is trial 00 is generated abberantly, so start on trial 1
for t=1:length(inblock) 
    % The read_trial gives us more info than read_run alone
    s = tosca_read_trial(Params,Data,t);
    if ~isempty(s)
        % TODO: what are the dimensions of this time vector
        % pulls out the tosca generated timestamps for each trial
        Tosca_times{t}=s.Time_s;
        start_time=Tosca_times{1,1}(1,1);
        licks{t,:}=s.Lickometer;
        states{t}=[0 (diff(s.State_Change)>0)];
        if states{t}(:,:)~=1
            StateChange(t,:)=1;
            check_for_error=1000;
        else
            StateChange(t,:)=find(states{1,t}(:,:)~=0, 1, 'first');
        end

        % find the time (in Tosca units) for the new sound
        for y=1:length(Tosca_times) 
            n=StateChange(y,1);
            New_sound_times(y)=Tosca_times{1,y}(1,n);
        end
        
        % Get CS+/CS- results
        if isequal(Data{t}.Result,'Hit')
            b_Outcome{t}=1;
            if setup.stim_protocol == 7
                targetFreq=s.cue.Signal.Waveform.Frequency_kHz;%pull out the target frequency
            end
            trialType{t}=1;
        elseif isequal(Data{t}.Result,'Miss')
            b_Outcome{t}=0;
            trialType{t}=1;
            targetFreq=s.cue.Signal.Waveform.Frequency_kHz;
        elseif isequal(Data{t}.Result,'Withhold')
            b_Outcome{t}=3;
            trialType{t}=0;
        elseif isequal(Data{t}.Result,'False Alarm')
            b_Outcome{t}=4;
            trialType{t}=0;;
        else
            b_Outcome{t}=NaN;
            trialType{t}=NaN;
        end
    end
end

A=exist('targetFreq');
if A==0
    targetFreq=NaN;
end

%% Extract stimulus-specific variables
% TODO: What is happening here?

% What are these variables?
% Is it always 2 or can this be multivariate?
V1 = [];
V2 = [];

for m = 1:length(Data)
    if setup.stim_protocol==1
        V1 = 0;
        V2 = 0;
        break
    elseif setup.stim_protocol == 2
        V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
        V2(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
    elseif setup.stim_protocol == 5
        V1(1,m)  = Data{m}.Sound.Signal.SAM.Rate_Hz;
        V2(1,m)  = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
    elseif setup.stim_protocol == 3 %FM sweep
        V1(1,m)  = Data{m}.Sound.Signal.FMSweep.Rate_oct_s;
        V2(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
    elseif setup.stim_protocol == 6
        V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
        V2(1,m)  = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
    elseif setup.stim_protocol == 7
        V1(1,m)  = Data{m}.cue.Signal.Waveform.Frequency_kHz;
        V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
    elseif setup.stim_protocol == 8
        V1(1,m)  = Data{m}.cue.CurrentSource.Level.dB_re_1_Vrms;
        V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
    end
end

%% Check for tosca trials that are errors, and remove them from the data

error_trials = {};
% Loop over blocks
% length(Data) = numBlocks
for j=1:length(Data)
    if isequal(Data{j}.Result,'Error')
        error_trials{j}=Data{j}.trial;
    else
        error_trials{j}=NaN;
    end
end
error_trials=cell2mat(error_trials);
% TODO: Remove this?
~isnan(error_trials);
k = find(error_trials>0);
if ~isempty(k)
    warning('Error trials found in Tosca data')
end

% TODO: What is happening here?
New_sound_times(:,k)=[];
if setup.stim_protocol>=2
    V1(:,k)=[];
    V2(:,k)=[];
end

Var1=[Var1,V1];
Var2=[Var2,V2];

%% Pull out loco info
% TODO: What are the dimensions of the data here?

Tosca_Run_number = num2str(setup.Tosca_run);
Tosca_Session = num2str(setup.Tosca_session);
mouseID = char(setup.mousename);
loco_data = dlmread([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']); %locomotor data
loco_times = loco_data(:,1)-start_time; %I am only looking at column 1
loco_times = loco_times(:,1)+abs(loco_times(1,1));
loco_activity = (abs(loco_data(:,3)));

%% Save everything to block
%Format used to be: data.([mouseID]).(['ImagingBlock' Imaging_Num]).VARIABLE
%And: data.([mouseID]).parameters

block.New_sound_times = New_sound_times;
block.start_time = start_time;
block.lick_time = licks;
block.Tosca_times = Tosca_times;
block.Outcome =  cell2mat(b_Outcome);
block.trialType = cell2mat(trialType);
block.TargetFreq = targetFreq;
block.parameters.variable1 = Var1; %index of variable1 (e.g. frequency)
block.parameters.variable2 = Var2; %index of variable 2 (e.g. level)
block.loco_data = loco_data;
block.loco_activity = loco_activity;
block.loco_times = loco_times;
block.setup = setup;

end