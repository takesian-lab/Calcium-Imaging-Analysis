function [block] = define_behavior_singleblock(block)

% DOCUMENTATION IN PROGRESS
%
% This function accesses the stim and locomotor data from the Tosca folder
% and stores it in block
%
% Argument(s):
%   block (struct)
%
% Returns:
%   block (struct)
%
% Notes:
%
% Variables needed from block.setup:
% -Tosca_path
% -Tosca_session
% -Tosca_run
% -mousename
% -stim_protocol
%
% Dependencies: sort_nat
%
% TODO: Remove magic numbers
% Search 'TODO'

%% Skip this function if Tosca data is not available

if ismissing(block.setup.Tosca_path)
    disp('Skipping Tosca data...');
    return
end

disp('Pulling out Tosca data...');

%% Go to Tosca folder and pull out files related to setup.Tosca_run

setup = block.setup;
cd(setup.Tosca_path{1})
allfiles=dir('*Run*');

countblock=1;
counttrials=1;
for n=1:length(allfiles)
    if numel(allfiles(n).name)<=32 %Get all run summaries %27 for other %cgs: numel is the number of array elements
        behaveblock{countblock}=allfiles(n).name;
        countblock=countblock+1;
    elseif numel(allfiles(n).name)>=37 %Get all trials
        trials{counttrials}=allfiles(n).name;
        counttrials=counttrials+1;
    end
end


behaveblock = sort_nat(behaveblock); %Replaced previous code that sorted behaveblock
trials=sort_nat(trials);
%% Pull out loco info, and check that loco trial markers match Tosca trial numbers

Tosca_Run_number = num2str(setup.Tosca_run);
Tosca_Session = num2str(setup.Tosca_session);
mouseID = char(setup.mousename);

%Find behaveblock run number that is equal to Tosca_Run_number
nDigits = length(Tosca_Run_number);
desiredString = ['Run' Tosca_Run_number '.txt'];
behaveblock_runNumbers = {};
for i = 1:length(behaveblock)
    s_index = strfind(behaveblock{i},desiredString); 
    if isempty(s_index)
        behaveblock_runNumbers{i} = nan;
    else
        behaveblock_runNumbers{i} = behaveblock{i}(1,s_index+3:s_index+3+nDigits-1);
    end
end
    
b = find(strcmp(behaveblock_runNumbers,Tosca_Run_number));
if isempty(b) || length(b) > 1 %Run not found or more than one run equals run #
    error('Check Tosca run number.')
end


[Data,Params] = tosca_read_run(behaveblock{b}); %Load block meta-data%%%HACK!!!!!
try
    %     loco_data = dlmread([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']); %locomotor data
    loco_data = tosca_read_loco([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']); %locomotor data
    tloco = loco_data.t(loco_data.ch > 0); % times of trial markers in locomotion data
    mark_loco = find(loco_data.ch>0);
    ntr = length(Data);% number of trials
    %Make sure number of trials from Data matches number of Tosca text files
    %Added by Maryse Oct. 7
    inblock=trials(contains(trials,['Run' Tosca_Run_number '-']));
    ntr_txt = length(inblock);
    if ntr_txt < ntr
        warning([num2str(ntr - ntr_txt) ' less Tosca trial(s) than recorded in Data.']);
        ntr = ntr_txt;
    end
    %
    ttr = NaN(ntr, 1);
    for k = 1:ntr
        tr = tosca_read_trial(Params, Data, k);
        ttr(k) = tr.Time_s(1);
    end
    % The two sets of time stamps should be identical within a few
    % milliseconds. Here, we'll check each locomotion marker and see if there
    % is a real trial starting within 50 ms. If so, keep that locomotion
    % marker. (note, Ken suggested 200ms, but Carolyn changed it on 10/6/20
    % to better correct for a loco error).
    
    ikeep = false(size(tloco));
    for k = 1:length(tloco)
        minDiff = min(abs(tloco(k) - ttr));
        if minDiff < 0.055
            ikeep(k) = true;
        end
    end
    
    tloco = tloco(ikeep); %updated start times for loco data, 9/24/20 cgs
    mark_loco= mark_loco(ikeep);% loco trial start idx, with extra/error loco removed
    
    
catch
    warning('no loco data available')
    loco_data =NaN;
    %     loco_activity = NaN;
    %     loco_times=NaN;
end
%% Read data from the run
Var1=[]; Var2=[];
inblock=trials(contains(trials,['Run' Tosca_Run_number '-'])); %% added hyphen to eliminate double digit spurious entries...
        trialcount=0;
        if length(inblock)>length(Data)
            inblock=inblock(1:end-1);%Hypothesis is trial 00 is generated abberantly, so start on trial 1
        end
for t=1:length(inblock) %Hypothesis is trial 00 is generated abberantly, so start on trial 1
    s=tosca_read_trial(Params,Data,t);%the read_trial gives us more info than read_run alone
    if ~isempty(s)
        Tosca_times{t}=s.Time_s; %pulls out the tosca generated timestamps for each trial
        start_time(t)=Tosca_times{1,t}(1,1);
        end_time(t) =Tosca_times{1,t}(1,end);
        zero_times{t}=Tosca_times{1,t}(1,:)-start_time(t);%set the start of each trial to zero and normalize
        licks{t,:}=s.Lickometer;
        rxn_time(t) = s.Rxn_time_ms;
        states{t}=[0 (diff(s.State_Change)>0)];
        if states{t}(:,:)~=1
            StateChange(t,:)=1;
        else
            StateChange(t,:)=find(states{1,t}(:,:)~=0, 1, 'first');
        end
        
        for y=1:length(zero_times) % find the time (in Tosca units) for the new sound
            n=StateChange(y,1);
            New_sound_times(y)=zero_times{1,y}(1,n);
            New_sound_idx(y) = n;
        end

        if setup.stim_protocol == 13
            holdingPeriod(t) = s.Script.output; %Variable for Maryse behavior stim
        end
                
        %Get CS+/CS- results
        try
            if isequal(Data{t}.Result,'Hit')
                b_Outcome{t}=1;
                if setup.stim_protocol == 7
                    targetFreq=s.cue.Signal.Waveform.Frequency_kHz;%pull out the target frequency
                elseif setup.stim_protocol == 13
                    targetFreq(t) = Data{t}.Target_kHz;
                end
                trialType{t}=1;
            elseif isequal(Data{t}.Result,'Miss')
                b_Outcome{t}=0;
                trialType{t}=1;
                if setup.stim_protocol == 13
                    targetFreq(t) = Data{t}.Target_kHz;
                elseif setup.stim_protocol == 9
                    try
                        targetFreq = s.cue.Signal.FMSweep.Rate_oct_s;
                    catch
                        targetFreq = nan;
                    end
                elseif setup.stim_protocol == 7
                    targetFreq=s.cue.Signal.Waveform.Frequency_kHz;
                else
                    targetFreq = s.cue.Signal.Waveform.Frequency_kHz;
                end
            elseif isequal(Data{t}.Result,'Withhold')
                b_Outcome{t}=3;
                trialType{t}=0;
                if setup.stim_protocol == 13
                    targetFreq(t) = Data{t}.Target_kHz;
                end
            elseif isequal(Data{t}.Result,'False Alarm')
                b_Outcome{t}=4;
                trialType{t}=0;
                if setup.stim_protocol == 13
                    targetFreq(t) = Data{t}.Target_kHz;
                end
            else
                b_Outcome{t}=NaN;
                trialType{t}=NaN;
                if setup.stim_protocol == 13
                    targetFreq(t) = Data{t}.Target_kHz;
                end
                
            end
        catch
            warning('Sound Pav tone trials')
            b_Outcome{t}=s.Result;
            trialType{t}=('SoundPav');
        end
        
        % find the loco times that are closest to the trial start times.
        % Use this information to find which locomotor timestamps
        % correspond to each trial
        
        %         when does each loco trial start:
        %         t_starts = find(loco_data(:,2)==1); %trial starts
        
        try
            locTrial_idx{t} = mark_loco(t)+1 : mark_loco(t+1)-1;
        catch
            locTrial_idx{t} = mark_loco(t)+1 : length(loco_data.t);
        end
        
        %each trial's locomotor trials, corrected by zeroing out the start
        %of each trial.
        zero_loc{t} = loco_data.t(locTrial_idx{t}(:)) - start_time(t);
        % divide the loco activity by trials to use in
        % define_sound_singleblock
        activity_trial{t} = loco_data.speed(locTrial_idx{t}(:));
    end
end

% now that all the loco times are corrected per trial, put them
% back together to get a loc trace that is on a correct timescale
loco_trace_times = [];
loco_trace_activity =[];

% added by Maryse 3/26/21 - do the same for licks/trial time
concatenated_trial_times = [];
concatenated_lick_times = [];

for j = 1:length(zero_loc)
    if j == 1
        loc_add = zero_loc{1,j}(:);
        loco_trace_times = [loco_trace_times; loc_add];
        activity_add = activity_trial{1,j}(:);
        loco_trace_activity =[loco_trace_activity;activity_add];

        %added by Maryse
        trial_add = zero_times{1,j}(:);
        lick_add = licks{1,j}(:);
        concatenated_trial_times = [concatenated_trial_times; trial_add];
        concatenated_lick_times = [concatenated_lick_times; lick_add];
    else
        loc_add = zero_loc{1,j}(:) + loc_add(end);
        loco_trace_times = [loco_trace_times; loc_add];
        activity_add = activity_trial{1,j}(:);
        loco_trace_activity =[loco_trace_activity;activity_add];

        %added by Maryse
        trial_add = zero_times{1,j}(:) + trial_add(end);
        lick_add =  licks{j,1}(:);
        concatenated_trial_times = [concatenated_trial_times; trial_add];
        concatenated_lick_times = [concatenated_lick_times; lick_add];
    end
    loco_trace_activity=abs(loco_trace_activity);
end

A=exist('targetFreq');
if A==0
    targetFreq=NaN;
end

%% Extract stimulus-specific variables

V1 = [];
V2 = [];

for m = 1:length(Data)
    if setup.stim_protocol==1 %noiseburst
        V1 = 0;
        V2 = 0;
        break
    elseif setup.stim_protocol == 2 %Receptive field
        V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
        V2(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
    elseif setup.stim_protocol == 4 %widefield,RF
        V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
        V2(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
    elseif setup.stim_protocol == 5 %SAM
        try
            V1(1,m)  = Data{m}.Sound.Signal.SAM.Rate_Hz;
            V2(1,m)  = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
        catch
            %blank trials
            V1(1,m)  = nan;
            V2(1,m)  = nan;
        end
    elseif setup.stim_protocol == 3 %FM sweep
        V1(1,m)  = Data{m}.Sound.Signal.FMSweep.Rate_oct_s;
        V2(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
    elseif setup.stim_protocol == 6 %SAM freq
        try
            V1(1,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
            V2(1,m)  = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
        catch
            %blank trials
            V1(1,m)  = nan;
            V2(1,m)  = nan;
        end
    elseif setup.stim_protocol == 7 %Behavior go/nogo freq. disc.
        try
            V1(1,m)  = Data{m}.cue.Signal.Waveform.Frequency_kHz;
            V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;
        catch
            V1(1,m)  = Params.Output_States(2).StimChans.Stimulus.Waveform.Tone.Frequency_kHz;
            V2(1,m)  = Params.Output_States(2).StimChans.Stimulus.Level.Level;
         end
    elseif setup.stim_protocol == 8 %Behavior ABI
        V1(1,m)  = Data{m}.cue.CurrentSource.Level.dB_re_1_Vrms;
        V2(1,m)  = Data{m}.cue.Signal.Level.dB_SPL;            
    elseif setup.stim_protocol == 9 || setup.stim_protocol == 11 %Random H20 or Air Puffs
        if strcmp(Data{m}.Type, 'CS+')
            type = 1;
        elseif strcmp(Data{m}.Type, 'CS-')
            type = 0;
        else
            type = nan;
        end
        V1(1,m)  = type;
        V2 = 0;
    elseif setup.stim_protocol == 10 %Noiseburst_ITI
        V1(1,m)  = Data{m}.Sound.Signal.Level.dB_SPL; %0dB for no stim, 70dB for stim
        if m == 1 %Stim interval
            %             V2(1,m)  = New_sound_times(m) - start_time;
            V2(1,m)  = nan;
        else
            V2(1,m)  = New_sound_times(m) - New_sound_times(m-1);
        end
    elseif setup.stim_protocol == 12 %Spontaneous
        V1 = 0;
        V2 = 0;
        break
    elseif setup.stim_protocol == 13 %Maryse behavior
        V1(1,m) = Data{1,m}.Standard_kHz;
        V2(1,m) = Data{1,m}.Target_kHz;
        try
            stim_level = Params.Output_States(2).StimChans(1).Stimulus.Level.Level;
        catch
            stim_level = Params.Tosca.Flowchart(3).State.SigMan.Channels(1).Channel.Level.Level;
        end
    else %stim_protocol doeesn't match any of the above
        warning(['stim_protocol ' num2str(setup.stim_protocol) ' does not exist yet'])
        break;
    end
end

%% Check for tosca trials that are errors, and remove them from the data

error_trials = {};
for j=1:length(Data)
    if isequal(Data{j}.Result,'Error')
        error_trials{j}=Data{j}.trial;
    else
        error_trials{j}=NaN;
    end
end
error_trials=cell2mat(error_trials);
k = find(error_trials>0);
k(k > ntr) = []; %Cannot remove trials that don't exist
if ~isempty(k)
    start_time(:,k) = [];
    Tosca_times(:,k)  = [];
    New_sound_times(:,k) = [];
    New_sound_idx(:,k) = [];
    licks(k,:) = [];
    b_Outcome(:,k) = [];
    trialType(:,k) = [];
    if ~isnan(targetFreq)
        targetFreq(:,k) = [];
    end
    zero_loc(:,k) = [];
    activity_trial(:,k) = [];
    rxn_time(:,k) = [];
    locTrial_idx(:,k) = [];
    if setup.stim_protocol>=2
        V1(:,k)=[];
        V2(:,k)=[];
    end
    if setup.stim_protocol == 13
    	holdingPeriod(:,k) = [];
    end
end

Var1=[Var1,V1];
Var2=[Var2,V2];

%% Save everything to block
%Format used to be: data.([mouseID]).(['ImagingBlock' Imaging_Num]).VARIABLE
%And: data.([mouseID]).parameters

block.start_time = start_time;
block.Tosca_times = Tosca_times;
block.errors = k;
block.New_sound_times = New_sound_times;
block.New_sound_idx = New_sound_idx;
block.lick_time = licks;
block.concat_times = concatenated_trial_times;
block.concat_licks = concatenated_lick_times;
block.Outcome =  cell2mat(b_Outcome);
block.trialType = cell2mat(trialType);
block.TargetFreq = targetFreq;
block.parameters.variable1 = Var1; %index of variable 1 (e.g. frequency)
block.parameters.variable2 = Var2; %index of variable 2 (e.g. level)
block.loco_data = loco_data; %raw loco data
block.loco_activity = loco_trace_activity; %trace of velocity
block.loco_times = loco_trace_times; %time-corrected, should match the velocity trace
block.loc_Trial_times = zero_loc; %timestamps for loco by each trial
block.loc_Trial_activity = activity_trial; % time-corrected velocity for each trial
block.rxn_time = rxn_time;
block.setup = setup;
block.locIDX = locTrial_idx;
if setup.stim_protocol == 13
    block.holdingPeriod = holdingPeriod;
    block.stim_level = stim_level;
end
end
