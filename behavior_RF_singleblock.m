function [data] = behavior_RF_singleblock(setup)

display('Pulling out Tosca data...');

%Needed from setup:
%Tosca_path
%Tosca_session
%Tosca_run
%Imaging_set
%Mousename ??

cd(setup.Tosca_path)

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
       
%The above method of indexing files does not work if Run # >9
%This is because MATLAB alphabetizes: Run 1, Run 10, Run 11, Run 2, etc.
%Block name should be format: MouseID-Session##-Run#.txt
%Find Run #, which could be an arbitrary number of digits:
TrialString = '.txt';
N = length(TrialString); %N characters to remove from end of block name
runs = nan(size(behaveblock));
for f = 1:length(behaveblock)
    tempRun = '';
    isNumber = 1;
    A = 0;
    while isNumber == 1 %Check each character at the end of the 
        possibleNum = behaveblock{f}(end-N-A);
        if ~isnan(str2double(possibleNum))
            tempRun = strcat(possibleNum, tempRun); %combine digits in the number
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
%Reorder behaveblock based on Run #
[~,sortedIndex] = sort(runs);
behaveblock = behaveblock(sortedIndex);

bl=0;
%a
Var1=[]; Var2=[]; isLocoSound=[];sound_v1=[]; sound_v2=[];
Imaging_Block_String = num2str(Imaging_Block(i));
Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
b=setup.Tosca_Runs{a,1}(i);
display(['Processing ' mouseID ': Block ' Imaging_Block_String ' Session ' Tosca_Session ' Run ' num2str(b)]);

bl=bl+1;
clear binfo inblock
[Data,Params] = tosca_read_run(behaveblock{b}); %Load block meta-data%%%HACK!!!!!
inblock=trials(contains(trials,['Run' num2str(b) '-'])); %% added hyphen to eliminate double digit spurious entries...
trialcount=0;

if length(inblock)>length(Data)
    inblock=inblock(1:end-1);
end

for t=1:length(inblock) %Hypothesis is trial 00 is generated abberantly, so start on trial 1
    s=tosca_read_trial(Params,Data,t);%the read_trial gives us more info than read_run alone
    if ~isempty(s)
        Tosca_times{t}=s.Time_s; %pulls out the tosca generated timestamps for each trial
        start_time=Tosca_times{1,1}(1,1);
        %             StateChange(t,:)=find(~s.State_Change,1,'first'); %state change from 1 to zero = when the sound is played
        states{t}=[0 (diff(s.State_Change)>0)];
        if states{t}(:,:)~=1
            StateChange(t,:)=1;
            check_for_error=1000;
        else
            StateChange(t,:)=find(states{1,t}(:,:)~=0, 1, 'first');
        end

        for y=1:length(Tosca_times) % find the time (in Tosca units) for the new sound
            n=StateChange(y,1);
            New_sound_times(y)=Tosca_times{1,y}(1,n);
        end
    end
end



%for RF
%for q=1:length(setup.Tosca_Runs{a,i}(i))
    for m = 1:length(Data)
        if setup.stim_protocol==1
            V1 = 0;
            V2 = 0;
            break
        elseif setup.stim_protocol == 2
            V1(i,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
            V2(i,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
        elseif setup.stim_protocol == 5
            V1(i,m)  = Data{m}.Sound.Signal.SAM.Rate_Hz;
            V2(i,m)  = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
        elseif setup.stim_protocol == 3 %FM sweep
            V1(i,m)  = Data{m}.Sound.Signal.FMSweep.Rate_oct_s;
            V2(i,m)  = Data{m}.Sound.Signal.Level.dB_SPL;
        elseif setup.stim_protocol == 6
            V1(i,m)  = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
            V2(i,m)  = Data{m}.Sound.Signal.SAM.Depth_0_minus1;
        end
    end
%end


%for noiseburst


%ERROR: Check for tosca trials that are errors, and remove them from the data
error_trials = {};
for j=1:length(Data)
    if isequal(Data{j}.Result,'Error')
        error_trials{j}=Data{j}.trial;
    else
        error_trials{j}=NaN;
    end
end
error_trials=cell2mat(error_trials);
~isnan(error_trials);
k = find(error_trials>0);
if ~isempty(k)
    warning('Error trials found in Tosca data')
end
clear error_trials;



New_sound_times(:,k)=[];
if setup.stim_protocol==2
    V1(:,k)=[];
    V2(:,k)=[];
end
data.([mouseID]).(['ImagingBlock' Imaging_Num]).New_sound_times=New_sound_times;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).start_time=start_time;

size(V1);

%make an 'else' for noiseburst?



Var1=[Var1,V1];
Var2=[Var2,V2];
clear V1 V2


%pull out loco info
Tosca_Run_number = num2str(setup.Tosca_Runs{a,1}(i));
loco_data = dlmread([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']);%locomotor data
loco_times = loco_data(:,1)-start_time;% I am only looking at column 1
loco_times = loco_times(:,1)+abs(loco_times(1,1));
loco_activity = (abs(loco_data(:,3)));
data.([mouseID]).(['ImagingBlock' Imaging_Num]).loco_activity=loco_activity;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).loco_times=loco_times;







end
data.([mouseID]).parameters.variable1=Var1;%index of variable1 (frequency)
data.([mouseID]).parameters.variable2=Var2;%index of variable 2 (level)

    %locomotor data
    %     for i=1:length(Imaging_Block(a,:))
    %         Tosca_Run_number = num2str(setup.Tosca_Runs(a,i));
    %         loco_data = dlmread([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']);%locomotor data
    %         loco_times = loco_data(:,1)-start_time;% I am only looking at column 1
    %         loco_times = loco_times(:,1)+abs(loco_times(1,1));
    %         loco_activity = (abs(loco_data(:,3)));
    %         data.([mouseID]).(['ImagingBlock' Imaging_Num]).loco_activity=loco_activity;
    %         data.([mouseID]).(['ImagingBlock' Imaging_Num]).loco_times=loco_times;
    %     end
end
end

