function [Var1,Var2,New_sound_times,start_time,loco_times,loco_activity] = behavior_RF(stim_protocol,a,username,mouseID,Tosca_Session,Tosca_Runs,path_name,Tosca_Run_number,Var1,Var2)

 
%  Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder 
%   folder = sprintf([path_name username '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]);
%     addpath (folder); %navigate to this folder
%   cd(folder)
allfiles=dir('*Run*'); 
   
countblock=1;
counttrials=1; 
for i=1:length(allfiles)
    if numel(allfiles(i).name)<=32 %Get all run summaries %27 for other %cgs: numel is the number of array elements
        behaveblock{countblock}=allfiles(i).name;
        countblock=countblock+1;
    elseif numel(allfiles(i).name)>=37 %Get all trials  
        trials{counttrials}=allfiles(i).name;
        counttrials=counttrials+1;
    end

end
bl=0;
for b=Tosca_Runs%Loop through blocks of interest and extract trial data + run data, boi was defined above 
    bl=bl+1; 
    clear binfo inblock 
    [Data,Params] = tosca_read_run(behaveblock{b}); %Load block meta-data%%%HACK!!!!! 
    inblock=trials(contains(trials,['Run' num2str(b) '-'])); %% added hyphen to eliminate double digit spurious entries... 
    trialcount=0;
    
    if length(inblock)>length(Data)
        inblock=inblock(1:end-1);
    end
   
          for t=1:length(inblock) %Hypothesis is trial00 is generated abberantly, so start on trial 1
            s=tosca_read_trial(Params,Data,t);%the read_trial gives us more info than read_run alone
            if isempty(s)~=1
                
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
        
end

%for RF
for i=1:length(Tosca_Runs)
    for m = 1:length(Data)
        if stim_protocol==1
            V1 =0
            V2  = 0
            if stim_protocol == 2
                V1 (i,m) = Data{m}.Sound.Signal.Waveform;
                V2 (i, m) = Data{m}.Sound.Signal.Level;
            end
        end
    end
end

 %for noiseburst
        
 
 %ERROR: Check for tosca trials that are errors, and remove them from the data
for j=1:length(Data)
    if isequal(Data{j}.Result,'Error')
                 error_trials{j}=Data{j}.trial;
                  else
                error_trials{j}=NaN;
    end
end
error_trials=cell2mat(error_trials);
~isnan(error_trials);
k = find(error_trials>0)
clear error_trials;


New_sound_times(:,k)=[];
if stim_protocol==2
V1(k)=[];
V2(:,k)=[];
Var1=[Var1;V1]
Var2=[Var2;V2];

end
%locomotor data

 loco_data = dlmread([mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt']);%locomotor data  
 loco_times = loco_data(:,1)-loco_data(1,1);% I am only looking at column 1
 loco_activity = (abs(loco_data(:,3)));

end

 