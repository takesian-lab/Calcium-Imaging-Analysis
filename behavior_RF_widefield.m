function [parameters] = behavior_RF_widefield(parameters,setup)

for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    Tosca_Session=setup.Session{(a)}
    date=setup.expt_date{(a)};
    Imaging_Block=setup.BOT_maps(a,:)
    Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder
    folder = sprintf([setup.path_name setup.username '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]);
    cd(folder)
    
    
    
    %  Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder
    folder = sprintf([setup.path_name setup.username '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]); %direct to specific Tosca folder within a
    %     addpath (folder); %navigate to this folder
    cd(folder)
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
    behaveblock=sort_nat(behaveblock);
    trials=sort_nat(trials);

    bl=0;
    for b=setup.Tosca_Runs%Loop through blocks of interest and extract trial data + run data, boi was defined above
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
                Tosca_times_looop{t}=s.Loop_time_s;
                parameters.start_time=Tosca_times{1,1}(1,1);
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
                    parameters.New_sound_times(y)=Tosca_times{1,y}(1,n);
                    
                end
                
                
            end
        end
        
    end
    
    %for RF
    if setup.ReceptiveField ==1
        for i=1:length(setup.Tosca_Runs)
            for m = 1:length(Data)
                parameters.frequency_list(i,m) = Data{m}.Sound.Signal.Waveform.Frequency_kHz;
                parameters.level_list (i, m) = Data{m}.Sound.Signal.Level.dB_SPL;
            end
        end
    end
    
    %for water
    if setup.run_water==1;
        for i=1:length(Tosca_Runs)
            for m = 1:length(Data)
                if isequal(Data{m}.Type,'CS+')
                    water_trial(m)=1;
                elseif isequal(Data{m}.Type,'CS-')
                    water_trial(m)=0;
                end
            end
        end
    end
    
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
    
    if setup.ReceptiveField ==1
        parameters.frequency_list(k)=[];
        parameters.level_list(:,k)=[];
    end
    parameters.New_sound_times(:,k)=[];
end
end
