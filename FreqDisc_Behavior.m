%%Behavioral analysis for frequency discrimination (Coronavirus updates 2020)


%%this intro section is here for the purposes of developing this funciton.
%%Eventually, it will exist within the Behavior_RF function as an option to
%%run or not 

stim_protocol=2;

%Create setup variable for files corresponding to stim_protocol using Info.mat
setup = struct;
setup.username = ''; %'Carolyn'
%setup.path_name = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
setup.path_name = 'D:\2P analysis\2P local data';
setup.stim_protocol = stim_protocol;
setup.run_redcell = 0;
cd(setup.path_name)
load('Info.mat')
setup = fillSetupFromInfoTable(setup, Info);
setup.Info = Info;

%% Experiment mice 
mouseID=setup.mousename{1,1};
sid=setup.Session{};
folder = sprintf(['D:/Behavior analysis/Behavior local data/' mouseID '/' sid]);
boi=[3:14];%blocks of interest (which ones are we actually analyzing)
Day=[('27')];
 
%% Load behavior data

%Read in all data files
%cd(tosca_dir)
cd(folder)
allfiles=dir('*Run*'); %dir will list folder contents with the name Run
countblock=1;
counttrials=1; 
for i=1:length(allfiles)
    if numel(allfiles(i).name)<=32 %Get all run summaries %27 for other %cgs: numel is the number of array elements
        behaveblock{countblock}=allfiles(i).name;
        countblock=countblock+1;
    elseif numel(allfiles(i).name)>=38 %Get all trials  
        trials{counttrials}=allfiles(i).name;
        counttrials=counttrials+1;
    end

end
behaveblock=sort_nat(behaveblock); %cgs-this line puts the trial in
%natural order such that it goes 1,2...10 vs 1, 10, 2... but it already
%seems to do that, and it was causing an error, so I took it out
%trials=sort_nat(trials);


%% Read inclear  blocks of interest.
%to run this section you need the following files " create_valid_varname,
%parse_ini_config, tosca_plot_trial, tosca_read_run, tosca_read_trial,
%tosca_trace. Ask how to direct back to these files. In the mean time, I
%copy them in at this point. 
bl=0;
for b=boi%Loop through blocks of interest and extract trial data + run data, boi was defined above 
    bl=bl+1; 
    bl
    clear binfo inblock 
    [binfo,params] = tosca_read_run(behaveblock{b}); %Load block meta-data%%%HACK!!!!! 
    inblock=trials(contains(trials,['Run' num2str(b) '-'])); %% added hyphen to eliminate double digit spurious entries... 
    trialcount=0;
    
    if length(inblock)>length(binfo)
        inblock=inblock(1:end-1)
    end
%    
        for t=1:length(inblock) %Hypothesis is trial00 is generated abberantly, so start on trial 1
            
            %[block{bl}.trial{t}.Timems, block{bl}.trial{t}.Looptimems, block{bl}.trial{t}.TrialChange, block{bl}.trial{t}.StateChange, block{bl}.trial{t}.RepTrigger, block{bl}.trial{t}.Lickometer, block{bl}.trial{t}.Lick] =
            s=tosca_read_trial(params,binfo,t);
            if isempty(s)~=1
            block{bl}.trial{t}.time=s.Loop_time_s-s.Loop_time_s(1); %align trial times to beginning of trial
            
            block{bl}.trial{t}.licktimes=block{bl}.trial{t}.time(find(s.Lickometer==1));
            block{bl}.trial{t}.StateChange= [0 block{bl}.trial{t}.time(diff(s.State_Change)>0)];
            block{bl}.trial{t}.freq= binfo{t}.cue.Signal.Waveform.Frequency_kHz;
           
            %Get CS+/CS- results
            if isequal(binfo{t}.Result,'Hit')
                block{bl}.trial{t}.result=1;
                target=block{bl}.trial{t}.freq;%pull out the target frequency
            elseif isequal(binfo{t}.Result,'Miss')
                block{bl}.trial{t}.result=0;
            elseif isequal(binfo{t}.Result,'Withhold')
                block{bl}.trial{t}.result=3;
            elseif isequal(binfo{t}.Result,'False Alarm')
                block{bl}.trial{t}.result=4;
            else
                block{bl}.trial{t}.result=NaN;
            end
            
            else 
                block{bl}.trial{t}.result=NaN;
            end
        end
        
    end
        



%% Analyze session 
session=[]; 
count=1;
for bl=1:length(block) 
for t=1:length(block{bl}.trial)
%plot(block{bl}.trial{t}.licktimes,ones(1,length(block{bl}.trial{t}.licktimes))*t,'r*')
%hold on; plot(block{bl}.trial{t}.licktimes2,(ones(1,length(block{bl}.trial{t}.licktimes2))*t)+.5,'bo')
 
    %[FREQUNCY RESULT]
     session(count,:)=[block{bl}.trial{t}.freq block{bl}.trial{t}.result bl];
     count;
    count=count+1;

end 
end

 


%% Plot cued vs. uncued 



for i=1:length(session)
    if session(i,1)>0
        freq=session(i,1);
        log=log10(freq);
        logtarg=log10(target);
        diffFreq=abs(log-logtarg);
        percFreq=diffFreq/0.3010;
        percFreq=round(percFreq,1)
        session(i,1)= percFreq;
    end
    
end

session(session(:,1) < 0, :) = [];
% session(session(:,1)==16, 1 ) = 1 ;
% session(session(:,1)==11.3100, 1 ) = 0.5 ;
% session(session(:,1)==8, 1 ) = 0.01 ;
% session(session(:,1)==13.4500, 1 ) = 0.75 ;
% session(session(:,1)==9.5100, 1 ) = 0.25 ;
% session(session(:,1)==8.8700, 1 ) = 0.1 ;


css=unique(session(:,1));

%Cued==session
for f=1:length(css)
    pc_session(f)=mean(session(find(session(:,1)==css(f)),2));
    session_Length(f)=length(find(session(:,1)==css(f)))
    if f>1
        pc_session(f)=pc_session(f)-3;
    end
end


%Estimated variance by bootstrap
%SI 

for f=1:length(css)
    clear meansamp
    for i=1:10000
        clear samp tempmat
        samp=datasample(1:session_Length(f),20,'Replace',true);
        tempmat=session(find(session(:,1)==css(f)),:);
        meansamp(i)=nanmean(tempmat(samp,2));
    end
    var_session(f)=std(meansamp); 
end

css(css==0.01)=0;



figure;

errorbar(css,pc_session,var_session,'-bo','LineWidth',2)
% plot(css,pc_session,'-bo','LineWidth',2)
hold on 
%errorbar(ucss,pc_uncued,var_uncued,'-ko','LineWidth',2)
xlabel('Frequency (8 is go)')
ylabel('Go probability') 
set(gca,'FontSize',14)
%title(num'cued vs. uncued day 1') 
legend({'Cued','Uncued'}) 







%% Fit Psych Curve

% Fit psychometirc functions
targets = [0.25 0.5 0.75] % 25 50 75 % performance
weights = ones(1,length(css)) % No weighting

% Fit
[coeffspc_session, ~, curvepc_session, thresholdpc_session] = ...
    FitPsycheCurveLogit(css, pc_session, weights, targets);


% Plot psychometic curves
plot(curvepc_session(:,1), curvepc_session(:,2), 'LineStyle', '--')
legend('Performance', 'Fit');



% Save data
filename=[mouseID '_Day' Day '.mat']
folder = ['D:\Behavior analysis\Behavior local data\Analyzed_Data\' mouseID ];
cd(folder);
save([folder '\' filename]); 









