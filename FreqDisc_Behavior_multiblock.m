%% Find blocks of interest/day
info_path = 'Z:\Carolyn\Behavior\SSRI_mice\Info_sheets';
compiled_blocks_path = 'Z:\Carolyn\Behavior\SSRI_mice\compiled blocks\Cn0012621M3';
save_path = 'Z:\Carolyn\Behavior\SSRI_mice\Analyzed data';
info_filename = 'Info_Cn0012621M3';

stim_protocol = 7;
Hit_threshold = 0.7;
remove_dailyPrep = 1;

%right now all of our analysis is by date; however, I am putting this in as
%a way to potentially change the way we look at our data
by_date = 1;

% add in Maryse's method for loading data
cd(info_path)
Info = importfile(info_filename);
[data] = fillSetupFromInfoTable_v3(Info, compiled_blocks_path,stim_protocol);

%the save_path is where the data are
% block_path = 'Z:\Carolyn\Behavior\SSRI_mice\compiled blocks\Cn0012621M3';
% cd(block_path)
% allfiles=dir('*Compiled*');
%% fix the looping problem
mouselist = data.setup.mousename;
mice = unique(mouselist);
% date = convertStringsToChars(data.setup.expt_date);
date = string(data.setup.expt_date);
uniquedays = unique(date); %still not working because of dashes - need to be removed?
sessionlist = string(data.setup.Session);
sessions = unique(sessionlist);
Outcome = [];
Freq = [];
for i2 = 1:length(mice) %loop through unique mice
    mouseIDX = find(strcmp(mouselist, mice{i2}));
    
%     for k2 = 1:length(sessions) %loop through unique sessions
%         sesIDX = find(strcmp(sessionlist, sessions{k2}));
        for j2 = 1:length(uniquedays) %loop through unique days
            dayIDX = find(strcmp(date, uniquedays{j2}));
            
            RunsOfInterest = intersect(mouseIDX,dayIDX);
            Outcome = [];
            Freq = [];
            for m2=1:length(RunsOfInterest)
                runnum = num2str(data.setup.Tosca_Runs{(RunsOfInterest(m2))});
                sessionname = strcat('Session', sessionlist{RunsOfInterest(m2)}, '_Run', runnum);
                blockname = strcat('data.', mouselist{RunsOfInterest(m2)}, '.Session', sessionlist{RunsOfInterest(m2)}, '_Run', runnum)
                block = eval(blockname);
                if block.HitRate>Hit_threshold  & block.prepTrial ==0
                    Outcome = [Outcome,block.Outcome];
                    Freq = [Freq, block.parameters.variable1];
                    target = block.TargetFreq;
                end
            end
           


for i=1:length(Freq)
    if Freq(i)>0
        freq1=Freq(1,i);
        log=log10(freq1);
        logtarg=log10(target);
        diffFreq=abs(log-logtarg);
        percFreq=diffFreq/0.3010;
        percFreq=round(percFreq,1)
        Freq(1,i)= percFreq;
    end
    
end

Freq(Freq(:,1) < 0, :) = [];
css=unique(Freq(1,:));

%Cued==session
for f=1:length(css)
    pc_session(f)=nanmean(Outcome(find(Freq(:)==css(f))));
    pc_val = Outcome(find(Freq(:)==css(f)));
    session_Length(f,:)=length(find(Freq(1,:)==css(f)));
    if f>1
        pc_session(f)=pc_session(f)-3;
        pc_val = pc_val -3;
    end
    corrects = find(pc_val==1);
    ypsych_size(f) = length(corrects);
end

%Estimated variance by bootstrap
%SI

for f=1:length(css)
    clear meansamp
    for i=1:10000
        clear samp tempmat
        samp=datasample(1:session_Length(f),20,'Replace',true);
        tempmat=Freq(:,find(Freq(1,:)==css(f)));
        meansamp(i)=nanmean(tempmat(samp));
    end
    var_session(f)=std(meansamp);
end



css(css==0.01)=0;
figure;

pc_session100 = pc_session*100;
errorbar(css,pc_session100,var_session,'-bo','LineWidth',2)
% plot(css,pc_session,'-bo','LineWidth',2)
hold on
%errorbar(ucss,pc_uncued,var_uncued,'-ko','LineWidth',2)
xlabel('Frequency (8 is go)')
ylabel('Go probability')
set(gca,'FontSize',14)
%title(num'cued vs. uncued day 1')
legend({'Cued','Uncued'})

% Fit psychometirc functions
targets = [0.25 0.5 0.75] % 25 50 75 % performanc

weights = ones(1,length(css)) % No weighting

% Fit
[coeffspc_session, ~, curvepc_session, thresholdpc_session] = ...
    FitPsycheCurveLogit_cgs(css, ypsych_size, session_Length, targets);


% Plot psychometic curves
plot(curvepc_session(:,1), curvepc_session(:,2), 'LineStyle', '--')
legend('Performance', 'Fit');

            
            %store data
            datenumber = strcat('day',num2str(datenum(uniquedays{j2})));
            FrequencyDisc.([mice{i2}]).([datenumber]).Outcome = Outcome;
            FrequencyDisc.([mice{i2}]).([datenumber]).Frequencies = Freq;
            FrequencyDisc.([mice{i2}]).([datenumber]).targetFrequency = target;
            FrequencyDisc.([mice{i2}]).([datenumber]).threshold =thresholdpc_session;
            
            
        end
    end



%% pull out thresholds

for i = 1:length(mice)
    daycount = 0;
    for j = 1:length(uniquedays)
        daynumber = strcat('day',num2str(datenum(uniquedays{j})))
        try 
            th = FrequencyDisc.([mice{i}]).([daynumber]).threshold(1,2)
             daycount = daycount+1;
            FrequencyDisc.allThresh(i,daycount) = th;
        catch
            warning('no thresholds for date')
        end
    end
end
figure;
plot (FrequencyDisc.allThresh);



%%
% Save data
% cd(save_path);
% save([folder '\' filename]);









