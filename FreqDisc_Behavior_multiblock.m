%% Find blocks of interest/day
info_path = 'Z:\Carolyn\Behavior\SSRI_mice\Info_sheets';
compiled_blocks_path = 'Z:\Carolyn\Behavior\SSRI_mice\compiled blocks\Cn0012621M3';
save_path = 'Z:\Carolyn\2P Imaging data\5HT sensor\Analyzed Data\Cn0012621F1\Baseline RF';
info_filename = 'Info_Cn0012621M3';

stim_protocol = 7;
Hit_threshold = 0.5;
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
%% get the relevant info for the experiment
%loop through mice
mouselist = data.setup.mousename;
mice = unique(mouselist);
date = convertStringsToChars(data.setup.expt_date);
uniquedays = unique(date); %still not working because of dashes - need to be removed?

for i = 1:length(mice)
%loop through dates

end%mice
%%
% use Maryse's function to load data from info sheet?
bl=0
for i = 1:size(allfiles,1)
    bl=bl+1
    name = ({allfiles(bl).name});
    load(name{1});
    %pull out appropriate variables for analysis. This can be updated later
    %with more info, but it is to keep from carrying over all of the blocks.
    expt_date = block.setup.expt_date;
    expt_date = datenum(expt_date);
    block_date = num2str(expt_date);
    Block_number = block.setup.Tosca_run;
    Outcome = block.Outcome;
    HR = block.HitRate;
    FA = block.FARate;
    TarFreq = block.TargetFrequency;
    Frequencies = block.parameters.variable1;
    FreqDiscData.(['Day' block_date]).(['block' Block_number]).Outcome =Outcome;
    FreqDiscData.(['Day' block_date]).(['block' Block_number]).HitRate = HR;
    FreqDiscData.(['Day' block_date]).(['block' Block_number]).FARate = FA;
    FreqDiscData.(['Day' block_date]).(['block' Block_number]).TarFreq = TarFreq;
    FreqDiscData.(['Day' block_date]).(['block' Block_number]).Frequencies = Frequencies;
    
end
numBlocks = bl;
%% find blocks of interest
for i = 1:numBlocks
    Block_number = sprintf('%03d',i);
    prepTrials(i) = FreqDiscData.(['Day' block_date]).(['block' Block_number]).prepTrial;
    HitRate (i) = FreqDiscData.(['Day' block_date]).(['block' Block_number]).HitRate;
    dates (i) = FreqDiscData.(['Day' block_date]).(['block' Block_number]).setup.expt_date;
    HitRate;
    if HitRate(i) >= Hit_threshold
        includeBlock(i) = 1;
    else includeBlock(i) = 0;
        
    end
end
if remove_dailyPrep == 1;
    
    r = find(prepTrials);
    includeBlock(r) = 0;
    
end
boi = find(includeBlock);
all_dates = unique(dates);
%% now sort by date
% notes April 2021: I need to grab all the important info from each block
% (above section) and save in teh FreqDiscData variable - no need to save
% all of the block data even if they are relatively small. Then I will have
% it loop by date and do everything below. I will pull out the thresholds
% for a given day, and save as a single variable and then plot those once
% this section is done. I just need the date sorting worked out and I think
% it will work nicely.

%% Analyze session - this will be where FreqDisc will actually start!
% session=[];
% count=1;
% for bl=1:length(block)
% for t=1:length(block{bl}.trial)
% %plot(block{bl}.trial{t}.licktimes,ones(1,length(block{bl}.trial{t}.licktimes))*t,'r*')
% %hold on; plot(block{bl}.trial{t}.licktimes2,(ones(1,length(block{bl}.trial{t}.licktimes2))*t)+.5,'bo')
%
%     %[FREQUNCY RESULT]
%      session(count,:)=[block{bl}.trial{t}.freq block{bl}.trial{t}.result bl];
%      count;
%     count=count+1;
%
% end
% end




%% Plot cued vs. uncued


Outcome = block.Outcome;
session = block.parameters.variable1;
target=block.TargetFreq;
for i=1:length(session)
    if session(i)>0
        freq=session(1,i);
        log=log10(freq);
        logtarg=log10(target);
        diffFreq=abs(log-logtarg);
        percFreq=diffFreq/0.3010;
        percFreq=round(percFreq,1)
        session(1,i)= percFreq;
    end
    
end

session(session(:,1) < 0, :) = [];
% session(session(:,1)==16, 1 ) = 1 ;
% session(session(:,1)==11.3100, 1 ) = 0.5 ;
% session(session(:,1)==8, 1 ) = 0.01 ;
% session(session(:,1)==13.4500, 1 ) = 0.75 ;
% session(session(:,1)==9.5100, 1 ) = 0.25 ;
% session(session(:,1)==8.8700, 1 ) = 0.1 ;


css=unique(session(1,:));

%Cued==session
for f=1:length(css)
    f
    pc_session(f)=nanmean(Outcome(find(session(1,:)==css(f))));
    session_Length(f,:)=length(find(session(1,:)==css(f)));
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
        tempmat=session(:,find(session(1,:)==css(f)));
        meansamp(i)=nanmean(tempmat(samp));
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
    FitPsycheCurveLogit_cgs(css, pc_session, weights, targets);


% Plot psychometic curves
plot(curvepc_session(:,1), curvepc_session(:,2), 'LineStyle', '--')
legend('Performance', 'Fit');


%%
% Save data
filename=[mouseID '_Day' Day '.mat']
folder = ['D:\Behavior analysis\Behavior local data\Analyzed_Data\' mouseID ];
cd(folder);
save([folder '\' filename]);









