%% Find blocks of interest/day


Hit_threshold = 0.5;
remove_dailyPrep = 1;

%this currently only works for fequency discrimination, but we could make
%this part of the pipeline more generic, and have it combine appropriate
%blocks for all of the different types of behavior
stim_protocol = 7; 

%right now all of our analysis is by date; however, I am putting this in as
%a way to potentially change the way we look at our data
by_date = 1; 

%the save_path is where the data are
block_path = 'D:/2P analysis/2P local data/Carolyn/analyzed/Daily Imaging';
cd(block_path)
allfiles=dir('*Block*');

%loop through all of the blocks, and pull out the ones with the appropriate
%stim protocol
bl=0;
b=0;
for i = 1:size(allfiles,1)
    b=b+1;
    name = ({allfiles(b).name});
    load(name{1});
    
    if block.setup.stim_protocol == stim_protocol;
        bl=bl+1
        Behav_number = sprintf('%03d',bl);
        behavData.(['behavblock' Behav_number]) = block;
        
    end   
end
numBlocks = bl;
%% find blocks of interest
for i = 1:numBlocks
    Behav_number = sprintf('%03d',bl);
    prepTrials(i) = behavData.(['behavblock' Behav_number]).prepTrial;
    HitRate (i) = behavData.(['behavblock' Behav_number]).HitRate;
    dates (i) = behavData.(['behavblock' Behav_number]).setup.expt_date;
    Mousename (i) = behavData.(['behavblock' Behav_number]).setup.mousename;
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
all_mice = unique(Mousename);
%% now sort by mouse then by date
%find blocks where all_mice=Mousename
%find blocks where all_dates=dates - but for each mouse

for i = 1:length(all_mice)
    %should I make this a structure that will make it easy to loop through
    %mice and then loop through blocks, or is it easier to keep it in this
    %big matrix?
    mouse_blocks{i} = find(all_mice(i)==Mousename(:));
    for j = 1:length(mouse_blocks{i}
        %you need to loop through the days
    end
end


%% Analyze session - this will be where FreqDisc will actually start!
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









