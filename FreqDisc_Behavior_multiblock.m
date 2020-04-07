%% Find blocks of interest/day

Hit_threshold = 0.5;
remove_dailyPrep = 1;

info_path = 'D:/2P analysis/2P local data/Carolyn';
cd(info_path)
Info = importfile('Info');
save_path = 'D:/2P analysis/2P local data/Carolyn/analyzed/Daily Imaging';
cd(save_path)
allfiles=dir('*Block*');

bl=0
for i = 1:size(allfiles,1)
bl=bl+1
name = ({allfiles(bl).name});
load(name{1});
Block_number = sprintf('%03d',bl);
FreqDiscData.(['block' Block_number]) = block;

end
numBlocks = bl;

%% find blocks of interest
for i = 1:numBlocks
    Block_number = sprintf('%03d',i);
    prepTrials(i) = FreqDiscData.(['block' Block_number]).prepTrial;
    HitRate (i) = FreqDiscData.(['block' Block_number]).HitRate;
    HitRate;
    if HitRate(i) >= Hit_threshold
        includeBlock(i) = 1;
    else includeBlock(i) = 0;
        
    end
end
if remove_dailyPrep == 1;
    r = find(prepTrials)
    includeBlock(r) = 0;
end
boi = find(includeBlock);
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









