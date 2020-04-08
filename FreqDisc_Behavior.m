
%%setup.stim_protocol=7;
Hit_threshold = 0.5;
Ignore_prep_trial = 1;
by_date=1;

%% find prep trials
findPrep = block.trialType;
r = ismember(0, findPrep);
if r ==0
    block.prepTrial = 1;
else block.prepTrial = 0;
end

%% find hitrate; determine if it is above threshold
Hits = find(block.Outcome==1);
FA = find(block.Outcome==4);
HitRate = length(Hits)./length(block.Outcome);
FARate = length(FA)./length(block.Outcome);

if HitRate >= Hit_threshold
    includeBlock = 1;
else includeBlock = 0;
end
%% visualize the data?
%set sound time to zero
sound_time = block.New_sound_times;
Trial_time = block.Tosca_times;
for i = 1:length(Trial_time)
    center2sound{i} = Trial_time{1,i}(:) - sound_time(i);
    centerLoc (i) =  find(center2sound{i}==0); 
    lickWindow(i,:) = (centerLoc(i)-35):(centerLoc(i)+1000);
    centeredLicks(i,:) = block.lick_time{i,1}(1,lickWindow(1,:));
end

%cs+ vs cs-
CSplus = find(block.trialType==1);
CSminus = find(block.trialType==0);
for i = 1:length(CSplus)
    CSplusLicks(i,:) = centeredLicks(CSplus(i),:);
end
for i = 1:length(CSminus)
    CSminusLicks(i,:) = centeredLicks(CSminus(i),:);
end

timescale = center2sound{1,1}(lickWindow(1,:),1);

subplot(2,1,1)
imagesc(CSplusLicks)
title('Licks aligned to Sound, CS+')


subplot(2,1,2)
imagesc(CSminusLicks)
title('Licks aligned to Sound, CS-')

   
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









