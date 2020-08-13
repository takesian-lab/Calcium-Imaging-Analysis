%plot_behavior_maryse_and_esther

%% Get relevant data from block

%Block data
mousename = block.setup.mousename;
expt_date = block.setup.expt_date;

%Stimulus data
allFrequencies = block.TargetFreq;
table = tabulate(allFrequencies); %use the function tabulate to extract how many times each frequency was played
uniqueFrequencies = table(:,1)';
nRepsPerFrequency = table(:,2)'; %number of repetitions per frequency 
holdingPeriod = block.holdingPeriod;

%Behavioral data
outcomes = block.Outcome; %get outcomes (hit, miss, FP, withholds) from block
misses = outcomes == 0;
hits = outcomes == 1;
withholds = outcomes == 3;
FPs = outcomes == 4;

%Reaction times
raw_reactionTimes = block.rxn_time;
raw_reactionTimes(raw_reactionTimes < 0) = nan; %trials with no responses become NaN
raw_reactionTimes = raw_reactionTimes/1000; %convert to seconds
reaction_times = raw_reactionTimes - holdingPeriod; %subtract holding period

%% Plot 1 - Psychometric curve

%Calculate hit rate per frequency
hitsPerFrequency = nan(1,length(uniqueFrequencies)); %make empty vector to fill with data
for i = 1:length(uniqueFrequencies)
    freqResponse = hits(allFrequencies == uniqueFrequencies(i));
    hitsPerFrequency(i) = sum(freqResponse);
end
hitRatePerFrequency = hitsPerFrequency./nRepsPerFrequency;

%set up axes
x = 1:length(hitRatePerFrequency);
y = hitRatePerFrequency;

%make figure
figure; hold on
plot(x,y,'Linewidth',2) %plot psychometric curve
hline(0.5,'r') %plot horizontal red line at 0.5
ylabel('Hit Rate')
xlabel('Alternating frequency')
set(gca, 'XTick', x)
set(gca, 'XTickLabel', uniqueFrequencies)
title([mousename ' ' expt_date])

%% Plot 2 - 


