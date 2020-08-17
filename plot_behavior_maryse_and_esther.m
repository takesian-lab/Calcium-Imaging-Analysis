%plot_behavior_maryse_and_esther

%% Get information from block

%Block data
mousename = char(block.setup.mousename);
expt_date = char(block.setup.expt_date);
dB_level = block.stim_level;
plotTitle = strjoin({mousename, expt_date, num2str(dB_level), 'dB'});

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

%% Plot 1 - Psychometric curve (Maryse)

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
subplot(2,1,1)
plot(x,y,'Linewidth',2) %plot psychometric curve
line([0 x(end)+1], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
ylabel('Hit Rate')
xlabel('Alternating frequency')
set(gca, 'XTick', x)
set(gca, 'XTickLabel', uniqueFrequencies)
title(plotTitle)

%Maryse To Do:
%1. plot psychometric curve with adjustments for trials where the lick was too early
%2. plot S-shaped curve over psychometric curve

%% Plot 2 - Reaction time vs. frequency (Esther)

%Calculate average reaction time per frequency
reaction_timesPerFrequency_STD = nan(1,length(uniqueFrequencies));
reaction_timesPerFrequency = nan(1,length(uniqueFrequencies)); %make empty vector to fill with data
for i = 1:length(uniqueFrequencies)
    freqResponse = reaction_times(allFrequencies == uniqueFrequencies(i));
    reaction_timesPerFrequency(i) = nanmean(freqResponse);
    reaction_timesPerFrequency_STD(i) = nanstd(freqResponse);
end
%reaction_timesPerFrequency = reaction_timesPerFrequency./nRepsPerFrequency;

%set up axes
x = 1:length(reaction_timesPerFrequency);
y = reaction_timesPerFrequency;
upper = y + reaction_timesPerFrequency_STD;
lower = y - reaction_timesPerFrequency_STD;

%make figure
%figure; 
subplot(2,1,2); hold on
bar(x,y)
%errorbar(x,y',lower,upper) %plot psychometric curve

scatterx = nan(1,length(reaction_times));
for i = 1:length(uniqueFrequencies)
    scatterx(allFrequencies == uniqueFrequencies(i)) = i;
end
scatter(scatterx,reaction_times, 'k')

%line([0 x(end)+1], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
ylabel('Reaction Time')
xlabel('Alternating frequency')
set(gca, 'XTick', x)
set(gca, 'XTickLabel', uniqueFrequencies)
title(plotTitle)
xlim([x(1)-1 x(end)+1])