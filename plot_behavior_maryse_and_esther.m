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
    hits_and_FPs = hits + FPs;
    freqResponse = hits_and_FPs(allFrequencies == uniqueFrequencies(i));
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

% Plot 2 - Reaction time vs. frequency (Esther)

%Calculate average reaction time per frequency
reaction_timesPerFrequency_STD = nan(1,length(uniqueFrequencies));
reaction_timesPerFrequency = nan(1,length(uniqueFrequencies)); %make empty vector to fill with data
for i = 1:length(uniqueFrequencies)
    freqResponse = reaction_times(allFrequencies == uniqueFrequencies(i));
    reaction_timesPerFrequency(i) = nanmean(freqResponse);
    reaction_timesPerFrequency_STD(i) = nanstd(freqResponse);
end

%set up axes
x = 1:length(reaction_timesPerFrequency);
y = reaction_timesPerFrequency;
upper = y + reaction_timesPerFrequency_STD;
lower = y - reaction_timesPerFrequency_STD;

%make figure
%figure; 
subplot(2,1,2); hold on
bar(x,y)
%errorbar(x,y',lower,upper) %Plot standard deviation error bars
scatterx = nan(1,length(reaction_times));
for i = 1:length(uniqueFrequencies)
    scatterx(allFrequencies == uniqueFrequencies(i)) = i;
end
scatter(scatterx,reaction_times, 'k')

ylabel('Reaction Time')
xlabel('Alternating frequency')
set(gca, 'XTick', x)
set(gca, 'XTickLabel', uniqueFrequencies)
title(plotTitle)
xlim([x(1)-1 x(end)+1])

%% Plot 3 - Locomotion and licks

%Get locomotor activity
loco_activity = block.loco_activity;
loco_times = block.loco_times;
ymax = round(1.05*max(loco_activity)); %Determine how high the graph ylim should be

%Get sound and lick times
computer_time = block.Tosca_times{1,1}(1,1);
trial_times = nan(size(block.Tosca_times));
lick_times = [];
for i = 1:length(block.Tosca_times)
    %Get the time each trial starts and subtract computer time
    trial_times(i) = block.Tosca_times{1,i}(1,1) - computer_time;
    %Get the times where the mouse was licking per trial and subtract computer time
    trial_lick_times = block.Tosca_times{1,i}(1,block.lick_time{i,1}==1) - computer_time;
    %Concatenate with previous lick times
    lick_times = [lick_times, trial_lick_times];
end
alternating_sound_times = trial_times + holdingPeriod;
trial_type = block.trialType;

%Determine water times
%Mouse gets water at reaction time if trial is a hit
reaction_times_per_trial = trial_times + reaction_times + holdingPeriod;
water_times = reaction_times_per_trial(hits);

%Plot entire block
figure; hold on
plot(loco_times,loco_activity) %Loco activity is a line graph
xlim([0 loco_times(end)])
ylim([0 ymax])
xlabel('Time (s)')
ylabel('Running speed (cm/s)')
title(plotTitle)

%Plot a vertical line for every sound presentation
for i = 1:length(alternating_sound_times) 
    xsound = alternating_sound_times(i);
    if trial_type(i) == 1 %Target
        lineColor = 'g';
    elseif trial_type(i) == 0 %Non-target
        lineColor = 'r';
    end
        
    line([xsound xsound], [0, ymax], 'Color', lineColor)
end

%Plot a circle for every lick
scatter(lick_times,zeros(size(lick_times))+ymax,'m','Filled')
%Plot a circle for water
scatter(water_times,zeros(size(water_times))+ymax,'c','Filled')

%% Plot example trial(s)

%Adjust xlim and ylim to show only a selection of trials
nTrialsToPlot = 1;
trialToStartWith = 35;

xStart = trial_times(trialToStartWith);
xEnd = trial_times(trialToStartWith + nTrialsToPlot);
y_ind1 = find(loco_times > xStart, 1, 'first');
y_ind2 = find(loco_times > xEnd, 1, 'first');
ymax = round(1.05*max(loco_activity(y_ind1:y_ind2,1)));

figure; hold on
plot(loco_times,loco_activity) %Loco activity is a line graph
xlim([xStart xEnd])
ylim([0 ymax])
xlabel('Time (s)')
ylabel('Running speed (cm/s)')
title(plotTitle)

%Plot a vertical line for every sound presentation
for i = 1:length(alternating_sound_times) 
    xsound = alternating_sound_times(i);
    if trial_type(i) == 1 %Target
        lineColor = 'g';
    elseif trial_type(i) == 0 %Non-target
        lineColor = 'r';
    end
        
    line([xsound xsound], [0, ymax], 'Color', lineColor)
end

%Plot a circle for every lick
scatter(lick_times,zeros(size(lick_times))+ymax,'m','Filled')
%Plot a circle for water
scatter(water_times,zeros(size(water_times))+ymax,'c','Filled')


%% Plot lick rasterplot

%For each trial, extract AltSound - 1 second to AltSound + 3 seconds

toscaFrameRate = loco_times(end)/length(loco_times);
baselinePeriodInSeconds = 4;
holdingPeriodInSeconds = 4;
baselinePeriodInFrames = floor(baselinePeriodInSeconds/toscaFrameRate);
holdingPeriodInFrames = floor(holdingPeriodInSeconds/toscaFrameRate);
totalPeriodInFrames = baselinePeriodInFrames + holdingPeriodInFrames;

trialRaster = zeros(length(trial_times),totalPeriodInFrames + 1);
computer_time = block.Tosca_times{1,1}(1,1);
trial_times = nan(size(block.Tosca_times));
lick_times = [];
for i = 1:length(block.Tosca_times)
    %Get the time each trial starts and subtract computer time
    trial_times(i) = block.Tosca_times{1,i}(1,1) - computer_time;
    %Get the times where the mouse was licking per trial and subtract computer time
    trial_lick_times = block.Tosca_times{1,i}(1,block.lick_time{i,1}==1) - computer_time;
    %Concatenate with previous lick times
    lick_times = [lick_times, trial_lick_times];
end
alternating_sound_times = trial_times + holdingPeriod;
trial_type = block.trialType;
allFrequencies;

for i = 1:length(alternating_sound_times)
    T = alternating_sound_times(i);
    T0 = T - baselinePeriodInSeconds;
    T1 = T + holdingPeriodInSeconds;
    estimatedTimeInFrames = T0:(T1-T0)/80:T1;
    L = lick_times(lick_times > T0);
    L = L(L < T1);
    if isempty(L)
        continue;
    else
        for l = 1:length(L)
            [timeVal, lickInd] = min(abs(estimatedTimeInFrames - L(l)));
            if timeVal < toscaFrameRate
                trialRaster(i,lickInd) = 1;
            end
        end
    end
    
end

[sortedFreqs, sortInd] = sort(allFrequencies);
sortedTrialRaster = trialRaster(sortInd,:);
sortedTrialType = trial_type(sortInd);

%figure; hold on

figure
imagesc(flipud(sortedTrialRaster))
set(gca, 'YTick', 1:length(sortedFreqs))
set(gca, 'YTickLabel', fliplr(sortedFreqs));
ylabel('Frequency')
xlabel('Frames aligned to response window')

figure; hold all %scatterplot
for i = 1:size(sortedTrialRaster,1)
    currentPoints = find(sortedTrialRaster(i,:));
    scatter(currentPoints,zeros(1,length(currentPoints)) + i, 50, 'm', 'filled');
    if ~isempty(currentPoints) && sortedTrialType(i) == 1
        firstPoint = currentPoints(1);
        if firstPoint >= baselinePeriodInFrames
            scatter(firstPoint,zeros(1,length(firstPoint)) + i, 50, 'b', 'filled');
        end
    end
end

set(gca, 'YTick', 1:length(sortedFreqs))
set(gca, 'YTickLabel', sortedFreqs);
ylabel('Frequency')
xlabel('Frames aligned to response window')
