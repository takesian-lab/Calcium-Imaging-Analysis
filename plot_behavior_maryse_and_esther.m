%plot_behavior_maryse_and_esther

%% Get information from block

%Magic numbers
responseWindow = 1;%3; %seconds
dynamicTask = 1; %use for version of task where waiting period dynamically changes

%Block data
mousename = char(block.setup.mousename);
expt_date = char(block.setup.expt_date);
dB_level = block.stim_level;
plotTitle = strjoin({mousename, expt_date, num2str(dB_level), 'dB'});

%Behavioral data
trial_type = block.trialType;
outcomes = block.Outcome; %get outcomes (hit, miss, FP, witholds) from block
earlyLicks = isnan(outcomes); %OR timeouts in some versions of the task
outcomes_removeEL = outcomes(~earlyLicks);
nTrials = length(outcomes_removeEL);
misses      = outcomes_removeEL == 0;
hits        = outcomes_removeEL == 1;
witholds    = outcomes_removeEL == 3;
FPs         = outcomes_removeEL == 4;

%Stimulus data
allFrequencies = log2(block.TargetFreq(~earlyLicks)); %in log scale
repeatingFrequency = mode(allFrequencies);
allFrequencies = round(abs(allFrequencies - repeatingFrequency),3); %in octave difference
table = tabulate(allFrequencies); %use the function tabulate to extract how many times each frequency was played
uniqueFrequencies = table(:,1)';
nRepsPerFrequency = table(:,2)'; %number of repetitions per frequency 
holdingPeriod = block.holdingPeriod(~earlyLicks);

%I removed the holding period in later versions of the task
noHoldingPeriod = false;
if length(unique(holdingPeriod)) == 1
    if unique(holdingPeriod) <= 0.2
        noHoldingPeriod = true;
        holdingPeriod(:) = 0;
    end
end

%In dynamic task, the wait period differs based on mouse's licking behavior
%and extends the duration of the holding period
if dynamicTask
    waitPeriod = block.waitPeriod(~earlyLicks);
    if any(isnan(waitPeriod))
        error('NaNs present in waitPeriod')
    end
    holdingPeriod = holdingPeriod + waitPeriod;
end
    
%Reaction times from the start of the sound
raw_reactionTimes = block.rxn_time(~earlyLicks);
raw_reactionTimes(raw_reactionTimes < 0) = nan; %trials with no responses become NaN
raw_reactionTimes = raw_reactionTimes/1000; %convert to seconds
if noHoldingPeriod || nanmean(raw_reactionTimes) < 5 %In later versions of the Tosca program reaction time doesn't include holding period
    reaction_times = raw_reactionTimes;
else
    reaction_times = raw_reactionTimes - holdingPeriod; %subtract holding period
end

%Make different figures if block is operant or training
if length(unique(block.trialType)) == 1
    blockType = 'Operant';
elseif length(unique(block.trialType)) >= 2 && length(uniqueFrequencies) == 2
    blockType = 'Training';
elseif length(unique(block.trialType)) >= 2 && length(uniqueFrequencies) > 2
    blockType = 'Testing';
else
    error('Could not determine block type')
end

%% Get sound times and make trial rasters:
%  Baseline and response periods are aligned to the response window
%  so if there is a holding period, the baseline will be the repeating tone
%  and if there isn't, the baseline will be silence
    
%Sound times (including early licks) for plotting the full block
%This is when the sound starts during each trial, and is not == to trial start
soundTimes = nan(1,nTrials);
for t = 1:length(block.Tosca_times)
    soundTimes(t) = block.Tosca_times{1,t}(1,block.New_sound_idx(t)) - block.start_time(1);
end

%Convert to indices to prepare for making lick raster and remove early licks
soundTimes_removeEL = soundTimes(~earlyLicks);
responseWindowStart = nan(1,nTrials); 
for i = 1:length(soundTimes_removeEL)
    currentResponseWindowTime = soundTimes_removeEL(i) + holdingPeriod(i);
    [~, responseWindowStart(i)] = min(abs(block.concat_times - currentResponseWindowTime));
end

%How much time to look at
baselinePeriodInSeconds = 5; %seconds before response window
responsePeriodInSeconds = 5; %seconds after response window
validResponsePeriodInSeconds = responseWindow; %seconds where lick sounds as a hit

trialFrameRate = block.concat_times(end)/length(block.concat_times);
baselinePeriodInFrames = round(baselinePeriodInSeconds/trialFrameRate);
responsePeriodInFrames = round(responsePeriodInSeconds/trialFrameRate);
validResponsePeriodInFrames = (validResponsePeriodInSeconds/trialFrameRate);
totalPeriodInFrames = baselinePeriodInFrames + responsePeriodInFrames;

lickRaster = zeros(nTrials,totalPeriodInFrames);
for i = 1:length(responseWindowStart)
    sound_start = responseWindowStart(i);
    a = 1;
    b = totalPeriodInFrames;
    A = sound_start - baselinePeriodInFrames;
    B = sound_start + responsePeriodInFrames - 1;
    %instance where A is before block starts
    if A <= 0
        a = -A + 2;
        A = 1;
    end
    
    %instance where B is after block ends
    if B > length(block.concat_licks)
        B = length(block.concat_licks);
        b = length(A:B);
    end
    lickRaster(i,a:b) = block.concat_licks(A:B,1);
end

%% Operant Plots

if isequal(blockType, 'Operant')
    nTrials = length(hits);
    hitRate = sum(hits)/nTrials;
    missRate = 1 - hitRate;
    hit_Rxn = reaction_times(hits == 1);
        
    figure;
   
    subplot(2,2,1)
    pie([hitRate, missRate], {['Hits (' num2str(hitRate*100) '%)'], ['Misses (' num2str(missRate*100) '%)']})
    title([blockType ' - ' num2str(length(hits)) ' trials'])

    subplot(2,2,3)
    scatter(ones(1,nTrials), reaction_times)
    hline(nanmean(reaction_times),'r')
    ylabel('Time (s)')
    title('Reaction times')
    
    subplot(2,2,2)
    area(sum(lickRaster))
    vline(baselinePeriodInFrames,'r')
    xlabel('Time (ms)')
    ylabel('Licks')
    title('Lick PSTH')
    
    subplot(2,2,4); hold on
    imagesc(lickRaster)
    vline(baselinePeriodInFrames,'r')
    xlabel('Time (ms)')
    ylabel('Trials')
    ylim([1 nTrials])
    title('Lick raster')
    
    suptitle(plotTitle)
end

%% Training

if isequal(blockType, 'Training')
    nTrials = length(hits);
    hitRate = sum(hits)/(sum(hits) + sum(misses));
    FPrate = sum(FPs)/(sum(FPs) + sum(witholds));
    hit_Rxn = reaction_times(hits == 1);
    FP_Rxn = reaction_times(FPs == 1);
        
    figure;
   
    subplot(6,2,[1 3 5])
    bar([1,2], [hitRate, FPrate])
    set(gca, 'XTickLabel', {'Hit Rate', 'FP Rate'})
    ylim([0 1])
    title([blockType ' - ' num2str(length(hits)) ' trials'])

    subplot(6,2,[7 9 11]); hold on
    bar([1,2], [mean(hit_Rxn), mean(FP_Rxn)])
    scatter(ones(1,length(hit_Rxn)), hit_Rxn)
    scatter(ones(1,length(FP_Rxn)) + 1,FP_Rxn)
    set(gca, 'XTick', [1 2])
    set(gca, 'XTickLabel', {'Hit Rate', 'FP Rate'})
    ylabel('Reaction time (s)')
    
    subplot(6,2,[2 4])
    area(sum(lickRaster(hits,:)), 'Facecolor','g')
    vline(baselinePeriodInFrames,'r')
    ylabel('Licks')
    title('Lick PSTH: Hits vs. FPs')
    
    subplot(6,2,[6 8])
    area(sum(lickRaster(FPs,:)), 'Facecolor','r')
    vline(baselinePeriodInFrames,'r')
    ylabel('Licks')
    
    subplot(6,2,[10 12]); hold on
    lickRasterForPlotting = lickRaster;
    lickRasterForPlotting(FPs,:) = -lickRasterForPlotting(FPs,:);
    imagesc(lickRasterForPlotting)
    vline(baselinePeriodInFrames,'r')
    xlabel('Time (ms)')
    ylabel('Trials')
    ylim([1 nTrials])
    title('Lick raster')
    colormap([1 0 0; 1 1 1; 0 1 0] )
    
    suptitle(plotTitle)
end


%% Plot 1 - Psychometric curve (Maryse)

if isequal(blockType, 'Testing')
    
uniqueFrequencies = table(:,1)';
nRepsPerFrequency = table(:,2)';
    
%Calculate hit rate per frequency
hitsPerFrequency = nan(1,length(uniqueFrequencies)); %make empty vector to fill with data
stdPerFrequency = nan(1,length(uniqueFrequencies));
for i = 1:length(uniqueFrequencies)
    hits_and_FPs = hits + FPs;
    freqResponse = hits_and_FPs(allFrequencies == uniqueFrequencies(i));
    hitsPerFrequency(i) = sum(freqResponse);
    stdPerFrequency(i) = std(freqResponse);
end
hitRatePerFrequency = hitsPerFrequency./nRepsPerFrequency;

%REMOVE EASY FREQUENCIES
% hitRatePerFrequency(nRepsPerFrequency == 2) = [];
% hitsPerFrequency(nRepsPerFrequency == 2) = [];
% stdPerFrequency(nRepsPerFrequency == 2) = [];
% uniqueFrequencies(nRepsPerFrequency == 2) = [];
% nRepsPerFrequency(nRepsPerFrequency == 2) = [];

%set up axes
x = 1:length(hitRatePerFrequency);
y = hitRatePerFrequency;

%make figure
figure;
subplot(3,1,1); hold all
errorbar(x, y, stdPerFrequency/2, 'Linewidth', 2)
%plot(x,y,'Linewidth',2) %plot psychometric curve
line([0 x(end)+1], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
line([0 x(end)+1], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
ylabel('Hit Rate')
%xlabel('Frequency difference in octaves')
set(gca, 'XTick', x)
%set(gca, 'XTickLabel', uniqueFrequencies)
ylim([0 1])

subplot(3,1,3); hold all
errorbar(uniqueFrequencies, y, stdPerFrequency/2, 'Linewidth', 2)
%plot(uniqueFrequencies,y,'Linewidth',2) %plot psychometric curve
line([0 uniqueFrequencies(end)], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
ylabel('Hit Rate')
xlabel('Frequency difference in octaves')
ylim([0 1])

% Fit psychometric functions

%For weighted fits
weights = nRepsPerFrequency;
yy = hitsPerFrequency;

%For WH fit
%Mean (u): The mean value of the distribution representing subject bias.
%Standard deviation (v): The variation of the distribution representing the subjects discrimination sensitivity.
%Guess rate (g) and lapse rate (l): Two additional parameters representing the subjects fallibility (ie. potential inability to ever reach 100% performance) at each end of the distribution/stimulus spectrum. 
u = mean(hitRatePerFrequency);
v = std(hitRatePerFrequency);
SPs = [0.1, 0.1,  inf, inf; % Upper limits for g, l, u ,v
       0.01, 0.05, u,  v;  % Start points for g, l, u ,v
       0,    0,    0,  0]; % Lower limits for g, l, u ,v

subplot(3,1,1)
targets = [0.25 0.5 0.75]; % 25 50 75 performance
[~, ~, fit_curve, fit_threshold] = FitPsycheCurveLogit_cgs(x, y, ones(size(uniqueFrequencies)), targets); %FOR GRAPH ONLY
[~, ~, fit_curve2, fit_threshold2] = FitPsycheCurveLogit_cgs(x, yy, weights, targets); %FOR GRAPH ONLY
[~, fit_curve3] = FitPsycheCurveWH(x, y, SPs); %FOR GRAPH ONLY
fit_threshold3 = nan(1,length(targets));
for f = 1:length(targets)
    [~, min_ind] = min(abs(fit_curve3(:,2) - targets(f)));
    fit_threshold3(f) = fit_curve3(min_ind,1);
end

% Plot psychometic curves
%plot(fit_curve(:,1), fit_curve(:,2), 'Linewidth', 2, 'LineStyle', '--', 'Color', 'g')
%plot(fit_curve2(:,1), fit_curve2(:,2)./100, 'Linewidth', 2, 'LineStyle', '--', 'Color', 'm')    
plot(fit_curve3(:,1), fit_curve3(:,2), 'Linewidth', 2, 'LineStyle', '--', 'Color', 'g')    
%scatter(fit_threshold, targets, 'x', 'g')
%scatter(fit_threshold2, targets, 'x', 'm')
scatter(fit_threshold3, targets, 'x', 'k')
%title(['Threshold = ' num2str(fit_threshold3(2))])

subplot(3,1,3)
targets = [0.25 0.5 0.75]; % 25 50 75 performance
[~, ~, fit_curve, fit_threshold] = FitPsycheCurveLogit_cgs(uniqueFrequencies, y, ones(size(uniqueFrequencies)), targets); %FOR GRAPH ONLY
[~, ~, fit_curve2, fit_threshold2] = FitPsycheCurveLogit_cgs(uniqueFrequencies, yy, weights, targets); %FOR GRAPH ONLY
[coeffs, fit_curve3] = FitPsycheCurveWH(uniqueFrequencies, y, SPs); %FOR GRAPH ONLY
fit_threshold3 = nan(1,length(targets));
for f = 1:length(targets)
    [~, min_ind] = min(abs(fit_curve3(:,2) - targets(f)));
    fit_threshold3(f) = fit_curve3(min_ind,1);
end

% THIS DOESNT WORK BECAUSE FIT IS NOT STRETCHED TO LINEAR AXIS
% %Change curve x from octave to linear scale
% fit_curve_percent = fit_curve3(:,1)./max(fit_curve3(:,1));
% fit_curve_linear = 1 + fit_curve_percent*(length(uniqueFrequencies)-1);

% Plot psychometic curves
%plot(fit_curve(:,1), fit_curve(:,2), 'Linewidth', 2, 'LineStyle', '--', 'Color', 'g')
%plot(fit_curve2(:,1), fit_curve2(:,2)./100, 'Linewidth', 2, 'LineStyle', '--', 'Color', 'm')    
plot(fit_curve3(:,1), fit_curve3(:,2), 'Linewidth', 2, 'LineStyle', '--', 'Color', 'g')    
%scatter(fit_threshold, targets, 'x', 'g')
%scatter(fit_threshold2, targets, 'x', 'm')
scatter(fit_threshold3, targets, 'x', 'k')
title(['Threshold = ' num2str(fit_threshold3(2))])

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
yy = reaction_timesPerFrequency;
upper = yy + reaction_timesPerFrequency_STD;
lower = yy - reaction_timesPerFrequency_STD;

%make figure
%figure; 
subplot(3,1,2); hold on
bar(x,yy)
for t = 1:length(x)
    text(x(t)-0.2,  0.9, num2str(nRepsPerFrequency(t)))
end
%errorbar(x,y',lower,upper) %Plot standard deviation error bars
scatterx = nan(1,length(reaction_times));
for i = 1:length(uniqueFrequencies)
    scatterx(allFrequencies == uniqueFrequencies(i)) = i;
end
scatter(scatterx,reaction_times, 'k')

ylabel('Reaction Time')
%xlabel('Frequency difference in octaves')
set(gca, 'XTick', x)
set(gca, 'XTickLabel', uniqueFrequencies)
ylim([0 1])
xlim([x(1)-1 x(end)+1])

% Plot psychometric curve with bins of size 2

% %set up axes
% A = hitRatePerFrequency;
% y = mean([A(1:2:end-1);A(2:2:end)]);
% x = 1:length(y);
% 
% %make figure
% subplot(3,2,5); hold all
% plot(x,y,'Linewidth',2) %plot psychometric curve
% line([0 x(end)+1], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
% ylabel('Hit Rate')
% xlabel('Frequency difference in octaves binned by 2')
% set(gca, 'XTick', x)
% set(gca, 'XTickLabel', uniqueFrequencies(2:2:end))
% 
% % Fit psychometric functions
% targets = [0.25 0.5 0.75]; % 25 50 75 performance
% weights = ones(1,length(y)); % No weighting
% [fit_coeffs, fit_curve_bin2, fit_threshold_bin2] = fitPsycheCurveLogit(x, y, weights, targets);
% 
% % Plot psychometic curves
% plot(fit_curve_bin2(:,1), fit_curve_bin2(:,2), 'Linewidth', 2, 'LineStyle', '--', 'Color', 'g')
% %legend('Performance', 'Fit');
% scatter(fit_threshold_bin2, targets, 'x', 'k')
% 
% [~, ~, actual_threshold_bin2] = fitPsycheCurveLogit(uniqueFrequencies(2:2:end), y, weights, targets)

suptitle(plotTitle)

end

%% Plot 3 - Locomotion and licks

%Get locomotor activity
loco_activity = block.loco_activity;
loco_times = block.loco_times;

ymax = round(1.05*max(loco_activity)); %Determine how high the graph ylim should be

%Plot entire block
figure; hold on
plot(loco_times,loco_activity) %Loco activity is a line graph
xlim([0 loco_times(end)])
ylim([0 ymax])
xlabel('Time (s)')
ylabel('Running speed (cm/s)')
title(plotTitle)

% %Plot a black vertical line for every sound presentation
% for i = 1:length(soundTimes)
%     xsound = soundTimes(i);
%     line([xsound xsound], [0, ymax], 'Color', 'k')
% end

%Plot a green or red vertical line (target/nontarget) for each response windows
%Plot licks and water delivery
alternating_sound_times = soundTimes_removeEL + holdingPeriod;
trial_type_removeEL = trial_type(~earlyLicks);
lick_times = block.concat_times(block.concat_licks == 1);
water_times = [];
for i = 1:length(alternating_sound_times) 
    xsound = alternating_sound_times(i);
    if trial_type_removeEL(i) == 1 %Target
        lineColor = 'g';
        if hits(i) == 1
            water_ind = find(lick_times > xsound, 1, 'first');
            water_times = [water_times, lick_times(water_ind)];
        end
    elseif trial_type_removeEL(i) == 0 %Non-target
        lineColor = 'r';
    end
        
    line([xsound xsound], [0, ymax], 'Color', lineColor)
end

%Plot a circle for every lick
scatter(lick_times,zeros(size(lick_times))+ymax,'m','Filled')
%Plot a circle for water
scatter(water_times,zeros(size(water_times))+ymax,'c','Filled')

%% Plot example trial(s)

plotExampleTrial = 0;

if plotExampleTrial
    
%Adjust xlim and ylim to show only a selection of trials
nTrialsToPlot = 3;
trialToStartWith = 5;

xStart = soundTimes(trialToStartWith);
xEnd = soundTimes(trialToStartWith + nTrialsToPlot);
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
    if trial_type_removeEL(i) == 1 %Target
        lineColor = 'g';
    elseif trial_type_removeEL(i) == 0 %Non-target
        lineColor = 'r';
    end
        
    line([xsound xsound], [0, ymax], 'Color', lineColor)
end

%Plot a circle for every lick
scatter(lick_times,zeros(size(lick_times))+ymax,'m','Filled')
%Plot a circle for water
scatter(water_times,zeros(size(water_times))+ymax,'c','Filled')
end

%% Plot #4 - Plot lick rasterplot

[sortedFreqs, sortInd] = sort(allFrequencies);
sortedTrialRaster = lickRaster(sortInd,:);
sortedHits = hits(sortInd);

figure
subplot(2,1,1)
imagesc(flipud(sortedTrialRaster))
set(gca, 'YTick', 1:length(sortedFreqs))
set(gca, 'YTickLabel', fliplr(sortedFreqs));
ylabel('Frequency')
xlabel('Frames aligned to response window')
title(plotTitle)
vline(baselinePeriodInFrames, 'w')

subplot(2,1,2); hold all %scatterplot
for i = 1:size(sortedTrialRaster,1)
    currentPoints = find(sortedTrialRaster(i,:));
    scatter(currentPoints,zeros(1,length(currentPoints)) + i, 50, 'm', 'filled');
    %Determine water delivery
    validPoints = currentPoints(currentPoints > baselinePeriodInFrames);
    if ~isempty(validPoints) && sortedHits(i) == 1
        firstPoint = validPoints(1);
        if firstPoint <= validResponsePeriodInFrames + baselinePeriodInFrames
            scatter(firstPoint,zeros(1,length(firstPoint)) + i, 50, 'b', 'filled');
        end
    end
end
xlim([1 totalPeriodInFrames])
vline(baselinePeriodInFrames, 'k')
set(gca, 'YTick', 1:length(sortedFreqs))
set(gca, 'YTickLabel', sortedFreqs);
ylabel('Frequency')
xlabel('Frames aligned to response window')

%% Plot frequency of hits and/or FPs over time (to see if mouse activity changes over a single session)

FPtimes = soundTimes(FPs);

figure;
subplot(2,1,1)
hist(FPtimes);
ylabel('Number of FPs')
xlabel('Time (seconds)')
xlim([0 loco_times(end)])
title(plotTitle)

subplot(2,1,2)
scatter(FPtimes, 1:length(FPtimes))
ylabel('FP count')
xlabel('Time (seconds)')
xlim([0 loco_times(end)])

%% Plot mouse's lick times during early lick trials

% if sum(earlyLicks) > 0
%     
% soundTimes_onlyEL = soundTimes(earlyLicks);
% earlyLick_times = nan(size(soundTimes_onlyEL));
% for i = 1:length(soundTimes_onlyEL)
%     earlyLick_ind = find(lick_times > soundTimes_onlyEL(i), 1, 'first');
%     if isempty(earlyLick_ind)
%         continue
%     end
%     earlyLick_times(i) = lick_times(earlyLick_ind) - soundTimes_onlyEL(i);
% end
% 
% figure;
% hist(earlyLick_times)
% title(['N early licks = ' num2str(length(earlyLick_times))])
% xlabel('Time after sound start (frames)')
% ylabel('N early licks')
% suptitle(plotTitle)
% 
% end

%% Plot mouse's lick times during holding period

%if sum(earlyLicks) == 0
    
%Convert soundTimes to indices
soundTimes_ind = nan(size(soundTimes_removeEL)); 
for i = 1:length(soundTimes_removeEL)
    [~, soundTimes_ind(i)] = min(abs(block.concat_times - soundTimes_removeEL(i)));
end

%How much time to look at
baselinePeriodInSeconds = 1; %seconds before sound start
responsePeriodInSeconds = 25; %seconds after sound start

baselinePeriodInFrames = round(baselinePeriodInSeconds/trialFrameRate);
responsePeriodInFrames = round(responsePeriodInSeconds/trialFrameRate);
totalPeriodInFrames = baselinePeriodInFrames + responsePeriodInFrames;

holdingPeriodRaster = zeros(length(soundTimes_ind),totalPeriodInFrames);
responseWindowStart_zeroed = nan(1,length(responseWindowStart));
for i = 1:length(soundTimes_ind)
    sound_start = soundTimes_ind(i);
    a = 1;
    b = totalPeriodInFrames;
    A = sound_start - baselinePeriodInFrames;
    B = sound_start + responsePeriodInFrames - 1;
    responseWindowStart_zeroed(i) = responseWindowStart(i) - sound_start;
    
    %instance where A is before block starts
    if A <= 0
        a = -A + 2;
        A = 1;
    end
    
    %instance where B is after block ends
    if B > length(block.concat_licks)
        B = length(block.concat_licks);
        b = length(A:B);
    end
       
    holdingPeriodRaster(i,a:b) = block.concat_licks(A:B,1);
end

% FIGURES

[sorted_responseWindows, sort_RW_ind] = sort(responseWindowStart_zeroed);
sorted_holdingPeriodRaster = holdingPeriodRaster(sort_RW_ind,:);

figure;
subplot(3,1,1); hold all
imagesc(sorted_holdingPeriodRaster)
ylim([1 nTrials])
xlim([1 totalPeriodInFrames])
xlabel('Time after sound start (frames)')
ylabel('Sorted trial')
title('Lick raster (full trials)')

for i = 1:length(sorted_responseWindows)
    xstart = sorted_responseWindows(i) + baselinePeriodInFrames;
    line([xstart xstart], [i - 0.5, i + 0.5], 'Color', 'r')
end

subplot(3,1,2)
area(sum(sorted_holdingPeriodRaster))
ylabel('Licks')
xlabel('Time after sound start (frames)')
xlim([1 totalPeriodInFrames])
title('Lick PSTH')

subplot(3,1,3)
hist(holdingPeriod)
xlabel('Hold duration (s)')
ylabel('N trials')
title('Holding period histogram')

suptitle(plotTitle)

%end