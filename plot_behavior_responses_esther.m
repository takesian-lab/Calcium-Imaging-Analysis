%plot_behavior_responses_esther

%% Load block and assign variables

%Load block file
% disp('Load block file')
% [block_filename,block_filepath] = uigetfile('.mat');
% load(block_filename)

%Block data
mousename = char(block.setup.mousename);
baseline_length = block.setup.constant.baseline_length; %Amount of time before sound onset in seconds
framerate = block.setup.framerate; %Frame rate of the 2p recording (# frames per second)
nBaselineFrames = baseline_length*framerate; %Amount of time before sound onset in frames
response_window = block.setup.constant.response_window; %Time where we expect a response to sound onset in seconds
response_window_in_frames = response_window*framerate; %Same as above in frames
sound_duration = 3; %seconds
sound_duration_in_frames = sound_duration*framerate;

%Stimulus data

stimInOv = 0;

if stimInOv
    stim_v1_kHz = unique(block.parameters.variable1');
    stim_v2_kHz = unique(block.parameters.variable2');
    stim_v1 = log2(block.parameters.variable1'); %repeating frequency per trial
    stim_v2 = log2(block.parameters.variable2'); %alternating frequency per trial
    holdingPeriod = block.holdingPeriod'; %holdiner period per trial
    repeating_log = unique(stim_v1);
    stim_v1 = stim_v1 - repeating_log;
    stim_v2 = round(stim_v2 - repeating_log,2);
    alternating = unique(stim_v2);
    repeating = unique(stim_v1);
else
    stim_v1 = block.parameters.variable1'; %repeating frequency per trial
    stim_v2 = block.parameters.variable2'; %alternating frequency per trial
    holdingPeriod = block.holdingPeriod'; %holdiner period per trial
    alternating = unique(stim_v2);
    repeating = unique(stim_v1);
end
%% Find sound-responsive neurons

STDlevel = 2; %# of standard deviations above baseline for a response to be considered active
AUC_F_level = 50; %Minimum area under the curve for a response to be considered active
sort_active = 0;  % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials

% identify +/- locomotor trials to remove
stim_v1_loco = stim_v1;
stim_v2_loco = stim_v2;
remove = [];

if sort_active ~= 0
    if sort_active==1
        remove = find(block.active_trials==1);%find active trials to remove
    elseif sort_active==2
        remove = find(block.active_trials==0);%find inactive trials to remove
    end
    stim_v1_loco(remove,:) = [];
    stim_v2_loco(remove,:) = [];
end

% Go through each cell and determine whether it is sound-responsive or not
cells = block.cell_number; %List of cell numbers determined by suite2p
responsiveCells = zeros(size(cells)); %We'll record which cells are responsive here
activityType = cell(size(cells));

for c = 1:length(cells)
    cellNumber = cells(c); %Suite2p label
    cellIndex = find(block.cell_number == cellNumber); %Row number
    
    %block.aligned_stim contains all of the fluorescence data for each neuron
    %It is arranged as N cells x trials x frames
    F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:)); %Get data for our cell only
    F7_baseline = F7(:,1:nBaselineFrames); %Find the baseline responses
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %Get DF/F = [response-mean(baseline)]/mean(baseline)

    %Remove +/- loco trials
    F7_df_f(remove,:) = [];

    %Average and smooth data
    avg_F7_df_f = smooth(nanmean(F7_df_f),3)';
    
    %Check if cell is active
    [active, tempActivity] = checkIfActive(avg_F7_df_f, nBaselineFrames, STDlevel, AUC_F_level, 0);
    if active
        responsiveCells(c) = 1;
        activityType{c} = tempActivity;
    end
end

%% Plot average response to all stimuli for all cells

cellResponsesToPlot = []; %Store responses here

for c = 1:length(cells)
    %Skip this cell if not responsive/activated
    if responsiveCells(c) ~= 1
        continue
    elseif strcmp('inhibited', activityType{c})
        continue
    end
    cellNumber = cells(c); %Suite2p label
    cellIndex = find(block.cell_number == cellNumber); %Row number    
    
    %Get fluorescence data
    F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:)); %Get data for our cell only
    F7_baseline = F7(:,1:nBaselineFrames); %Find the baseline responses
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %Get DF/F = [response-mean(baseline)]/mean(baseline)

    %Remove +/- loco trials
    F7_df_f(remove,:) = [];

    %Average data
    avg_F7_df_f = nanmean(F7_df_f);
    
    %Add to list of cellResponses
    cellResponsesToPlot = [cellResponsesToPlot; avg_F7_df_f];
end

%Plot graph
y = smooth(nanmean(cellResponsesToPlot),3)'; %Average and smooth responses

figure
plot(y)
vline(nBaselineFrames, '-') %plot vertical red line at sound onset
xlim([0 size(y,2)])
ylabel('DF/F')
xlabel('Time (frames)')
title(['Average of ' num2str(size(cellResponsesToPlot,1)) ' cells'])


%% Plot average response to all stimuli per cell (START HERE)
% Hint - copy the above script but put the plot inside the for loop

for c = 1:length(cells)
    %Skip this cell if not responsive/activated
    if responsiveCells(c) ~= 1
        continue
    elseif strcmp('inhibited', activityType{c})
        continue
    end
    cellNumber = cells(c); %Suite2p label
    cellIndex = find(block.cell_number == cellNumber); %Row number    
    
    %Get fluorescence data
    F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:)); %Get data for our cell only
    F7_baseline = F7(:,1:nBaselineFrames); %Find the baseline responses
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %Get DF/F = [response-mean(baseline)]/mean(baseline)

    %Remove +/- loco trials
    F7_df_f(remove,:) = [];

    %Average data
    avg_F7_df_f = nanmean(F7_df_f);
    
    %Add to list of cellResponses
   
    y = smooth(avg_F7_df_f,3)'; %Average and smooth responses
    figure
    plot(y)
    vline(nBaselineFrames, '-') %plot vertical red line at sound onset
    xlim([0 size(y,2)])
    ylabel('DF/F')
    xlabel('Time (frames)')
    title(['Average of cell #' num2str(cellNumber)])
end

%% Plot average response to stimulus type per cell
%% Use stim_v1_loco and stim_v2_loco to extract from block.aligned_stim
%% This is what I need to work on for HCRP final report (change title/legend)

cellResponsesToPlot = []; %Store responses here

store_cellNumber = nan(sum(responsiveCells),1);
neurometric_curve_max = nan(sum(responsiveCells),length(alternating));
neurometric_curve_mean = nan(sum(responsiveCells),length(alternating));
store_modIndexMax = nan(sum(responsiveCells),1);
store_modIndexMean = nan(sum(responsiveCells),1);
count = 1;

for c = 1:length(cells)
    %Skip this cell if not responsive/activated
    if responsiveCells(c) ~= 1
        continue
    elseif strcmp('inhibited', activityType{c})
        continue
    end
    cellNumber = cells(c); %Suite2p label
    store_cellNumber(count,1) = cellNumber;
    cellIndex = find(block.cell_number == cellNumber); %Row number    
    
    %Get fluorescence data
    F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:)); %Get data for our cell only
    F7_baseline = F7(:,1:nBaselineFrames); %Find the baseline responses
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %Get DF/F = [response-mean(baseline)]/mean(baseline)

    %Remove +/- loco trials
    F7_df_f(remove,:) = [];
    
    %defining new repeating/alternating tones
    repeating = unique(stim_v1_loco); 
    alternating = unique(stim_v2_loco);
    
    figure;
    %pull out only 1 alternating frequency at a time
    maxvals = nan(1,length(alternating));
    meanvals = nan(1,length(alternating));
    for a = 1:length(alternating)
        currentFrequency = alternating(a);
        currentFrequencyIndex = stim_v2_loco == currentFrequency; 
        F7_freq = F7_df_f(currentFrequencyIndex,:);

        %Average data
        avg_F7_freq = nanmean(F7_freq);

        %Add to list of cellResponses

        y = smooth(avg_F7_freq,3)'; %Average and smooth responses
        subplot(3,length(alternating),a)
        plot(y)
        vline(nBaselineFrames, '-') %plot vertical red line at sound onset
        xlim([0 size(y,2)])
        ylabel('DF/F')
        xlabel('Time (frames)')
        title([num2str(currentFrequency)])
        axis square
        
        maxvals(a) = max(y);
        meanvals(a) = nanmean(y(1,nBaselineFrames:nBaselineFrames+sound_duration_in_frames));
        
        if isnan(maxvals(a)) || isnan(meanvals(a))
            error('nan')
        end
    end
    
    neurometric_curve_max(count,:) = maxvals;
    neurometric_curve_mean(count,:) = meanvals;
    
    %Adjust max of each plot
    for a = 1:length(alternating)
        subplot(3,length(alternating),a)
        ylim([0 quantile(maxvals,0.9)])
    end

    %Modulation index
    modIndexMax = (nanmean(maxvals(1:3))-nanmean(maxvals(end-2:end)))/(nanmean(maxvals(1:3)) + nanmean(maxvals(end-2:end)));
    modIndexMean = (nanmean(meanvals(1:3))-nanmean(meanvals(end-2:end)))/(nanmean(meanvals(1:3)) + nanmean(meanvals(end-2:end)));
    
    store_modIndexMax(count,1) = modIndexMax;
    store_modIndexMean(count,1) = modIndexMean;
    
    %Plot neurometric curve
    subplot(3,length(alternating),[1+length(alternating) length(alternating)*2]); hold on
    plot(maxvals)
    ylabel('Max response')
    xlabel('Alternating frequency')
    set(gca,'XTick',1:2:length(alternating))
    set(gca,'XTickLabel',num2str(alternating(1:2:end)))
    title(['Modulation index: ' num2str(modIndexMax)])
    
    %Fit gaussian
    options = fitoptions('gauss1');
    options.Lower = [0 1 0];
    options.Upper = [inf 80 80];

    gauss_fit = fit((1:length(maxvals))', maxvals', 'gauss1',options);
    gauss_curve = gauss_fit(1:0.1:length(maxvals));
    gauss_width = gauss_fit.c1*2;
    [B,BI] = max(gauss_fit(1:length(maxvals))); 
    best_freq = alternating(BI);
    plot(1:0.1:length(maxvals),gauss_curve, 'r');

    subplot(3,length(alternating),[1+(2*length(alternating)) length(alternating)*3]); hold on
    plot(meanvals)
    ylabel('Mean response')
    xlabel('Alternating frequency')
    set(gca,'XTick',1:2:length(alternating))
    set(gca,'XTickLabel',num2str(alternating(1:2:end)))
    title(['Modulation index: ' num2str(modIndexMean)])  
    
    %Fit gaussian
    options = fitoptions('gauss1');
    options.Lower = [0 1 0];
    options.Upper = [inf 80 80];

    gauss_fit = fit((1:length(maxvals))', meanvals', 'gauss1',options);
    gauss_curve = gauss_fit(1:0.1:length(maxvals));
    gauss_width = gauss_fit.c1*2;
    [B,BI] = max(gauss_fit(1:length(maxvals))); 
    best_freq = alternating(BI);
    plot(1:0.1:length(maxvals),gauss_curve, 'r');
    
   suptitle(['cell #' num2str(cellNumber)])

   count = count + 1;
end

%%
data = struct;
data.mousename = mousename;
data.block = block.setup.block_filename;
data.STDlevel = STDlevel;
data.AUC_F_level = AUC_F_level;
data.sort_active = sort_active;
data.cell_number = nan(size(cells));
data.neurometric_curve_max = neurometric_curve_max;
data.neurometric_curve_mean = neurometric_curve_mean;
data.modIndexMax = store_modIndexMax;
data.modIndexMean = store_modIndexMean;

%% Figure

normMax = data.neurometric_curve_max./max(data.neurometric_curve_max,[],2);
normMean = data.neurometric_curve_mean./max(data.neurometric_curve_mean,[],2);

[~, max_I] = sort(data.modIndexMax);
sorted_max = data.neurometric_curve_max(max_I,:);
sorted_max = normMax(max_I,:);
[~, mean_I] = sort(data.modIndexMean);
sorted_mean = data.neurometric_curve_mean(mean_I,:);
sorted_mean = normMean(mean_I,:);


figure;

subplot(2,2,1)
imagesc(sorted_max)
ylabel('Sorted cells')
set(gca, 'XTick', 1:length(alternating))
set(gca, 'XTickLabel', alternating)
c = colorbar;
title('Max neurometric curve')
ylabel(c, 'Max df/f response')

subplot(2,2,3); hold on
scatter(ones(1,length(data.modIndexMax)), data.modIndexMax)
ylim([-1 1])
hline(0)
hline(nanmean(data.modIndexMax), 'k')
title(['Average mod index: ' num2str(nanmean(data.modIndexMax))])
ylabel('Modulation index')

subplot(2,2,2)
imagesc(sorted_mean)
ylabel('Sorted cells')
set(gca, 'XTick', 1:length(alternating))
set(gca, 'XTickLabel', alternating)
c = colorbar;
title('Mean neurometric curve')
ylabel(c,'Mean df/f response')

subplot(2,2,4); hold on
scatter(ones(1,length(data.modIndexMean)), data.modIndexMean)
ylim([-1 1])
hline(0)
hline(nanmean(data.modIndexMean), 'k')
title(['Average mod index: ' num2str(nanmean(data.modIndexMean))])
ylabel('Modulation index')


