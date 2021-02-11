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
    
    %defining new repeating/alternating tones
    repeating = unique(stim_v1_loco); 
    alternating = unique(stim_v2_loco);
    
    figure;
    %pull out only 1 alternating frequency at a time
    maxvals = nan(1,length(alternating));
    for a = 1:length(alternating)
        currentFrequency = alternating(a);
        currentFrequencyIndex = stim_v2_loco == currentFrequency; 
        F7_freq = F7_df_f(currentFrequencyIndex,:);

        %Average data
        avg_F7_freq = nanmean(F7_freq);

        %Add to list of cellResponses

        y = smooth(avg_F7_freq,3)'; %Average and smooth responses
        subplot(3,ceil(length(alternating)/3),a)
        plot(y)
        vline(nBaselineFrames, '-') %plot vertical red line at sound onset
        xlim([0 size(y,2)])
        ylabel('DF/F')
        xlabel('Time (frames)')
        title([num2str(currentFrequency)])
        axis square
        
        maxvals(a) = max(y);
    end
    
    %Adjust max of each plot
    for a = 1:length(alternating)
        subplot(3,ceil(length(alternating)/3),a)
        ylim([0 quantile(maxvals,0.9)])
    end
    
    suptitle(['cell #' num2str(cellNumber)])
    
    figure
    plot(maxvals)
    ylabel('Max response')
    xlabel('Alternating frequency')
    set(gca,'XTick',1:2:length(alternating))
    set(gca,'XTickLabel',num2str(alternating(1:2:end)))
    suptitle(['cell #' num2str(cellNumber)])

end

%% Plot "neurometric" curve
% Look at the Psychometric Curve section in the code plot_behavior_maryse_and_esther
% Can you plot a "neurometric" curve for each cell using the peak amplitude of its response to the alternating frequencies?


% %Calculate hit rate per frequency
% hitsPerFrequency = nan(1,length(uniqueFrequencies)); %make empty vector to fill with data
% for i = 1:length(uniqueFrequencies)
%     hits_and_FPs = hits + FPs;
%     freqResponse = hits_and_FPs(allFrequencies == uniqueFrequencies(i));
%     hitsPerFrequency(i) = sum(freqResponse);
% end
% hitRatePerFrequency = hitsPerFrequency./nRepsPerFrequency;
% 
% %set up axes
% x = 1:length(hitRatePerFrequency);
% y = hitRatePerFrequency;
% 
% %make figure
% figure;
% subplot(2,1,1); hold all
% plot(x,y,'Linewidth',2) %plot psychometric curve
% line([0 x(end)+1], [0.5, 0.5], 'Color', 'r') %plot horizontal red line at 50% hit rate (chance)
% ylabel('Hit Rate')
% xlabel('Alternating frequency')
% set(gca, 'XTick', x)
% set(gca, 'XTickLabel', uniqueFrequencies)
% title(plotTitle)
% 
% % Fit psychometric functions
% targets = [0.25 0.5 0.75]; % 25 50 75 performance
% weights = ones(1,length(y)); % No weighting
% [fit_coeffs, fit_curve, fit_threshold] = fitPsycheCurveLogit(x, y, weights, targets);
% 
% % Plot psychometic curves
% plot(fit_curve(:,1), fit_curve(:,2), 'Linewidth', 2, 'LineStyle', '--', 'Color', 'g')
% %legend('Performance', 'Fit');
% scatter(fit_threshold, targets, 'x', 'k')

% 
% %% Load block and assign variables
% 
% % %Load block file
% % disp('Load block file')
% % [block_filename,block_filepath] = uigetfile('.mat');
% % load(block_filename)
% 
% %Block data
% mousename = char(block.setup.mousename);
% baseline_length = block.setup.constant.baseline_length; %Amount of time before sound onset in seconds
% framerate = block.setup.framerate; %Frame rate of the 2p recording (# frames per second)
% nBaselineFrames = baseline_length*framerate; %Amount of time before sound onset in frames
% response_window = block.setup.constant.response_window; %Time where we expect a response to sound onset in seconds
% response_window_in_frames = response_window*framerate; %Same as above in frames
% 
% %Stimulus data
% stim_v1 = block.parameters.variable1'; %target frequency per trial
% stim_v2 = block.parameters.variable2'; %nontarget (alternating) frequency per trial
% %holdingPeriod = block.holdingPeriod'; %holdiner period per trial
% target = unique(stim_v1);
% nontargets = unique(stim_v2);
% 
% %% Find sound-responsive neurons
% 
% STDlevel = 2; %# of standard deviations above baseline for a response to be considered active
% AUC_F_level = 50; %Minimum area under the curve for a response to be considered active
% sort_active = 1;  % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials
% 
% % identify +/- locomotor trials to remove
% stim_v1_loco = stim_v1;
% stim_v2_loco = stim_v2;
% remove = [];
% 
% if sort_active ~= 0
%     if sort_active==1
%         remove = find(block.active_trials==1);%find active trials to remove
%     elseif sort_active==2
%         remove = find(block.active_trials==0);%find inactive trials to remove
%     end
%     stim_v1_loco(remove,:) = [];
%     stim_v2_loco(remove,:) = [];
% end
% 
% % Go through each cell and determine whether it is sound-responsive or not
% cells = block.cell_number; %List of cell numbers determined by suite2p
% responsiveCells = zeros(size(cells)); %We'll record which cells are responsive here
% activityType = cell(size(cells));
% 
% for c = 1:length(cells)
%     cellNumber = cells(c); %Suite2p label
%     cellIndex = find(block.cell_number == cellNumber); %Row number
%     
%     %block.aligned_stim contains all of the fluorescence data for each neuron
%     %It is arranged as N cells x trials x frames
%     F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:)); %Get data for our cell only
%     F7_baseline = F7(:,1:nBaselineFrames); %Find the baseline responses
%     F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %Get DF/F = [response-mean(baseline)]/mean(baseline)
% 
%     %Remove +/- loco trials
%     F7_df_f(remove,:) = [];
% 
%     %Average and smooth data
%     avg_F7_df_f = smooth(nanmean(F7_df_f),3)';
%     
%     %Check if cell is active
%     [active, tempActivity] = checkIfActive(avg_F7_df_f, nBaselineFrames, STDlevel, AUC_F_level, 0);
%     if active
%         responsiveCells(c) = 1;
%         activityType{c} = tempActivity;
%     end
% end
% 
% %% Plot average response to all stimuli for all cells
% 
% cellResponsesToPlot = []; %Store responses here
% 
% for c = 1:length(cells)
%     %Skip this cell if not responsive/activated
%     if responsiveCells(c) ~= 1
%         continue
%     elseif strcmp('inhibited', activityType{c})
%         continue
%     end
%     cellNumber = cells(c); %Suite2p label
%     cellIndex = find(block.cell_number == cellNumber); %Row number    
%     
%     %Get fluorescence data
%     F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:)); %Get data for our cell only
%     F7_baseline = F7(:,1:nBaselineFrames); %Find the baseline responses
%     F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %Get DF/F = [response-mean(baseline)]/mean(baseline)
% 
%     %Remove +/- loco trials
%     F7_df_f(remove,:) = [];
% 
%     %Average data
%     avg_F7_df_f = nanmean(F7_df_f);
%     
%     %Add to list of cellResponses
%     cellResponsesToPlot = [cellResponsesToPlot; avg_F7_df_f];
% end
% 
% %Plot graph
% y = smooth(nanmean(cellResponsesToPlot),3)'; %Average and smooth responses
% 
% figure
% plot(y)
% vline(nBaselineFrames, '-') %plot vertical red line at sound onset
% xlim([0 size(y,2)])
% ylabel('DF/F')
% xlabel('Time (frames)')
% title(['Average of ' num2str(size(cellResponsesToPlot,1)) ' cells'])
% 
% %% Plot average response to all stimuli per cell
% % Hint - copy the above script but put the plot inside the for loop
% 
% 
% %% Plot average response by stimulus type per cell
% % Use stim_v1_loco and stim_v2_loco to extract one stim type at a time from block.aligned_stim
% 
% 
% %% Plot "neurometric" curve
% % Look at the Psychometric Curve section in the code plot_behavior_maryse_and_esther
% % Can you plot a "neurometric" curve for each cell using the peak amplitude of
% % its response to the alternating frequencies?

