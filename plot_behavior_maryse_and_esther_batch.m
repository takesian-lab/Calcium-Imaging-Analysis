%plot_behavior_maryse_and_esther_batch

%% Load files

stim_protocol = 13; %Maryse behavior
stim_name = 'training'; %training or testing or behavior

saveData = 0;
loadPreviousData = 0;

if loadPreviousData
    %Load data
    [FileName,PathName] = uigetfile('.mat');
    load([PathName '/' FileName])
else
    PC_name = getenv('computername');

    switch PC_name
        case 'RD0366' %Maryse
            info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
            save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
            compiled_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\Compiled Blocks v2';
            info_filename = 'Info v2';    
        otherwise
            disp('Computer does not match known users')
            return
    end

    cd(info_path)
    Info = importfile(info_filename);

    [data] = fillSetupFromInfoTable_v3(Info, compiled_blocks_path, stim_protocol, stim_name);
    disp('Data structure created.')
end

%% Prepare variables

mice = unique(data.setup.mousename(:,1));
nMice = length(mice);
dates = cell(1,nMice);
nDatesPerMouse = nan(1,nMice);

for i = 1:nMice
    mouseInd = strcmp(data.setup.mousename, mice{i});
    dates{i} = unique([data.setup.expt_date{mouseInd}]);
    nDatesPerMouse(i) = length(dates{i});
end
    
    batch_hitRate = nan(nMice, max(nDatesPerMouse));
    batch_FPrate = nan(nMice, max(nDatesPerMouse));    
if strcmp(stim_name, 'testing')    ||  strcmp(stim_name, 'behavior') 
    batch_freqs = cell(nMice, max(nDatesPerMouse));
    batch_hitRatePerFreq = cell(nMice, max(nDatesPerMouse)); 
    batch_rxnTimesPerFreq = cell(nMice, max(nDatesPerMouse)); 
end

%% Extract data from blocks

for a = 1:nMice %Mice
    for b = 1:nDatesPerMouse(a) %Training days
        
       mouseID = mice{a};
       currentDate = dates{a}(b);
       dateInd = strcmp([data.setup.expt_date{a,:}], currentDate);
       matchingBlocks = data.setup.unique_block_names(a,dateInd);
       
       %Prepare variables to combine blocks together
       combinedFreqs = [];
       combinedHits = [];
       combinedMisses = [];
       combinedFPs = [];
       combinedWitholds = [];
       combinedRxnTimes = [];

       %Make complicated condition to only include 40dB blocks if they're
       %the only data for that day
       stimLevels = nan(1,length(matchingBlocks));
       for c = 1:length(matchingBlocks)
           unique_block_name = matchingBlocks{c};
           block = data.([mouseID]).([unique_block_name]);
           stimLevels(c) = block.stim_level;
       end
       
       if unique(stimLevels) == 40
           %continue;
       else
           matchingBlocks(stimLevels == 40) = [];
       end
           
       for c = 1:length(matchingBlocks)
           unique_block_name = matchingBlocks{c};
           block = data.([mouseID]).([unique_block_name]);

           %% Get information from block
           
            %Block data
            trial_type = block.trialType;
            outcomes = block.Outcome;
            earlyLicks = isnan(outcomes);
            outcomes_removeEL = outcomes(~earlyLicks);
            nTrials = length(outcomes_removeEL);
            misses      = outcomes_removeEL == 0;
            hits        = outcomes_removeEL == 1;
            witholds    = outcomes_removeEL == 3;
            FPs         = outcomes_removeEL == 4;

            allFrequencies = log2(block.TargetFreq(~earlyLicks)); %in log scale
            repeatingFrequency = mode(allFrequencies);
            allFrequencies = round(abs(allFrequencies - repeatingFrequency),2); %in octave difference
            holdingPeriod = block.holdingPeriod(~earlyLicks);

            %I removed the holding period in later versions of the task
            noHoldingPeriod = false;
            if length(unique(holdingPeriod)) == 1
                if unique(holdingPeriod) <= 0.2
                    noHoldingPeriod = true;
                    holdingPeriod(:) = 0;
                end
            end

            %Reaction times from the start of the sound
            raw_reactionTimes = block.rxn_time(~earlyLicks);
            raw_reactionTimes(raw_reactionTimes < 0) = nan; %trials with no responses become NaN
            raw_reactionTimes = raw_reactionTimes/1000; %convert to seconds
            if noHoldingPeriod
                reaction_times = raw_reactionTimes;
            else
                reaction_times = raw_reactionTimes - holdingPeriod; %subtract holding period
            end
            
            %Store data to combine with other blocks
            combinedFreqs    = [combinedFreqs, allFrequencies];
            combinedHits     = [combinedHits, hits];
            combinedMisses   = [combinedMisses, misses];
            combinedFPs      = [combinedFPs, FPs];
            combinedWitholds = [combinedWitholds, witholds];
            combinedRxnTimes = [combinedRxnTimes, reaction_times]; 
            
       end
       
        %Training data
        batch_hitRate(a,b) = sum(combinedHits)/(sum(combinedHits) + sum(combinedMisses));
        batch_FPrate(a,b) = sum(combinedFPs)/(sum(combinedFPs) + sum(combinedWitholds));
        
        %Pscyhometric data
        if strcmp(stim_name, 'testing') || strcmp(stim_name, 'behavior')
            if isempty(combinedFreqs) %No 60dB data
                continue
            end
            table = tabulate(combinedFreqs);
            uniqueFrequencies = table(:,1)';
            nRepsPerFrequency = table(:,2)';
            hitsPerFrequency = nan(1,length(uniqueFrequencies)); 
            rxnTimesPerFrequency = nan(1,length(uniqueFrequencies)); 
            for i = 1:length(uniqueFrequencies)
                hits_and_FPs = combinedHits + combinedFPs;
                freqResponse = hits_and_FPs(allFrequencies == uniqueFrequencies(i));
                hitsPerFrequency(i) = sum(freqResponse);
                rxnTimesPerFrequency(i) = nanmean(combinedRxnTimes(allFrequencies == uniqueFrequencies(i)));
            end
            hitRatePerFrequency = hitsPerFrequency./nRepsPerFrequency;

            batch_freqs{a,b} = uniqueFrequencies;
            batch_hitRatePerFreq{a,b} = hitRatePerFrequency; 
            batch_rxnTimesPerFreq{a,b} = rxnTimesPerFrequency;
        end
    end
end

%%  PLOTS

if strcmp(stim_name, 'training')
    
    %HIT RATE AND FP RATE
    figure;
    
    %subplot(1,2,1)
    plot(batch_hitRate'); hold on
    ylim([0 1])
    ylabel('Percent')
    xlabel('Training day')
    %title('Hit rate')
    legend(mice, 'Location', 'southwest')
    
    %subplot(1,2,2)
    plot(batch_FPrate')
    ylim([0 1])
    xlabel('Training day')
    title('FP rate')
    legend(mice, 'Location', 'southwest')
    
elseif strcmp(stim_name, 'testing') || strcmp(stim_name, 'behavior')
 
    %PSYCHOMETRIC PLOT
    
    %Combine all freqs into one x axis
    tempFreqs = [];
    for i = 1:numel(batch_freqs)
        tempFreqs = [tempFreqs, batch_freqs{i}];
    end
    x = unique(tempFreqs);
    
    %Obtain hit rate and reaction time per combined frequencies 
    y_mat = nan(max(nDatesPerMouse),length(x),nMice);
    rxn_mat = nan(max(nDatesPerMouse),length(x),nMice); 
    for a = 1:nMice
        for b = 1:nDatesPerMouse(a)
            for i = 1:length(x)
                current_x = x(i);
                freq_ind = find(batch_freqs{a,b} == current_x);
                if isempty(freq_ind)
                    continue
                end
                y_mat(b,i,a) = batch_hitRatePerFreq{a,b}(freq_ind);
                rxn_mat(b,i,a) = batch_rxnTimesPerFreq{a,b}(freq_ind);
            end
        end
    end

    %Estimate and plot psychometric functions
    batch_threshold = nan(nMice, max(nDatesPerMouse));
    batch_curves = nan(size(y_mat));
    binCurves = 0;

    figure
    c = parula(18); %Colormap
    
    for a = 1:nMice
        dateLegend = {};
        for b = 1:nDatesPerMouse(a)
            
            % Plot psychometric curve with bins of size 2
            y = y_mat(b,:,a);
            xx = x;
            xx(isnan(y)) = [];
            y(isnan(y)) = [];
            
            if isempty(y)
                continue
            end
            
            if binCurves
                binned_y = nanmean([y(1:2:end-1);y(2:2:end)]);
                binned_x = mean([xx(1:2:end-1);xx(2:2:end)]);
            else
                binned_y = y;
                binned_x = xx;
            end

            % Raw psychometric curves
            % Plot these on log scale
            subplot(4,nMice,a); hold all
            plot(1:length(binned_x),binned_y, 'Color', c(b,:), 'Linewidth', 2)
            ylabel('Hit rate')
            set(gca, 'XTick', 1:length(binned_x))
            set(gca, 'XTickLabel', x(1:floor(length(x)/length(binned_x)):length(x)))
            title(mice{a})
            xlim([1 length(binned_x)])
            dateLegend{b} = ['Day ' num2str(b)]; 
            
            % Fit psychometric functions
            targets = [0.25 0.5 0.75]; % 25 50 75 performance
            weights = ones(1,length(binned_y)); % No weighting
            [~, fit_curve, fit_threshold] = fitPsycheCurveLogit(1:length(binned_x), binned_y, weights, targets); %FOR GRAPH ONLY
            [~, ~, actual_threshold] = fitPsycheCurveLogit(binned_x, binned_y, weights, targets);

            if actual_threshold(2) > binned_x(end)
                batch_threshold(a,b) =  binned_x(end);
                warning(['Threshold ' num2str(actual_threshold(2)) ' out of bounds'])
            elseif actual_threshold(2) < binned_x(1)
                batch_threshold(a,b) = binned_x(1);
                warning(['Threshold ' num2str(actual_threshold(2)) ' out of bounds'])
            else
                batch_threshold(a,b) = actual_threshold(2);
            end
            
            subplot(4,nMice,nMice + a); hold all
            plot(fit_curve(:,1),fit_curve(:,2), 'Color', c(b,:), 'Linewidth', 2)
            %scatter(fit_threshold(2), targets(2), 'x', 'k')
            ylabel('Hit rate')
            xlabel('deltaF (octaves)')
            xlim([1 length(binned_x)])
            set(gca, 'XTick', 1:length(binned_x))
            set(gca, 'XTickLabel', x(1:floor(length(x)/length(binned_x)):length(x)))
        end
        %legend(dateLegend, 'Location', 'northwest')
    end

    % PLOT THRESHOLDS
    batch_threshold(batch_threshold == 0) = nan;
    subplot(4,nMice,[2*nMice + 1: 3*nMice]); hold all
    c = bone(5);
    for a = 1:nMice
        plot(batch_threshold(a,:), 'Color', c(a,:))       
    end
    legend(mice)
    plot(nanmean(batch_threshold),'LineWidth', 2, 'Color', 'k')
    ylim([0 0.5])
    xlim([1 10])
    ylabel('Threshold (octaves)')
    title('Discrimination threshold')

    % PLOT HIT RATE AND FP RATE
    subplot(4,nMice,[3*nMice + 1: 4*nMice]); hold all
    c1 = autumn(8);
    c2 = cool(8);
    for a = 1:nMice
        plot(batch_hitRate(a,:), 'Color', c1(a,:))
        plot(batch_FPrate(a,:), 'Color', c2(a,:))  
    end
    plot(nanmean(batch_hitRate),'LineWidth', 2, 'Color', 'k')
    plot(nanmean(batch_FPrate),'LineWidth', 2, 'Color', 'k')
    ylim([0 1])
    xlim([1 10])
    ylabel('Percent')
    xlabel('Training day')
    title('Hit rate and FP rate')
end
    
%% Save data

if loadPreviousData
    cd(PathName) %Save in the same place you loaded data from
    %save([FileName(1:end-4) '_reload'])
else
    if saveData
        cd(save_path)
        d = datestr(now,'yyyymmdd-HHMMSS');
        save(['Data_' d '.mat'], 'data', '-v7.3');
    end
end