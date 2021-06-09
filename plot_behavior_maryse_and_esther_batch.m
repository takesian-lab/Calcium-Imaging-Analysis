%plot_behavior_maryse_and_esther_batch

%% Load files

stim_protocol = 13; %Maryse behavior
stim_name = 'testing'; %training or testing or behavior
dB_to_include = 1; %0 = both, 1 = 60 only, 2 = 40 only, 3 = 60 unless no 60 for that day, 4 = 40 unless no 40 for that day
binCurves = 0; %Bin by 2 = 1
removeBadDays = 0; %don't include sessions where mouse performed badly on easy frequencies
easyFreqThreshold = 0.2; %easy frequency octaves
easyPerfThreshold = 0.8; %minimum hit rate

%%
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
            info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior';
            save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior';
            compiled_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior\Compiled Blocks v2';
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
    
    batch.hitRate = nan(nMice, max(nDatesPerMouse));
    batch.FPrate = nan(nMice, max(nDatesPerMouse));    
if strcmp(stim_name, 'testing')    ||  strcmp(stim_name, 'behavior') 
    batch.freqs = cell(nMice, max(nDatesPerMouse));
    batch.hitRatePerFreq = cell(nMice, max(nDatesPerMouse)); 
    batch.rxnTimesPerFreq = cell(nMice, max(nDatesPerMouse)); 
end

%% Extract data from blocks

for a = 1:nMice %Mice
    for b = 1:nDatesPerMouse(a) %Training days
        
       mouseID = mice{a};
       currentDate = dates{a}(b);
       dateInd = strcmp([data.setup.expt_date{a,:}], currentDate);
       matchingBlocks = data.setup.unique_block_names(a,dateInd);
       
       %Prepare variables to combine blocks together
       combinedLevels = [];
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
       
       switch dB_to_include
           case 1
               matchingBlocks(stimLevels == 40) = [];
               
           case 2
               matchingBlocks(stimLevels == 60) = [];
               
           case 3
               if unique(stimLevels) ~= 60
                   matchingBlocks(stimLevels == 40) = [];
               end
               
           case 4
               if unique(stimLevels) ~= 40
                   matchingBlocks(stimLevels == 40) = [];
               end
               
           otherwise
               %do nothing
       end

           
       for c = 1:length(matchingBlocks)
           unique_block_name = matchingBlocks{c};
           block = data.([mouseID]).([unique_block_name]);
           disp(unique_block_name)
           disp(block.stim_level)

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
            allFrequencies = round(abs(allFrequencies - repeatingFrequency),3); %in octave difference
            holdingPeriod = block.holdingPeriod(~earlyLicks);
            try
                holdingPeriod = holdingPeriod + block.waitPeriod(~earlyLicks);
            catch
            end

            %I removed the holding period in later versions of the task
            noHoldingPeriod = false;
            if length(unique(holdingPeriod)) == 1
                if unique(holdingPeriod) <= 0.2
                    noHoldingPeriod = true;
                    holdingPeriod(:) = 0;
                end
            end

            %If mouse is not responsive to easy frequencies, skip block
            if removeBadDays
                easyFrequencies = allFrequencies > easyFreqThreshold;
                easyPerformance = hits(easyFrequencies);
                if sum(easyPerformance)/length(easyPerformance) < easyPerfThreshold
                    continue
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
            combinedLevels   = [combinedLevels, block.stim_level];
            combinedFreqs    = [combinedFreqs, allFrequencies];
            combinedHits     = [combinedHits, hits];
            combinedMisses   = [combinedMisses, misses];
            combinedFPs      = [combinedFPs, FPs];
            combinedWitholds = [combinedWitholds, witholds];
            combinedRxnTimes = [combinedRxnTimes, reaction_times]; 
            
       end
       
        %Training data
        batch.levels{a,b} = combinedLevels;
        batch.hitRate(a,b) = sum(combinedHits)/(sum(combinedHits) + sum(combinedMisses));
        batch.FPrate(a,b) = sum(combinedFPs)/(sum(combinedFPs) + sum(combinedWitholds));
        
        %Pscyhometric data
        if strcmp(stim_name, 'testing') || strcmp(stim_name, 'behavior')
            if isempty(combinedFreqs) %No data
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

            batch.freqs{a,b} = uniqueFrequencies;
            batch.hitRatePerFreq{a,b} = hitRatePerFrequency; 
            batch.rxnTimesPerFreq{a,b} = rxnTimesPerFrequency;
        end
    end
end

%%  PLOTS

if strcmp(stim_name, 'training')
    
    %HIT RATE AND FP RATE
    figure;
    
    subplot(1,2,1)
    plot(batch.hitRate'); hold on
    ylim([0 1])
    ylabel('Percent')
    xlabel('Training day')
    title('Hit rate')
    legend(mice, 'Location', 'southwest')
    
    subplot(1,2,2)
    plot(batch.FPrate')
    ylim([0 1])
    xlabel('Training day')
    title('FP rate')
    legend(mice, 'Location', 'southwest')
    
elseif strcmp(stim_name, 'testing') || strcmp(stim_name, 'behavior')
 
    %PSYCHOMETRIC PLOT
    
    %Combine all freqs into one x axis
    tempFreqs = [];
    for i = 1:numel(batch.freqs)
        tempFreqs = [tempFreqs, batch.freqs{i}];
    end
    x = unique(tempFreqs);
    
    %Obtain hit rate and reaction time per combined frequencies 
    y_mat = nan(max(nDatesPerMouse),length(x),nMice);
    rxn_mat = nan(max(nDatesPerMouse),length(x),nMice); 
    for a = 1:nMice
        for b = 1:nDatesPerMouse(a)
            for i = 1:length(x)
                current_x = x(i);
                freq_ind = find(batch.freqs{a,b} == current_x);
                if isempty(freq_ind)
                    continue
                end
                y_mat(b,i,a) = batch.hitRatePerFreq{a,b}(freq_ind);
                rxn_mat(b,i,a) = batch.rxnTimesPerFreq{a,b}(freq_ind);
            end
        end
    end
    
    %Estimate and plot psychometric functions
    %FIRST all plots per mice to show goodness of fit [figType = 1]
    %SECOND all mice in same graph [figType = 2]
    batch.threshold = nan(nMice, max(nDatesPerMouse));
    batch.curves = nan(size(y_mat));
    
    c = parula(18); %Colormap

    for figType = 1:2
        if figType == 2; figure; end
    
    for a = 1:nMice
        if figType == 1; figure; end
        dateLegend = {};
        for b = 1:nDatesPerMouse(a)
            nPlotsX = ceil(nDatesPerMouse(a)/5);
            nPlotsY = ceil(nDatesPerMouse(a)/nPlotsX);
            
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
            if figType == 1
                subplot(nPlotsX, nPlotsY, b); hold all
                plot(1:length(binned_x),binned_y, 'Linewidth', 2)
                title(batch.levels{a,b})
            elseif figType == 2
                subplot(4,nMice,a); hold all
                plot(1:length(binned_x),binned_y, 'Color', c(b,:), 'Linewidth', 2)
                title(mice{a})
            end
            ylabel('Hit rate')
            set(gca, 'XTick', 1:length(binned_x))
            set(gca, 'XTickLabel', x(1:floor(length(x)/length(binned_x)):length(x)))
            xlim([1 length(binned_x)])
            ylim([0 1])
            dateLegend{b} = ['Day ' num2str(b)]; 
            
            % Fit psychometric functions
            targets = [0.25 0.5 0.75]; % 25 50 75 performance
            u = mean(hitRatePerFrequency);
            v = std(hitRatePerFrequency);
            SPs = [0.1, 0.1,  inf, inf; % Upper limits for g, l, u ,v
                   0.01, 0.05, u,  v;  % Start points for g, l, u ,v
                   0,    0,    0,  0]; % Lower limits for g, l, u ,v
               
            [~, fit_curve] = FitPsycheCurveWH(1:length(binned_x), binned_y, SPs); %FOR GRAPH ONLY
            fit_threshold = nan(1,length(targets));
            for f = 1:length(targets)
                [~, min_ind] = min(abs(fit_curve(:,2) - targets(f)));
                fit_threshold(f) = fit_curve(min_ind,1);
            end
            
            %Actual values
            [coeffs, fit_curve2] = FitPsycheCurveWH(binned_x, binned_y, SPs); %FOR GRAPH ONLY
            actual_threshold = nan(1,length(targets));
            for f = 1:length(targets)
                [~, min_ind] = min(abs(fit_curve2(:,2) - targets(f)));
                actual_threshold(f) = fit_curve2(min_ind,1);
            end
            %weights = ones(1,length(binned_y)); % No weighting
            %[~, fit_curve, fit_threshold] = fitPsycheCurveLogit(1:length(binned_x), binned_y, weights, targets); %FOR GRAPH ONLY
            %[~, ~, actual_threshold] = fitPsycheCurveLogit(binned_x, binned_y, weights, targets);

            if actual_threshold(2) > binned_x(end)
                batch.threshold(a,b) =  binned_x(end);
                warning(['Threshold ' num2str(actual_threshold(2)) ' out of bounds'])
            elseif actual_threshold(2) < binned_x(1)
                batch.threshold(a,b) = binned_x(1);
                warning(['Threshold ' num2str(actual_threshold(2)) ' out of bounds'])
            else
                batch.threshold(a,b) = actual_threshold(2);
            end
            
            if figType == 1
                subplot(nPlotsX, nPlotsY, b); hold all
                plot(fit_curve(:,1),fit_curve(:,2), 'Color', 'g', 'Linewidth', 2)
                hline(0.5, 'r')
                scatter(fit_threshold, targets, 'x', 'k')
            elseif figType == 2
                subplot(4,nMice,nMice + a); hold all
                plot(fit_curve(:,1),fit_curve(:,2), 'Color', c(b,:), 'Linewidth', 2)
            end
            ylabel('Hit rate')
            xlabel('deltaF (octaves)')
            xlim([1 length(binned_x)])
            ylim([0 1])
            set(gca, 'XTick', 1:length(binned_x))
            set(gca, 'XTickLabel', x(1:floor(length(x)/length(binned_x)):length(x)))
        end
        %legend(dateLegend, 'Location', 'northwest')
        if figType == 1
            suptitle(mice{a})
        end
    end

    end

    % PLOT THRESHOLDS
    batch.threshold(batch.threshold == 0) = nan;
    subplot(4,nMice,[2*nMice + 1: 3*nMice]); hold all
    c = bone(nMice);
    for a = 1:nMice
        plot(batch.threshold(a,:), 'Color', c(a,:))       
    end
    legend(mice)
    plot(nanmean(batch.threshold),'LineWidth', 2, 'Color', 'k')
    ylim([0 0.5])
    xlim([1 10])
    ylabel('Threshold (octaves)')
    title('Discrimination threshold')

    % PLOT HIT RATE AND FP RATE
    subplot(4,nMice,[3*nMice + 1: 4*nMice]); hold all
    c1 = autumn(8);
    c2 = cool(8);
    for a = 1:nMice
        plot(batch.hitRate(a,:), 'Color', c1(a,:))
        plot(batch.FPrate(a,:), 'Color', c2(a,:))  
    end
    plot(nanmean(batch.hitRate),'LineWidth', 2, 'Color', 'k')
    plot(nanmean(batch.FPrate),'LineWidth', 2, 'Color', 'k')
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