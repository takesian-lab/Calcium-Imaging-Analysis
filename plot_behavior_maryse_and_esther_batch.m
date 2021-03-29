%plot_behavior_maryse_and_esther_batch

%% Load files

stim_protocol = 13; %Maryse behavior
stim_name = 'training';

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

        case 'TAKESIANLAB2P' %2P computer
            info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
            save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\Compiled Blocks';
            info_filename = 'Info';    
        otherwise
            disp('Computer does not match known users')
            return
    end

    cd(info_path)
    Info = importfile(info_filename);

    [data] = fillSetupFromInfoTable_v3(Info, compiled_blocks_path, stim_protocol, stim_name);
    disp('Data structure created.')
end

%% 
estimated_threshold = [];
estimated_threshold_bin2 = [];
raw_x = {};
raw_y = {};
freqs = {};
raw_x_bin2 = {};
raw_y_bin2 = {};
freqs_bin2 = {};
curve = {};
curve_bin2 = {};

for a=1:size(data.setup.mousename,1) %Mice
    for b=1:size(data.setup.mousename,2) %ROIs
        
       if isempty(data.setup.mousename{a,b})
           continue;
       end
       mouseID = data.setup.mousename{a,b};
       unique_block_name = data.setup.unique_block_names{a,b}(1);
       block = data.([mouseID]).([unique_block_name]);

       %% Get information from block

        %Block data
        mousename = char(block.setup.mousename);
        expt_date = char(block.setup.expt_date);
        dB_level = block.stim_level;
        plotTitle = strjoin({mousename, expt_date, num2str(dB_level), 'dB'});

        %Behavioral data
        outcomes = block.Outcome; %get outcomes (hit, miss, FP, withholds) from block
        earlyLicks = isnan(outcomes);
        outcomes_removeEL = outcomes(~earlyLicks);
        misses = outcomes_removeEL == 0;
        hits = outcomes_removeEL == 1;
        withholds = outcomes_removeEL == 3;
        FPs = outcomes_removeEL == 4;

        %Stimulus data
        allFrequencies = log2(block.TargetFreq(~earlyLicks)); %in log scale
        repeatingFrequency = mode(allFrequencies);
        allFrequencies = round(abs(allFrequencies - repeatingFrequency),2); %in octave difference
        table = tabulate(allFrequencies); %use the function tabulate to extract how many times each frequency was played
        uniqueFrequencies = table(:,1)';
        nRepsPerFrequency = table(:,2)'; %number of repetitions per frequency 
        holdingPeriod = block.holdingPeriod(~earlyLicks);

        %Reaction times
        raw_reactionTimes = block.rxn_time(~earlyLicks);
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

        % Fit psychometric functions
        targets = [0.25 0.5 0.75]; % 25 50 75 performance
        weights = ones(1,length(y)); % No weighting
        [~, fit_curve, fit_threshold] = fitPsycheCurveLogit(x, y, weights, targets); %FOR GRAPH ONLY
        [~, ~, actual_threshold] = fitPsycheCurveLogit(uniqueFrequencies, y, weights, targets);

        estimated_threshold(a,b) = actual_threshold(2);
        raw_x{a,b} = x;
        raw_y{a,b} = y;
        freqs{a,b} = uniqueFrequencies;
        curve{a,b} = fit_curve;
        
        % Plot psychometric curve with bins of size 2

        %set up axes
        A = hitRatePerFrequency;
        y = mean([A(1:2:end-1);A(2:2:end)]);
        x = 1:length(y);
  
        % Fit psychometric functions
        targets = [0.25 0.5 0.75]; % 25 50 75 performance
        weights = ones(1,length(y)); % No weighting
        [fit_coeffs, fit_curve, fit_threshold] = fitPsycheCurveLogit(x, y, weights, targets);
        [~, ~, actual_threshold] = fitPsycheCurveLogit(uniqueFrequencies(2:2:end), y, weights, targets);
        
        estimated_threshold_bin2(a,b) = actual_threshold(2);
        raw_x_bin2{a,b} = x;
        raw_y_bin2{a,b} = y;
        freqs_bin2{a,b} = uniqueFrequencies(2:2:end);
        curve_bin2{a,b} = fit_curve;
        
    end
end

%% Plot

estimated_threshold(estimated_threshold == 0) = nan;
estimated_threshold_bin2(estimated_threshold_bin2 == 0) = nan;

figure; hold all %estimated threshold
mice = {};

for a=1:4%size(data.setup.mousename,1) %Mice
       mouseID = data.setup.mousename{a,1};
       mice{a} = mouseID;
       
%        figure; hold all %estimated threshold
%        plot([estimated_threshold{a,:}]);
%        plot([estimated_threshold_bin2{a,:}]);
%        legend({'Threshold Bin 1', 'Threshold Bin 2'})
%        title(mouseID)
         
        plot(estimated_threshold_bin2(a,:));        
end
legend(mice(1:4))
plot(nanmean(estimated_threshold_bin2(1:4,:)),'LineWidth', 2, 'Color', 'k')
ylim([0 0.5])
xlim([1 15])

%%

figure; hold all %estimated threshold
mice = {};

for a=5:6%size(data.setup.mousename,1) %Mice
       mouseID = data.setup.mousename{a,1};
       mice{a} = mouseID;
       
%        figure; hold all %estimated threshold
%        plot([estimated_threshold{a,:}]);
%        plot([estimated_threshold_bin2{a,:}]);
%        legend({'Threshold Bin 1', 'Threshold Bin 2'})
%        title(mouseID)
         
        plot(estimated_threshold_bin2(a,:));        
end
legend(mice(5:6))
plot(nanmean(estimated_threshold_bin2(5:6,:)),'LineWidth', 2, 'Color', 'k')
ylim([0 0.5])
xlim([1 15])

%%

for a=1:size(data.setup.mousename,1) %Mice
       mouseID = data.setup.mousename{a,1};
       mice{a} = mouseID;
       
       figure; hold all %curve
       nCurves = 0;
       for c = 1:size(curve_bin2,2)
           if isempty(curve_bin2{a,c})
               continue
           end
           
           plot(curve_bin2{a,c}(:,1),curve_bin2{a,c}(:,2))
           nCurves = nCurves + 1;
       end
       legend(num2str(1:nCurves)')
       title(mouseID)
end

%% Save data

if loadPreviousData
    cd(PathName) %Save in the same place you loaded data from
    save([FileName(1:end-4) '_reload'])
else
    cd(save_path)
    d = datestr(now,'yyyymmdd-HHMMSS');
    save(['Data_' d '.mat'], 'data', '-v7.3');
end