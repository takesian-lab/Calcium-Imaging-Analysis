% Extract cell-specific data from blocks and save in excel format
% v3 works with Info spreadsheets

%stim protocol code is:
%noiseburst = 1
%ReceptiveField = 2
%FM sweep = 3
%SAM = 5
%widefield = 4
%SAM freq = 6
%Behavior = 7 and 8
%Random H20 = 9
%Noiseburst_ITI = 10
%Random air puff = 11

%% Set up environment

clear all
close all

PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
        info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis';
        blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\CompiledBlocks_SFN';
        save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\ExtractedData';
        info_filename = 'Combined Info SFN';
       
        stimProtocol = 2; %See stimProtocol list above
        use_zscore = 1; %use z-score instead of df/f
        STDlevel = 2; %minimum # standard deviations above baseline (or Z-scores) to be considered significant
        AUC_F_level = 5; %minimum area under the curve to be considered significant for df/f traces
        AUC_S_level = 10; %minimum area under the curve to be considered significant for spike traces
        sort_active = 1;  %0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 0; %plot stim responses for each cell
        plot_tuning = 0; %plot tuning for each cell (RF only)
        save_data = 1;
        save_figures = 0;
        
    case 'RD0332' %Carolyn
        info_path = 'Z:\Carolyn\2P Imaging data\SSRI response stimuli pilot\VxDD062420M3';
        blocks_path = 'Z:\Carolyn\2P Imaging data\SSRI response stimuli pilot\compiled blocks';
        save_path = 'Z:\Carolyn\2P Imaging data\SSRI response stimuli pilot\VxDD062420M3\extracted';
        info_filename = 'VxDD062420M3_postFLX';

        stimProtocol = 2; %See stimProtocol list above
        use_zscore = 1; %use z-score instead of df/f
        STDlevel = 2; %minimum # standard deviations above baseline (or Z-scores) to be considered significant
        AUC_F_level = 5; %minimum area under the curve to be considered significant for df/f traces
        AUC_S_level = 10; %minimum area under the curve to be considered significant for spike traces
        sort_active = 1;  %0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 1; %plot stim responses for each cell
        plot_tuning = 1; %plot tuning for each cell (RF only)
        save_data = 1;
        save_figures = 1;
        
    case 'RD-6-TAK1' %Anne
        info_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
        blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
        save_path = 'Y:\Anne\VIPvsNDNF_Analysis_FINAL2';
        info_filename = 'ResponsiveCells';
 
        stimProtocol = 2; %See stimProtocol list above
        use_zscore = 1; %use z-score instead of df/f
        STDlevel = 2; %minimum # standard deviations above baseline (or Z-scores) to be considered significant
        AUC_F_level = 5; %minimum area under the curve to be considered significant for df/f traces
        AUC_S_level = 10; %minimum area under the curve to be considered significant for spike traces
        sort_active = 1;  %0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 1; %plot stim responses for each cell
        plot_tuning = 1; %plot tuning for each cell (RF only)
        save_data = 1;
        save_figures = 1;
        
    otherwise
        error('Computer does not match known users')
end

%% Load data and initialize variables for final spreadsheet

%COLUMN HEADERS FOR FINAL SPREADSHEET
nominalColumnHeaders = {'Group', 'Mouse ID', 'FOV', 'Data type', 'Block', 'Cell Number', 'F-Activity', 'S-Activity'};
numericalColumnHeaders = {'Peak Amplitude', 'P1', 'Peak Latency', 'Peak Width', 'Peak Area', 'Trough Amplitude', 'T1', 'Trough Latency', 'Trough Width', 'Trough Area',...
    'Combined Amplitude', 'Combined Latency', 'Combined Peak Latency', 'Combined Width', 'Combined Area'};

cd(info_path)
Info = importfile(info_filename);

%Create data structure for files corresponding to stim_protocol
disp('Loading blocks to include...')
[blockData] = fillSetupFromInfoTable_v3(Info, blocks_path, stimProtocol);

%Make cell list with nominal info to store and preallocate space for variables
cellList = {};
unique_block_names = {};
fsList = [];
constants = [];
all_v1 = [];
all_v2 = [];

for a = 1:size(blockData.setup.mousename,1) %Mice
        mouseID = blockData.setup.mousename{a};
    for f = 1:size(blockData.setup.mousename,2) %FOVs
        for b = 1:size(blockData.setup.unique_block_names{a,f},2) %Blocks
            unique_block_name = blockData.setup.unique_block_names{a,f}(b);
            block = blockData.([mouseID]).([unique_block_name]);
            nCells = length(block.cell_number);
            v1 = block.parameters.variable1;
            v2 = block.parameters.variable2;

            clear tempCellList;
            clear temp_unique;
            clear temp_fsList;
            clear temp_constants;

            [tempCellList(1:nCells,1)] = deal(block.setup.expt_group);
            [tempCellList(1:nCells,2)] = deal(mouseID);
            [tempCellList(1:nCells,3)] = deal(block.setup.FOV);
            [tempCellList(1:nCells,4)] = deal(blockData.setup.stim);
            [tempCellList(1:nCells,5)] = deal(block.setup.block_filename);
            [tempCellList(1:nCells,6)] = block.cell_number;
            [temp_unique(1:nCells,1)] = deal(unique_block_name);
            [temp_fsList(1:nCells,1)] = deal(block.setup.framerate);
            [temp_constants(1:nCells,1)] = deal(block.setup.constant.baseline_length);
            [temp_constants(1:nCells,2)] = deal(block.setup.constant.after_stim);

            cellList = [cellList; tempCellList];
            unique_block_names = [unique_block_names; temp_unique];
            fsList = [fsList; temp_fsList];
            constants = [constants; temp_constants];
            all_v1 = [all_v1, v1];
            all_v2 = [all_v2, v2];
            
        end
    end
end

%Determine what size the raster will be
%If blocks have different frame rates or trial durations (e.g. baseline)
%then we will have a problem because the traces won't match up.
if length(unique(constants(:,1))) ~= 1 || length(unique(constants(:,2))) ~= 1
    error('Blocks have different trial durations. Check and recompile blocks to match.')
end

if length(unique(fsList)) == 1 
    fs = unique(fsList);
    baseline = unique(constants(:,1));
    after = unique(constants(:,2));
    baseline_inFrames = round(baseline*fs);
    after_inFrames = round(after*fs);
    nFrames = baseline_inFrames + after_inFrames;
else
    error('Blocks have different framerates.')
end

%Units for figures
if use_zscore
    units = {'z-score', 'spikes'};
else
    units = {'df/f', 'spikes'};
end

figures = cell(size(cellList,1),21); %Store figures (if save_figures == 1)
activity = cell(size(cellList,1),2); %Auto-determined activity (activated/prolonged/suppressed)
[data_F, data_S] = deal(nan(size(cellList,1),length(numericalColumnHeaders))); %Numerical data
[raster_F, raster_S] = deal(nan(size(cellList,1),nFrames));

%% Prepare stim-specific variables

%RF will contain data unique to each stim protocol
RF = struct;
st = struct;

if stimProtocol == 2 %RF
    RF.stim = {'Intensity', 'Frequency'}; %corresponds to v1 and v2
    RF.ColumnHeaders = {'Sparseness', 'BF', 'BF_I', 'CF', 'CF_I', 'BestInt', 'Int_RMS', 'Int_halfwidth',...
        'BW_20_RMS', 'BW_20_halfwidth', 'BW_BF_RMS', 'BW_BF_halfwidth', 'ISI', 'dPrime', 'Type'};
    st.FRA = {};
elseif stimProtocol == 3 %FM
    RF.stim = {};
    RF.ColumnHeaders = {'TBD'};
elseif stimProtocol == 5 %SAM
    RF.stim = {};
    RF.ColumnHeaders = {'TBD'};
elseif stimProtocol == 6 %SAMfreq
    RF.stim = {};
    RF.ColumnHeaders = {'TBD'};
end

%Preallocate space to save "receptive field" data for all stim types
%Find all possible stim combinations (not including 0dB 'blank' trials)

switch stimProtocol
    case {2,3} %RF, FM
        all_v2(all_v2 == 0) = [];

    case {9,10,11} %water, NoiseITI, air
        all_v1(all_v1 == 0) = [];

    case {5,6} %SAM, SAMfreq %NAN trials instead of zeros           
        all_v1(isnan(all_v1)) = [];

    case {1,12} %Noiseburst, spontaneous
%         [stim_v1, stim_v2, stim_v0] = deal(ones(1,nTrials));

    otherwise
        error('Stim type is currently not compatible with removing 0dB trials')
end

unique_all_v1 = unique(all_v1)';
unique_all_v2 = unique(all_v2)';
n1 = length(unique_all_v1);
n2 = length(unique_all_v2);

%If RF, swap v1/v2 and order intensities from highest to lowest
if stimProtocol == 2
    [unique_all_v2, unique_all_v1] = deal(unique_all_v1, unique_all_v2);  
    unique_all_v1 = flipud(unique_all_v1);
end

[st.response, st.isResponsive, st.peak, st.peakavg, st.area, RF.stim_v1, RF.stim_v2] = deal(nan(n1,n2,size(cellList,1)));
st.data = nan(size(cellList,1), length(RF.ColumnHeaders));
st.sp_by_v1 = nan(size(cellList,1),n1);
st.sp_by_v2 = nan(size(cellList,1),n2);
RF.F = st;
RF.S = st;
    
%% Loop through all cells in each block

disp('Extracting data...')
oldBlock = '';

for c = 1:size(cellList,1)
    mouseID = cellList{c,2};
    unique_block_name = unique_block_names{c,1};
    block = blockData.([mouseID]).([unique_block_name]);
    if ~strcmp(oldBlock, cellList{c,5})
        oldBlock = cellList{c,5};
        disp(oldBlock)
    end
    cellNumber = str2double(cellList{c,6});
    cellIndex = find(block.cell_number == cellNumber);
    
    baseline_length = block.setup.constant.baseline_length; %seconds
    framerate = block.setup.framerate;
    nBaselineFrames = baseline_length*framerate; %frames
    
    stim_v1 = block.parameters.variable1';
    stim_v2 = block.parameters.variable2';
    
    %Identify loco trials to remove
    remove = [];
    if sort_active == 1
        remove = find(block.active_trials == 1); %active trials
    elseif sort_active == 2
        remove = find(block.active_trials == 0); %inactive trials
    end

    %Separate 0dB 'blank' sound trials
    %depending on stim type, 0dB trials are stored in variable1 or variable2
    switch stimProtocol
        case {2,3} %RF, FM
            stim_v0 = stim_v2;
            
        case {9,10,11} %water, NoiseITI, air
            stim_v0 = stim_v1;
            stim_v2 = zeros(size(stim_v1));
            
        case {5,6} %SAM, SAMfreq %NAN trials instead of zeros           
            stim_v0 = stim_v1;
            stim_v0(isnan(stim_v0)) = 0;

        case {1,12} %Noiseburst, spontaneous
            if length(stim_v1) > 1 || length(stim_v2) > 1
                error('Check stim parameters')
            end
            nTrials = length(block.Sound_Time);
            [stim_v1, stim_v2, stim_v0] = deal(ones(1,nTrials));

        otherwise
            error('Stim type is currently not compatible with removing 0dB trials')
    end
    
    %Separate blank and stim trials
    blankTrials = stim_v0 == 0; %0dB trials
    stimTrials = ~blankTrials;
    
    %Get list of all possible stim without blanks prior to removing loco trials
    unique_stim_v1 = unique(stim_v1(stim_v0 ~= 0));
    unique_stim_v2 = unique(stim_v2(stim_v0 ~= 0));
    
    %If RF, swap v1/v2 and order intensities from highest to lowest
    if stimProtocol == 2 %RF
        [stim_v2, stim_v1] = deal(stim_v1, stim_v2);
        [unique_stim_v2, unique_stim_v1] = deal(unique_stim_v1, unique_stim_v2);
        unique_stim_v1 = flipud(unique_stim_v1);
    end
    
    %Find matching indices of unique_all_v1 and unique_all_v2
    ind_v1 = nan(1,length(unique_stim_v1));
    ind_v2 = nan(1,length(unique_stim_v2));
    for i = 1:length(unique_stim_v1)
        ind_v1(i) = find(unique_all_v1 == unique_stim_v1(i));
    end
    for i = 1:length(unique_stim_v2)
        ind_v2(i) = find(unique_all_v2 == unique_stim_v2(i));
    end 
    
    %Remove loco trials
    blankTrials(remove,:) = 0;
    stimTrials(remove,:) = 0;
    
    %Get trial indices
    blankTrials = find(blankTrials);
    stimTrials = find(stimTrials);
    
    %Only keep stim values for non-blank trials
    stim_v1 = stim_v1(stimTrials);
    stim_v2 = stim_v2(stimTrials);
              
    %Pull out all stim-aligned traces for this cell
    F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:));
    F7_baseline = F7(:,1:nBaselineFrames); %baseline for each trial
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %compute df/f: (total-mean)/mean
    F7_zscore = (F7-nanmean(F7_baseline,2))./nanstd(F7_baseline,0,2); %compute zscore: (total-mean)/std
    spks = squeeze(block.aligned_stim.spks_stim(cellIndex,:,:));

    %Remove trials with infinite values [this was a bug in a small number of blocks]
    [inf_rows,~] = find(isinf(F7_df_f));
    remove_inf = unique(inf_rows);
    if ~isempty(remove_inf)
        stim_v1(stimTrials == remove_inf) = [];
        stim_v2(stimTrials == remove_inf) = [];
        stimTrials(stimTrials == remove_inf) = [];
        blankTrials(blankTrials == remove_inf) = [];
    end
    
    %Separate stim and blank trials
    if use_zscore
        F = F7_zscore(stimTrials,:);
        F_blanks = F7_zscore(blankTrials,:);
    else
        F = F7_df_f(stimTrials,:);
        F_blanks = F7_df_f(blankTrials,:);
    end
    S = spks(stimTrials,:);
    S_blanks = spks(blankTrials,:);

    %GET AVERAGED RESPONSES
    %check if each condition is active, then concatenate and keep only active conditions
    analyze_by_stim_condition = 1; %I assume we'll almost always want to do this
    if analyze_by_stim_condition

        F_by_Stim = [];
        S_by_Stim = [];

        for v = 1:length(unique_stim_v1)
            for vv = 1:length(unique_stim_v2)
                vi = ind_v1(v);
                vvi = ind_v2(vv);
                stim_rows = intersect(find(stim_v1 == unique_stim_v1(v)), find(stim_v2 == unique_stim_v2(vv)));
                
                if isempty(stim_rows) %Some stim combinations might not exist due to loco
                    continue
                end

                %Plot individual figures of each stim condition
                plotStimFigures = 0; %Use for troubleshooting
                if plotStimFigures; figure; end

                for i = 1:2 % F and S
                    if i == 1
                        A = 'F';
                        nFramesToAvg = framerate; %1s
                        trials = F(stim_rows,:);
                        blank_trials = F_blanks;
                        AUClevel = AUC_F_level;
                        zflag = use_zscore;
                    elseif i == 2
                        A = 'S';
                        nFramesToAvg = round(framerate/2); %0.5s
                        trials = S(stim_rows,:);
                        blank_trials = S_blanks;
                        AUClevel = AUC_S_level;
                        zflag = 0;
                    end

                    if plotStimFigures; subplot(1,2,i); end
                    
                    [isActive, tempActivity, pk_temp, tr_temp] = checkIfActive_v2(trials, nBaselineFrames, framerate, STDlevel, AUClevel, plotStimFigures, units{i}, zflag);
                    %[isActiveComparedToBlanks, pValue] = compare_to_blank_trials(trials, blank_trials);
                    if isActive
                        if i == 1
                            F_by_Stim = [F_by_Stim; trials];
                        elseif i == 2
                            S_by_Stim = [S_by_Stim; trials];
                        end
                    end

                    %RF data
                    y = smooth(nanmean(trials,1),3)';
                    response = nanmean(y(1,nBaselineFrames:nBaselineFrames+nFramesToAvg));

                    %Store either peak or trough data
                    if strcmp(tempActivity, 'suppressed')
                        combinedData = tr_temp;
                    else
                        combinedData = pk_temp;
                    end

                    %Compute average of peak or trough region for receptive field
                    if any(combinedData) %skip if nans
                        c1 = combinedData(2) + nBaselineFrames; %latency and correct for baseline
                        c2 = c1 + combinedData(4); %t1 + width
                        peakavg = nanmean(y(1,c1:c2));
                    else
                        peakavg = nan;
                    end

                    RF.stim_v1(vi,vvi,c) = unique_stim_v1(v);
                    RF.stim_v2(vi,vvi,c) = unique_stim_v2(vv);
                    RF.(A).response(vi,vvi,c) = response;
                    RF.(A).isResponsive(vi,vvi,c) = isActive;                            
                    RF.(A).peak(vi,vvi,c) = combinedData(1);
                    RF.(A).peakavg(vi,vvi,c) = peakavg;
                    RF.(A).area(vi,vvi,c) = combinedData(5);
                end
                   if plotStimFigures; suptitle(['Cell ' num2str(c) ': Stim ' num2str(unique_stim_v1(v)) '/' num2str(unique_stim_v2(vv))]); end
            end
        end

        %Unless none of the stim conditions end up being significant,
        %use significant conditions only for average
        if ~isempty(F_by_Stim) 
            F = F_by_Stim;
        end

        if ~isempty(S_by_Stim)
            S = S_by_Stim;
        end                   
    end

    %FILL RASTER
    avg_F = nanmean(F,1)';
    avg_S = nanmean(S,1)';

    %Data might be empty if too many trials were removed due to loco activity 
    if isempty(avg_F)
        avg_F = nan(1,nFrames);
    end
    if isempty(avg_S)
        avg_S = nan(1,nFrames);
    end
    raster_F(c,1:length(avg_F)) = avg_F;
    raster_S(c,1:length(avg_S)) = avg_S;


    %COMPUTE LATENCIES, AMPLITUDES, AND WIDTHS
    if plot_graphs; figures{c,1} = figure; hold on; end

    for i = 1:2 %F and S
        if i == 1
            trials = F;
            blankTrials = F_blanks;
            AUClevel = AUC_F_level;
            zflag = use_zscore;
        elseif i == 2
            trials = S;
            blankTrials = S_blanks;
            AUClevel = AUC_S_level;
            zflag = 0;
        end

        if plot_graphs; subplot(1,2,i);end

        [~, tempActivity, peak_data, trough_data] = checkIfActive_v2(trials, nBaselineFrames, framerate, STDlevel, AUClevel, plot_graphs, units{i}, zflag);

        %Store activity
        activity{c,i} = tempActivity;

        %Convert to seconds and store data
        peak_data(1:4) = peak_data(1:4)./framerate;
        trough_data(1:4) = trough_data(1:4)./framerate;

        %For storing combined data, only keep trough OR peak data in one column
        combined_data = nan(1,5);
        if isequal(tempActivity, 'activated') || isequal(tempActivity, 'prolonged')
            combined_data = peak_data;
        elseif isequal(tempActivity, 'suppressed')
            combined_data = trough_data;
        end

        if i == 1 %GCaMP
            data_F(c,:) = [peak_data, trough_data, combined_data];
        elseif i == 2 %Spikes
            data_S(c,:) = [peak_data, trough_data, combined_data];
        end

    end

    if plot_graphs; suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber))); end 

    %RECEPTIVE FIELD ONLY
    if stimProtocol == 2 %RF
        
        %Compute tuning sparseness
        [RF.F.data(c,1), RF.F.sp_by_v2(c,ind_v1), RF.F.sp_by_v1(c,ind_v2), figures{c,2}] = compute_tuning_sparseness_v2(RF.F.peak(ind_v1,ind_v2,c), plot_tuning, 'df/f',  unique_stim_v2, unique_stim_v1, cellNumber);
        [RF.S.data(c,1), RF.S.sp_by_v2(c,ind_v1), RF.S.sp_by_v1(c,ind_v2), figures{c,3}] = compute_tuning_sparseness_v2(RF.S.peak(ind_v1,ind_v2,c), plot_tuning, 'spikes', unique_stim_v2, unique_stim_v1, cellNumber);  

        %Compute frequency tuning
        if any(any(RF.F.isResponsive(ind_v1,ind_v2,c)))
            [RF.F.FRA{c,1}, RF.F.data(c,2:end), figures{c,4:8}] = ... 
            compute_frequency_tuning_v2(RF.F.isResponsive(ind_v1,ind_v2,c), RF.F.response(ind_v1,ind_v2,c),...
            plot_tuning, unique_stim_v2, unique_stim_v1, cellNumber);
        end
        
        if any(any(RF.S.isResponsive(ind_v1,ind_v2,c)))
            [RF.S.FRA{c,1}, RF.S.data(c,2:end), figures{c,9:13}] = ... 
            compute_frequency_tuning_v2(RF.S.isResponsive(ind_v1,ind_v2,c), RF.S.response(ind_v1,ind_v2,c),...
            plot_tuning, unique_stim_v2, unique_stim_v1, cellNumber);
        end

        if plot_tuning
            %Plot RF using various measures
            figures{c,14} = figure('units','normalized','outerposition',[0 0 0.5 1]); hold on

            subplot(5,2,1)
            imagesc(RF.F.response(ind_v1,ind_v2,c))
            title('df/f: mean response')
            ylabel('dB')

            subplot(5,2,2)
            imagesc(RF.S.response(ind_v1,ind_v2,c))
            title('spikes: mean response')

            subplot(5,2,3)
            imagesc(RF.F.peak(ind_v1,ind_v2,c))
            title('df/f: peak amplitude');
            ylabel('dB')

            subplot(5,2,4)
            imagesc(RF.S.peak(ind_v1,ind_v2,c))
            title('spikes: peak amplitude')

            subplot(5,2,5)
            imagesc(RF.F.peakavg(ind_v1,ind_v2,c))
            title('df/f: average around peak');
            ylabel('dB')

            subplot(5,2,6)
            imagesc(RF.S.peakavg(ind_v1,ind_v2,c))
            title('spikes: average around peak')

            subplot(5,2,7)
            imagesc(RF.F.area(ind_v1,ind_v2,c))
            title('df/f: response area')
            ylabel('dB')

            subplot(5,2,8)
            imagesc(RF.S.area(ind_v1,ind_v2,c))
            title('spikes: response area')

            subplot(5,2,9)
            imagesc(RF.F.isResponsive(ind_v1,ind_v2,c))
            title('df/f: is responsive (0 or 1)')
            xlabel('Frequencies')
            ylabel('dB')

            subplot(5,2,10)
            imagesc(RF.S.isResponsive(ind_v1,ind_v2,c))
            title('spikes: is responsive (0 or 1)')
            xlabel('Frequencies')

            suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber)))
            set(findobj(gcf,'type','axes'),'XTickLabel',unique_stim_v2,'YTickLabel', unique_stim_v1,...
                   'XTick',1:length(unique_stim_v2), 'YTick',1:length(unique_stim_v1));

            %visualize_cell plots
            [figures{c,15:21}] = visualize_cell_AT_v2(block,cellNumber,[1,2,4]);
        end
    end

    if save_figures
        mouse = num2str(block.setup.mousename);
        block_filename = num2str(block.setup.block_filename);
        cellPath = num2str(cellNumber);

        cd(save_path)
        [~, ~, ~] = mkdir(save_path,mouse);
        [~, ~, ~] = mkdir(mouse,block_filename);
        cd([mouse '\' block_filename]);
        [~, ~, ~] = mkdir(cellPath);
        cd(cellPath); 

        figureNames = {'detect_activity'...
                      'sparseness_F'...
                      'sparseness_S'...
                      'intensityfit_F'...
                      'gaussfitCF_F'...
                      'gaussfitBW20_F'...
                      'gaussfitBF_F'...
                      'FRA_F'...
                      'intensityfit_S'...
                      'gaussfitCF_S'...
                      'gaussfitBW20_S'...
                      'gaussfitBF_S'...
                      'FRA_S'...
                      'RFs'...
                      'trace'...
                      'trace_avg'...
                      'trace_by_freq'...
                      'trace_by_int'...
                      'trace_RF_F'...
                      'trace_RF_S'...
                      'RF_simple'};           

        for f = 1:size(figures,2)           

            if isempty(figures{c,f}); continue; end

            h_name = [num2str(f) '_' mouse '_' cellPath '_' figureNames{f}];
            saveas(figures{c,f}, h_name, 'fig')
            saveas(figures{c,f}, h_name, 'jpg')               
        end

        close all;
    end
end
    
ExtractedData = struct;
ExtractedData.StimProtocl = stimProtocol;
ExtractedData.StimType = blockData.setup.stim;
ExtractedData.STDlevel = STDlevel;
ExtractedData.AUC_F_level = AUC_F_level;
ExtractedData.AUC_S_level = AUC_S_level;
ExtractedData.Sort_Active = sort_active;
ExtractedData.Use_ZScore = use_zscore;
%ExtractedData.BlockData = blockData; %this will make files big, but we could include it
ExtractedData.ColumnHeaders = [nominalColumnHeaders, strcat('F-',numericalColumnHeaders), strcat('S-',numericalColumnHeaders)];
ExtractedData.Activity = activity;
ExtractedData.NominalData = cellList;
ExtractedData.NumericalData = [data_F,data_S];
ExtractedData.Calcium_Raster = raster_F;
ExtractedData.Spikes_Raster = raster_S;
if stimProtocol == 2 %RF
    ExtractedData.RF = RF;
end
if save_figures
    ExtractedData.FigureNames = figureNames;
    ExtractedData.Figures = figures;
end

%% Save extracted data
if save_data == 1
    cd(save_path)
    d = datestr(now,'yyyymmdd-HHMMSS');
    filename = ['extractedData_' blockData.setup.stim '_' d '.mat'];
    ExtractedData.Filename = filename;
    save(filename, 'ExtractedData');
end

%% Plot sorted rasters (This section should be able to run on its own after just loading ExtractedData)
    
%Plot by activity type and peak/trough amplitude
%Sorting based on df/f data instead of spikes

sortType = 1; %sort by F activity (1) or S activity (2)
activityList = [ExtractedData.Activity(:,sortType)];
groupList = ExtractedData.NominalData(:,1);
Groups = unique(groupList);
Activity = {'activated', 'prolonged', 'suppressed'};
groupRasters = struct;
groupAverages = struct;
groupActivity = {};

for g = 1:length(Groups)
    currentCells = strcmpi(groupList,Groups(g));
        
    resorted_raster_F = [];
    resorted_raster_S = [];
    average_F = [];
    average_S = [];
    activity_type = [];
        
    for i = 1:length(Activity)       
        activeRows = strcmpi(activityList, Activity{i});
        currentRows = and(currentCells, activeRows);
        cellNumbers = find(currentRows); %cell numbers to use for comparison in locomotor 

        %find rows in current activity type and sort by amplitude and/or latency
        current_raster_F = ExtractedData.Calcium_Raster(currentRows,:);
        current_raster_S = ExtractedData.Spikes_Raster(currentRows,:);

        %Store smoothed average for plots
        average_F(i,:) = smooth(nanmean(current_raster_F,1));
        average_S(i,:) = smooth(nanmean(current_raster_S,1));
           
        %Use combined peak data to sort by peaks for activated/prolonged and troughs for suppressed
        current_F_amplitude = ExtractedData.NumericalData(currentRows,10); %MAGIC NUMBER
        current_S_amplitude = ExtractedData.NumericalData(currentRows,26); %MAGIC NUMBER
         
        [~, F_sort_ind] = sort(current_F_amplitude, 'descend');
        [~, S_sort_ind] = sort(current_S_amplitude, 'descend');

        %Sort based on sortType
        if sortType == 1 %Based on F
            sort_ind = F_sort_ind;
        elseif sortType == 2 %Based on S
            sort_ind = S_sort_ind;
        end

        resorted_raster_F = [resorted_raster_F; current_raster_F(sort_ind,:)];
        resorted_raster_S = [resorted_raster_S; current_raster_S(sort_ind,:)];

        temp_activity = cell(size(current_raster_F,1),1);
        [temp_activity{:}] = deal(Activity(i));
        activity_type = [activity_type; temp_activity];
    
        %Save cell order (makeValidName gets rid of invalid characters like '-' or '/')
        ExtractedData.CellsInOrder.(matlab.lang.makeValidName([Groups{g}])).([Activity{i}]).F = cellNumbers(sort_ind); 
        ExtractedData.CellsInOrder.(matlab.lang.makeValidName([Groups{g}])).([Activity{i}]).S = cellNumbers(sort_ind);
    end
                
    groupRasters.F{g} = resorted_raster_F;
    groupRasters.S{g} = resorted_raster_S;
    groupAverages.F{g} = average_F;
    groupAverages.S{g} = average_S;
    groupActivity{g} = activity_type;
end
        

%% PLOT FIGURE BY GROUP

figure; hold on
suptitle(ExtractedData.StimType)

for g = 1:length(Groups)

    subplot(6,length(Groups),g); hold all
    plot(groupAverages.F{g}', 'LineWidth', 2)
    legend(Activity)
    xlabel('Frames')
    ylabel('DF/F')
    title([Groups{g} ' DF/F'])

    subplot(6,length(Groups),g+3*length(Groups)); hold all
    plot(groupAverages.S{g}', 'LineWidth', 2)
    legend(Activity)
    xlabel('Frames')
    ylabel('Spikes')
    title([Groups{g} ' Spikes'])

    subplot(6,length(Groups),[g+length(Groups),g+2*length(Groups)])
    imagesc(groupRasters.F{g}(:,1:end-1))
    ylabel('Cells')
    h = colorbar;
    set(get(h,'label'),'string','DF/F');
    caxis([-0.5, 5]);

    subplot(6,length(Groups),[g+4*length(Groups),g+5*length(Groups)])
    imagesc(groupRasters.S{g}(:,1:end-1))
    ylabel('Cells')
    xlabel('Frames')
    h = colorbar;
    set(get(h,'label'),'string','Spikes');
    caxis([0, 100]);
end

%% PLOT FIGURE BY ACTIVITY TYPE

figure; hold on
suptitle(ExtractedData.StimType)

for a = 1:length(Activity)

    subplot(6,length(Activity),a); hold all
    for g = 1:length(Groups)
        plot(groupAverages.F{g}(a,:), 'LineWidth', 2)
    end
    legend(Groups)
    xlabel('Frames')
    ylabel('DF/F')
    title([Activity{a} ' DF/F'])

    subplot(6,length(Activity),a+3*length(Activity)); hold all
    for g = 1:length(Groups)
        plot(groupAverages.S{g}(a,:), 'LineWidth', 2)
    end
    legend(Groups)
    xlabel('Frames')
    ylabel('Spikes')
    title([Activity{a} ' Spikes'])

    subplot(6,length(Activity),[a+length(Activity),a+2*length(Activity)])
    grp_Raster = [];
    for g = 1:length(Groups)
        grp_Activity  = strcmp([groupActivity{g}{:}], Activity{a})';
        grp_Raster = [grp_Raster; groupRasters.F{g}(grp_Activity,1:end-1)];
    end
    imagesc(grp_Raster)
    ylabel('Cells')
    h = colorbar;
    set(get(h,'label'),'string','DF/F');
    caxis([-0.5, 5]);

    subplot(6,length(Activity),[a+4*length(Activity),a+5*length(Activity)])
    grp_Raster = [];
    for g = 1:length(Groups)
        grp_Activity  = strcmp([groupActivity{g}{:}], Activity{a})';
        grp_Raster = [grp_Raster; groupRasters.S{g}(grp_Activity,1:end-1)];
    end
    imagesc(grp_Raster)
    ylabel('Cells')
    xlabel('Frames')
    h = colorbar;
    set(get(h,'label'),'string','Spikes');
    caxis([0, 100]);
end

%% Save new ExtractedData with cellOrder

if save_data == 1
    cd(save_path)
    filename = ExtractedData.Filename;
    save(filename, 'ExtractedData');
end

%% Helper script

removeOutlier = 0;

outlier = [25];

if removeOutlier
    ExtractedData.AutoActivity(outlier,:) = [];
    ExtractedData.NominalData(outlier,:) = [];
    ExtractedData.NumericalData(outlier,:) = [];
    ExtractedData.Calcium_Raster(outlier,:) = [];
    ExtractedData.Spikes_Raster(outlier,:) = [];
end