% Extract cell-specific data from blocks and save in excel format
% v3 works with Info spreadsheets

%stim protocol code is:
%noiseburst = 1
%ReceptiveField = 2
%FM sweep = 3
%SAM = 6s
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
        info_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study';
        blocks_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\CompiledBlocks_V2';
        save_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\ExtractedData';
        info_filename = 'Info_NxDG021621F2';
       
        stimProtocol = 3; %See stimProtocol list above
        STDlevel = 2; %minimum # standard deviations above baseline to be considered significant
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
        STDlevel = 2; %minimum # standard deviations above baseline to be considered significant
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
        STDlevel = 2; %minimum # standard deviations above baseline to be considered significant
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
for a = 1:size(blockData.setup.mousename,1) %Mice
    mouseID = blockData.setup.mousename{a};
    
    for b = 1:size(blockData.setup.unique_block_names{a},2) %FOVs
        unique_block_name = blockData.setup.unique_block_names{a}(b);
        block = blockData.([mouseID]).([unique_block_name]);
        nCells = length(block.cell_number);
        clear tempCellList;
        clear temp_unique;
        [tempCellList(1:nCells,1)] = deal(block.setup.expt_group);
        [tempCellList(1:nCells,2)] = deal(mouseID);
        [tempCellList(1:nCells,3)] = deal(block.setup.FOV);
        [tempCellList(1:nCells,4)] = deal(blockData.setup.stim);
        [tempCellList(1:nCells,5)] = deal(block.setup.block_filename);
        [tempCellList(1:nCells,6)] = block.cell_number;
        [temp_unique(1:nCells,1)] = deal(unique_block_name);
        
        cellList = [cellList; tempCellList];
        unique_block_names = [unique_block_names; temp_unique];
    end
end

figures = cell(size(cellList,1),21); %Store figures (if save_figures == 1)
activity = cell(size(cellList,1),2); %Auto-determined activity (activated/prolonged/suppressed)
[data_F, data_S] = deal(nan(size(cellList,1),length(numericalColumnHeaders))); %Numerical data
[raster_F, raster_S] = deal(nan(size(cellList,1),76)); %76 is a magic number

if stimProtocol == 2 %RF
    RF.ColumnHeaders = {'Sparseness', 'BF', 'BF_I', 'CF', 'CF_I', 'BestInt', 'Int_RMS', 'Int_halfwidth',...
        'BW_20_RMS', 'BW_20_halfwidth', 'BW_BF_RMS', 'BW_BF_halfwidth', 'ISI', 'Type'};
    st = struct; %temp variable
    st.FRA = {};
    [st.response, st.isResponsive, st.peak, st.peakavg, st.area, RF.stim_v1, RF.stim_v2] = deal(nan(8,8,size(cellList,1))); %8x8 is magic number for RF size
    [st.sp_by_int, st.sp_by_freq] = deal(nan(size(cellList,1),8));
    st.data = nan(size(cellList,1), length(RF.ColumnHeaders));
    RF.F = st;
    RF.S = st;
end

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
        case {3,2} %{'FM','RF'}
            stim_v0 = stim_v2;
            
        case {10,9,11} %{'NoiseITI', 'water', 'air'}
            stim_v0 = stim_v1;
            stim_v2 = zeros(size(stim_v1));
            
        case {6} %{'SAM', 'SAMfreq'} %NAN trials instead of zeros           
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
    F = F7_df_f(stimTrials,:);
    F_blanks = F7_df_f(blankTrials,:);
    S = spks(stimTrials,:);
    S_blanks = spks(blankTrials,:);


    %GET AVERAGED RESPONSES
    %check if each condition is active, then concatenate and keep only active conditions
    nStimConditions = size(unique([stim_v1,stim_v2],'rows'),1); %skip if there's only one stim condition (e.g. NoiseITI)
    analyze_by_stim_condition = 1; %I assume we'll almost always want to do this
    if analyze_by_stim_condition &&  nStimConditions > 1 

        F_by_Stim = [];
        S_by_Stim = [];

        unique_stim_v1 = unique(stim_v1);
        unique_stim_v2 = unique(stim_v2);

        %If RF, order intensities from highest to lowest
        if stimProtocol == 2 %RF
            unique_stim_v2 = flipud(unique_stim_v2);
        end

        for v = 1:length(unique_stim_v1)
            for vv = 1:length(unique_stim_v2)
                stim_rows = intersect(find(stim_v1 == unique_stim_v1(v)), find(stim_v2 == unique_stim_v2(vv)));

                %Plot individual figures of each stim condition
                plotStimFigures = 0;
                if plotStimFigures; figure; end

                for i = 1:2 % F and S
                    if i == 1
                        A = 'F';
                        nFramesToAvg = framerate; %1s
                        trials = F(stim_rows,:);
                        blank_trials = F_blanks;
                        AUClevel = AUC_F_level;
                        units = 'df/f';
                    elseif i == 2
                        A = 'S';
                        nFramesToAvg = round(framerate/2); %0.5s
                        trials = S(stim_rows,:);
                        blank_trials = S_blanks;
                        AUClevel = AUC_S_level;
                        units = 'spikes';
                    end

                    if isempty(trials); continue; end

                    if plotStimFigures; subplot(1,2,i); end

                    [isActive, tempActivity, pk_temp, tr_temp] = checkIfActive_v2(trials, nBaselineFrames, STDlevel, AUClevel, plotStimFigures, units);

                    if isActive
                        if i == 1
                            F_by_Stim = [F_by_Stim; trials];
                        elseif i == 2
                            S_by_Stim = [S_by_Stim; trials];
                        end
                    end

                    if stimProtocol == 2 %RF
                        y = smooth(nanmean(trials,1),3)';
                        response = nanmean(y(1,nBaselineFrames:nBaselineFrames+nFramesToAvg));

                        %Store either peak or trough data
                        if strcmp(tempActivity, 'suppressed')
                            combinedData = tr_temp;
                        else
                            combinedData = pk_temp;
                        end

                        %Compute average of peak or trough region for receptive field
                        if any(combinedData) %skip if all nans
                            c1 = combinedData(2) + nBaselineFrames; %latency and correct for baseline
                            c2 = c1 + combinedData(4); %t1 + width
                            peakavg = nanmean(y(1,c1:c2));
                        else
                            peakavg = nan;
                        end

                        RF.stim_v1(vv,v,c) = unique_stim_v1(v);
                        RF.stim_v2(vv,v,c) = unique_stim_v2(vv);
                        RF.(A).response(vv,v,c) = response;
                        RF.(A).isResponsive(vv,v,c) = isActive;                            
                        RF.(A).peak(vv,v,c) = combinedData(1);
                        RF.(A).peakavg(vv,v,c) = peakavg;
                        RF.(A).area(vv,v,c) = combinedData(5);
                    end
                end
                   if plotStimFigures; suptitle([' Cell ' num2str(c) ': Stim ' num2str(unique_stim_v1(v)) '/' num2str(unique_stim_v2(vv))]); end
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
        avg_F = nan(1,76); %magic number
    end
    if isempty(avg_S)
        avg_S = nan(1,76); %magic number
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
            units = 'df/f';
        elseif i == 2
            trials = S;
            blankTrials = S_blanks;
            AUClevel = AUC_S_level;
            units = 'spikes'; 
        end

        if plot_graphs; subplot(1,2,i);end

        [~, tempActivity, peak_data, trough_data] = checkIfActive_v2(trials, nBaselineFrames, STDlevel, AUClevel, plot_graphs, units);

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
        [RF.F.data(c,1), RF.F.sp_by_int(c,:), RF.F.sp_by_freq(c,:), figures{c,2}] = compute_tuning_sparseness_v2(RF.F.peak(:,:,c), plot_tuning, 'df/f',  unique_stim_v1, unique_stim_v2, cellNumber);
        [RF.S.data(c,1), RF.S.sp_by_int(c,:), RF.S.sp_by_freq(c,:), figures{c,3}] = compute_tuning_sparseness_v2(RF.S.peak(:,:,c), plot_tuning, 'spikes', unique_stim_v1, unique_stim_v2, cellNumber);  

        %Compute frequency tuning
        if any(any(RF.F.isResponsive(:,:,c)))
            [RF.F.FRA{c,1}, RF.F.data(c,2:end), figures{c,4:8}] = ... 
            compute_frequency_tuning_v2(RF.F.isResponsive(:,:,c), RF.F.response(:,:,c),...
            plot_tuning, unique_stim_v1, unique_stim_v2, cellNumber);
        end

        if any(any(RF.S.isResponsive(:,:,c)))
            [RF.S.FRA{c,1}, RF.S.data(c,2:end), figures{c,9:13}] = ... 
            compute_frequency_tuning_v2(RF.S.isResponsive(:,:,c), RF.S.response(:,:,c),...
            plot_tuning, unique_stim_v1, unique_stim_v2, cellNumber);
        end

        if plot_tuning
            %Plot RF using various measures
            figures{c,14} = figure('units','normalized','outerposition',[0 0 0.5 1]); hold on

            subplot(5,2,1)
            imagesc(RF.F.response(:,:,c))
            title('df/f: mean response')
            ylabel('dB')

            subplot(5,2,2)
            imagesc(RF.S.response(:,:,c))
            title('spikes: mean response')

            subplot(5,2,3)
            imagesc(RF.F.peak(:,:,c))
            title('df/f: peak amplitude');
            ylabel('dB')

            subplot(5,2,4)
            imagesc(RF.S.peak(:,:,c))
            title('spikes: peak amplitude')

            subplot(5,2,5)
            imagesc(RF.F.peakavg(:,:,c))
            title('df/f: average around peak');
            ylabel('dB')

            subplot(5,2,6)
            imagesc(RF.S.peakavg(:,:,c))
            title('spikes: average around peak')

            subplot(5,2,7)
            imagesc(RF.F.area(:,:,c))
            title('df/f: response area')
            ylabel('dB')

            subplot(5,2,8)
            imagesc(RF.S.area(:,:,c))
            title('spikes: response area')

            subplot(5,2,9)
            imagesc(RF.F.isResponsive(:,:,c))
            title('df/f: is responsive (0 or 1)')
            xlabel('Frequencies')
            ylabel('dB')

            subplot(5,2,10)
            imagesc(RF.S.isResponsive(:,:,c))
            title('spikes: is responsive (0 or 1)')
            xlabel('Frequencies')

            suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber)))
            set(findobj(gcf,'type','axes'),'XTickLabel',unique_stim_v1,'YTickLabel', unique_stim_v2,...
                   'XTick',1:length(unique_stim_v1), 'YTick',1:length(unique_stim_v2));

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

for g = 1:length(Groups)
    currentCells = strcmpi(groupList,Groups(g));
        
    resorted_raster_F = [];
    resorted_raster_S = [];
    average_F = [];
    average_S = [];
        
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
           
        %Save cell order (makeValidName gets rid of invalid characters like '-' or '/')
        ExtractedData.CellsInOrder.(matlab.lang.makeValidName([Groups{g}])).([Activity{i}]).F = cellNumbers(sort_ind); 
        ExtractedData.CellsInOrder.(matlab.lang.makeValidName([Groups{g}])).([Activity{i}]).S = cellNumbers(sort_ind);
    end
                
    groupRasters.F{g} = resorted_raster_F;
    groupRasters.S{g} = resorted_raster_S;
    groupAverages.F{g} = average_F;
    groupAverages.S{g} = average_S;
end
        

%% PLOT FIGURE

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