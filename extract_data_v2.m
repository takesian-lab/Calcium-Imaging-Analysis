% Extract cell-specific data from blocks and save in excel format
% cellList is a spreadsheet with the following variables:
% 1. Group
% 2. Mouse name
% 3. FOV
% 4. Stim type
% 5. Block filename
% 6. Cell number

%% Set up

clear all

%COLUMN HEADERS FOR FINAL SPREADSHEET
nominalColumnHeaders = {'Group', 'Mouse ID', 'FOV', 'Data type', 'Block', 'Cell Number', 'F-Activity', 'S-Activity'};
numericalColumnHeaders = {'Peak Amplitude', 'P1', 'Peak Latency', 'Peak Width', 'Trough Amplitude', 'T1', 'Trough Latency', 'Trough Width',...
    'Combined Amplitude', 'Combined L1', 'Combined L2', 'Combined Width', 'Peak AUC', 'Trough AUC', 'Combined AUC'};

%% Environment

PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
        cellList_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
        blocks_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
        save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData MET';
        cellList_filename = 'Responsive cells v2';
        %cellList_path = 'Z:\Maryse\2p analysis\Presentations\Joint lab meeting April 7, 2021';
        %blocks_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\CompiledBlocks_v2';
        %save_path = 'Z:\Maryse\2p analysis\Presentations\Joint lab meeting April 7, 2021';
        %cellList_filename = 'Behavior_NDNF';
        %cellList_filename = 'Behavior_Thy1';
        
        stimType = 'RF'; %To look at one stim type at a time
        STDlevel = 2;
        AUC_F_level = 5;
        AUC_spks_level = 10;
        sort_active = 1;  % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 0;
        save_data = 0;
        analyze_by_stim_condition = 1; %determine if cell is active based on individual stim conditions
        
    case 'RD0332' %Carolyn
        cellList_path = 'Z:\Carolyn\2P Imaging data\SSRI response stimuli pilot\VxDD062420M3';
        blocks_path = 'Z:\Carolyn\2P Imaging data\SSRI response stimuli pilot\compiled blocks';
        save_path = 'Z:\Carolyn\2P Imaging data\SSRI response stimuli pilot\VxDD062420M3\extracted';
        cellList_filename = 'VxDD062420M3_postFLX';
        
        stimType = 'FM'; %To look at one stim type at a time
        STDlevel = 2;
        AUC_F_level = 0.05;
        AUC_spks_level = 5;
        sort_active = 0; % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 0;
        save_data = 0;
        analyze_by_stim_condition = 0; %determine if cell is active based on individual stim conditions
        
    otherwise
        disp('Computer does not match known users')
        return
end

%% Load data

cd(cellList_path)
cellList = importfile(cellList_filename);
cellList(1,:) = []; %Remove header

%% Loop through blocks

cd(blocks_path)

%Look at only one stim type at a time
stimTypes = [cellList{:,4}]';
if ~isempty(stimType) %Allows for stimType column to be empty (assuming user has only included one type of stim)
    cellList = cellList(strcmpi(stimTypes, stimType),:);
end

%INITIALIZE VARIABLES FOR FINAL SPREADSHEET
data = {}; %Nominal data
activity = cell(size(cellList,1),2); %Auto-determined activity (activated/prolonged/suppressed)
[data_F, data_S] = deal(nan(size(cellList,1),length(numericalColumnHeaders))); %Numerical data
[raster_F, raster_S] = deal(nan(size(cellList,1),76)); %76 is a magic number

%Loop through all cells in each block
blocks = [cellList{:,5}]';
uniqueBlocks = unique(blocks);
count = 1;

for b = 1:length(uniqueBlocks)
    currentBlock = uniqueBlocks{b};
    block_cellList = cellList(strcmpi(blocks, currentBlock),:);
    load(currentBlock);
    baseline_length = block.setup.constant.baseline_length; %seconds
    framerate = block.setup.framerate;
    nBaselineFrames = baseline_length*framerate; %frames
    
    data = [data; block_cellList(:,1:6)];
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
    switch stimType
        case {'FM','RF'}
            stim_v0 = stim_v2;
            
        case {'NoiseITI', 'water', 'air'}
            %Regular noise is not accounted for yet
            stim_v0 = stim_v1;
            stim_v2 = zeros(size(stim_v1));
            
        case {'SAM', 'SAMfreq'} %NAN trials instead of zeros           
            stim_v0 = stim_v1;
            stim_v0(isnan(stim_v0)) = 0;
            
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

    %Store stim values in case they change for some cells
    store_stim_v1 = stim_v1;
    store_stim_v2 = stim_v2;
    store_stimTrials = stimTrials;
    store_blankTrials = blankTrials;
    
    for c = 1:size(block_cellList,1)
        cellNumber = block_cellList{c,6};
        if cellNumber == 'NaN'
            count = count + 1;
            continue
        end

        %when we remove inf below stim might change so refresh it with original stim list
        stim_v1 = store_stim_v1;
        stim_v2 = store_stim_v2;
        stimTrials = store_stimTrials;
        blankTrials = store_blankTrials;
            
        cellIndex = find(block.cell_number == cellNumber);
        if isempty(cellIndex)
            warning(['Cell number ' num2str(cellNumber) ' not found.'])
            count = count + 1;
            continue;
        end
                
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

        
        %GET AVERAGED AND SMOOTHED RESPONSES
        %check if each condition is active, then concatenate and keep only active conditions
        nStimConditions = size(unique([stim_v1,stim_v2],'rows'),1); %skip if there's only one stim condition (e.g. NoiseITI)
        if analyze_by_stim_condition &&  nStimConditions > 1 

            F_by_Stim = [];
            S_by_Stim = [];
            
            unique_stim_v1 = unique(stim_v1);
            unique_stim_v2 = unique(stim_v2);
            
            for v = 1:length(unique_stim_v1)
                for vv = 1:length(unique_stim_v2)
                    stim_rows = intersect(find(stim_v1 == unique_stim_v1(v)), find(stim_v2 == unique_stim_v2(vv)));
                    F_temp = F(stim_rows,:);
                    S_temp = S(stim_rows,:);

                    F_temp_smoothed = smooth(nanmean(F_temp,1),3)';
                    S_temp_smoothed = smooth(nanmean(S_temp,1),3)';
                    [F_active, ~] = checkIfActive_v2(F_temp_smoothed, nBaselineFrames, STDlevel, AUC_F_level, 0);
                    [S_active, ~] = checkIfActive_v2(S_temp_smoothed, nBaselineFrames, STDlevel, AUC_F_level, 0);
                    
                    if F_active
                        F_by_Stim = [F_by_Stim; F_temp];
                    end
                    
                    if S_active
                        S_by_Stim = [S_by_Stim; S_temp];
                    end
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

        avg_F = smooth(nanmean(F,1),3)';
        avg_S = smooth(nanmean(S,1),3)';
        
        %FILL RASTER
        %Data might be empty if too many trials were removed due to loco activity
        if isempty(avg_F)
            avg_F = nan(1,76); %magic number
        end
        if isempty(avg_S)
            avg_S = nan(1,76); %magic number
        end
        raster_F(count,1:length(avg_F)) = avg_F;
        raster_S(count,1:length(avg_S)) = avg_S;
                
        %COMPUTE LATENCIES, AMPLITUDES, AND WIDTHS
        if plot_graphs == 1; figure; hold on; end
        
        for i = 1:2 %F and S
            if i == 1
                y = avg_F;
                AUC_level = AUC_F_level;
                units = 'DF/F';               
            elseif i == 2
                y = avg_S;
                AUC_level = AUC_spks_level;
                units = 'Deconvolved spikes';
            end
                    
            baseline = y(1,1:nBaselineFrames);
            response = y(1,nBaselineFrames+1:end);
            peak_threshold = nanmean(baseline) + STDlevel*std(baseline);
            trough_threshold = nanmean(baseline) - STDlevel*std(baseline);
                                
            %PEAK COMPUTATIONS
            peak_data = nan(1,4);
            [peak, peak_latency] = max(response);
            if peak >= peak_threshold && any(response) %only store data if peak is above threshold
                [p1_latency] = find(response >= peak_threshold, 1, 'first');
                [p2_latency] = find(response(1, peak_latency:end) <= peak_threshold, 1, 'first') - 2;
                p1 = response(p1_latency);
                p2 = response(peak_latency + p2_latency);

                %AUC
                if isempty(p2_latency)
                    p2_latency_temp = length(response);
                else
                    p2_latency_temp = p2_latency + peak_latency;
                end
                peak_trace = response(1,p1_latency:p2_latency_temp);
                peak_trace(peak_trace < peak_threshold) = peak_threshold;
                peak_trace_no_nan = peak_trace(~isnan(peak_trace)); %trapz function does not work on nans
                aup = trapz(abs(peak_trace_no_nan - peak_threshold)); %Area under peak above threshold

                %Adjust for baseline
                peak_latency = peak_latency + nBaselineFrames;
                p1_latency = p1_latency + nBaselineFrames;
                p2_latency = p2_latency + peak_latency;

                %Width
                peak_width = p2_latency - p1_latency;
            else
                [peak, p1, p2, p1_latency, p2_latency, peak_latency, peak_width, aup] = deal(nan);
            end

            %Store
            if ~isempty(peak);          peak_data(1) = peak;            else;   peak = nan;         end
            if ~isempty(p1_latency);    peak_data(2) = p1_latency;      else;   p1_latency = nan;   end
            if ~isempty(peak_latency);  peak_data(3) = peak_latency;    else;   peak_latency = nan; end
            if ~isempty(peak_width);    peak_data(4) = peak_width;      else;   peak_width = nan;   end
                    
            %TROUGH COMPUTATIONS
            trough_data = nan(1,4);
            [trough, trough_latency] = min(response);
            if trough <= trough_threshold %only store data if trough is below threshold
                [t1_latency] = find(response <= trough_threshold, 1, 'first');
                [t2_latency] = find(response(1, trough_latency:end) >= trough_threshold, 1, 'first') - 2;
                t1 = response(t1_latency);
                try
                    t2 = response(trough_latency + t2_latency);
                    t2_latency_temp = t2_latency + trough_latency;
                catch
                    t2_latency_temp = length(response);
                end

                %AUC
                %                 if isempty(t2_latency)
                %                     t2_latency_temp = length(response);
                %                 else
                %                     t2_latency_temp = t2_latency + trough_latency;
                %                 end
                
                trough_trace = response(1,t1_latency:t2_latency_temp);
                trough_trace(trough_trace > trough_threshold) = trough_threshold;
                trough_trace_no_nan = trough_trace(~isnan(trough_trace));
                aat = trapz(abs(trough_trace_no_nan - trough_threshold)); %Area above trough and below threshold

                %Adjust for baseline
                trough_latency = trough_latency + nBaselineFrames;
                t1_latency = t1_latency + nBaselineFrames;
                t2_latency = t2_latency + trough_latency;

                %Width
                trough_width = t2_latency - t1_latency;
            else
                [trough, t1, t2, t1_latency, t2_latency, trough_latency, trough_width, aat] = deal(nan);
            end
                    
            %Store
            if ~isempty(trough);            trough_data(1) = trough;            else;   trough = nan;           end
            if ~isempty(t1_latency);        trough_data(2) = t1_latency;        else;   t1_latency = nan;       end
            if ~isempty(trough_latency);    trough_data(3) = trough_latency;    else;   trough_latency = nan;   end
            if ~isempty(trough_width);      trough_data(4) = trough_width;      else;   trough_width = nan;     end

            %Auto-determined activity (suppressed/prolonged/activated)
            if ~isnan(aup)&& aup >= AUC_level; aup_pass = true; else; aup_pass = false; end
            if ~isnan(aat)&& aat >= AUC_level; aat_pass = true; else; aat_pass = false; end

            
            tempActivity = 'none';

            if isnan(peak) && isnan(trough)
                tempActivity = 'none';
            elseif ~aat_pass && ~aup_pass
                tempActivity = 'none';
            elseif isnan(peak) && ~isnan(trough) && aat_pass
                tempActivity = 'suppressed';
            elseif ~isnan(peak) && isnan(trough) && aup_pass
                if peak_latency > 40 || isempty(p2_latency)
                    tempActivity = 'prolonged';
                else
                    tempActivity = 'activated';
                end
            elseif ~isnan(peak) && ~isnan(trough)
                if (trough_latency < peak_latency) && aat_pass
                    tempActivity = 'suppressed';
                elseif (peak_latency < trough_latency) && aat_pass && ~aup_pass
                    tempActivity = 'suppressed';
                elseif aup_pass
                    if peak_latency > 40 || isempty(p2_latency)
                        tempActivity = 'prolonged';
                    else
                        tempActivity = 'activated';
                    end
                else
                    tempActivity = 'none';
                end
            else
                tempActivity = 'none';
            end

            activity{count,i} = tempActivity;
                    
            %Plot
            if plot_graphs == 1 && any(y)
                subplot(1,2,i); hold on
                plot(y)
                hline(nanmean(baseline), 'k')
                hline(peak_threshold, 'r')
                hline(trough_threshold, 'c')
                scatter(peak_latency, peak, 'o', 'r')
                scatter(p1_latency, p1, 'o', 'r')
                scatter(p2_latency, p2, 'o', 'r')
                scatter(trough_latency, trough, 'o', 'c')
                scatter(t1_latency, t1, 'o', 'c')
                scatter(t2_latency, t2, 'o', 'c')
                vline(nBaselineFrames, 'k')
                xlabel('Frames')
                ylabel(units)
                if strcmpi(tempActivity, 'activated') || strcmpi(tempActivity, 'prolonged')
                    title([tempActivity ' -  AUC: ' num2str(aup)])
                    plot(p1_latency:(nBaselineFrames + p2_latency_temp), peak_trace, 'g')
                elseif strcmpi(activity, 'suppressed')
                    title([tempActivity ' - AUC: ' num2str(aat)])
                    plot(t1_latency:(nBaselineFrames + t2_latency_temp), trough_trace, 'g')
                else
                    title(tempActivity)
                end

                if i == 2
                    suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber)))
                end
            end

            %Convert to seconds and store data
            peak_data = peak_data./framerate;
            trough_data = trough_data./framerate;

            %For storing combined data, only keep trough OR peak data in one column
            if isequal(tempActivity, 'activated') || isequal(tempActivity, 'prolonged')
                combined_data = peak_data;
                combined_auc = aup;
            elseif isequal(tempActivity, 'suppressed')
                combined_data = trough_data;
                combined_auc = aat;
            else
                combined_data = nan(1,4);
                combined_auc = nan;
            end

            if i == 1 %GCaMP
                data_F(count,:) = [peak_data, trough_data, combined_data, aup, aat, combined_auc];
            elseif i == 2 %Spikes
                data_S(count,:) = [peak_data, trough_data, combined_data, aup, aat, combined_auc];
            end

        end
        count = count + 1;
    end
end
    
ExtractedData = struct;
ExtractedData.StimType = stimType;
ExtractedData.STDlevel = STDlevel;
ExtractedData.AUC_F_level = AUC_F_level;
ExtractedData.AUC_spks_level = AUC_spks_level;
ExtractedData.Sort_Active = sort_active;
ExtractedData.Analyze_By_Stim_Condition = analyze_by_stim_condition;
ExtractedData.ColumnHeaders = [nominalColumnHeaders, strcat('F-',numericalColumnHeaders), strcat('S-',numericalColumnHeaders)];
ExtractedData.Activity = activity;
ExtractedData.NominalData = data;
ExtractedData.NumericalData = [data_F,data_S];
ExtractedData.Calcium_Raster = raster_F;
ExtractedData.Spikes_Raster = raster_S;

%% Save extracted data
if save_data == 1
    cd(save_path)
    d = datestr(now,'yyyymmdd-HHMMSS');
    save(['extractedData_' stimType '_' d '.mat'], 'ExtractedData');
%   save(['extractedData_' stimType '.mat'], 'ExtractedData');
end
%   save(['motorSort_list_' dataType '.mat'], 'CellsInOrder');

%% Plot sorted rasters
    
%Plot by activity type and peak/trough amplitude
    
Groups = {'VIP', 'NDNF'};
tempActivity = {'activated', 'prolonged', 'suppressed'};
    
cellList = [ExtractedData.NominalData{:,1}]'; %Modify
activityList = [ExtractedData.Activity(:,1)]; %Modify
cellOrder = [];

for g = 1:length(Groups)
    currentCells = strcmpi(cellList,Groups{g});
        
    resorted_raster_F = [];
    resorted_raster_S = [];
    average_F = [];
    average_S = [];
        
    for i = 1:length(tempActivity)       
        currentActivity = tempActivity{i};
        activeRows = strcmpi(activityList, currentActivity);
        currentRows = and(currentCells, activeRows);
        cellNumbers = find(currentRows); %cell numbers to use for comparison in locomotor 

        %find rows in current activity type and sort by amplitude and/or latency
        current_raster_F = ExtractedData.Calcium_Raster(currentRows,:);
        current_raster_S = ExtractedData.Spikes_Raster(currentRows,:);

        %Store average for plots
        average_F(i,:) = nanmean(current_raster_F,1);
        average_S(i,:) = nanmean(current_raster_S,1);

        if i == 1 || i == 2 %Peaks for activated and prolonged
            current_F_amplitude = ExtractedData.NumericalData(currentRows,1);
            current_S_amplitude = ExtractedData.NumericalData(currentRows,16);
        elseif i == 3 %Troughs for suppressed
            current_F_amplitude = ExtractedData.NumericalData(currentRows,5);
            current_S_amplitude = ExtractedData.NumericalData(currentRows,20);
        end
            
        [~, F_sort_ind] = sort(current_F_amplitude, 'descend');
        [~, S_sort_ind] = sort(current_S_amplitude, 'descend');

        %SORT BOTH BY GCAMP (so they match)
        resorted_raster_F = [resorted_raster_F; current_raster_F(F_sort_ind,:)];
        resorted_raster_S = [resorted_raster_S; current_raster_S(F_sort_ind,:)];

        %Save cell order
        cellOrder = [cellOrder; cellNumbers(F_sort_ind)];
        CellsInOrder.([Groups{g}]).([tempActivity{i}]) = cellNumbers(F_sort_ind);
    end
            
    if g == 1
        VIP_raster_F = resorted_raster_F;
        VIP_raster_S = resorted_raster_S;
        VIP_average_F = average_F;
        VIP_average_S = average_S;
        %VIP_cell_order = cellOrder;
    elseif g == 2
        NDNF_raster_F = resorted_raster_F;
        NDNF_raster_S = resorted_raster_S;
        NDNF_average_F = average_F;
        NDNF_average_S = average_S;
        %NDNF_cell_order = cellOrder;
    end
    
end
        
       
smoothCurves = 1;

if smoothCurves
    for s = 1:size(VIP_average_F,1)
        VIP_average_F(s,:) = smooth(VIP_average_F(s,:));
    end
    for s = 1:size(NDNF_average_F,1)
        NDNF_average_F(s,:) = smooth(NDNF_average_F(s,:));
    end
end

figure; hold on
suptitle(ExtractedData.StimType)

subplot(6,2,1); hold all
plot(VIP_average_F', 'LineWidth', 2)
legend(tempActivity)
xlabel('Frames')
ylabel('DF/F')
title('VIP DF/F')

subplot(6,2,2); hold all
plot(VIP_average_S', 'LineWidth', 2)
legend(tempActivity)
xlabel('Frames')
ylabel('Spikes')
title('VIP Spikes')

subplot(6,2,[3,5])
imagesc(VIP_raster_F(:,1:end-1))
ylabel('Cells')
%xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','DF/F');
caxis([-0.5, 5]);

subplot(6,2,[4,6])
imagesc(VIP_raster_S(:,1:end-1))
ylabel('Cells')
%xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','Spikes');
caxis([0, 100]);

subplot(6,2,7)
plot(NDNF_average_F', 'LineWidth', 2)
legend(tempActivity)
xlabel('Frames')
ylabel('DF/F')
title('NDNF DF/F')

subplot(6,2,8)
plot(NDNF_average_S', 'LineWidth', 2)
legend(tempActivity)
xlabel('Frames')
ylabel('Spikes')
title('NDNF Spikes')

subplot(6,2,[9,11])
imagesc(NDNF_raster_F(:,1:end-1))
ylabel('Cells')
xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','DF/F');
caxis([-0.5, 5]);

subplot(6,2,[10,12])
imagesc(NDNF_raster_S(:,1:end-1))
ylabel('Cells')
xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','Spikes');
caxis([0, 100]);

%% Save cell order
if save_data == 1
    cd(save_path)
    d = datestr(now,'yyyymmdd-HHMMSS');
    save(['CellsInOrder_' stimType '_' d '.mat'], 'CellsInOrder');
    %save(['extractedData_' dataType '.mat'], 'ExtractedData');
end
    %save(['motorSort_list_' dataType '.mat'], 'CellsInOrder');

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