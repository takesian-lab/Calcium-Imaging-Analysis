% Extract cell-specific data from blocks and save in excel format
% cellList is a spreadsheet with the following variables:
% 1. Group
% 2. Mouse name
% 3. FOV
% 4. Data type
% 5. Block filename
% 6. Cell number
% 7. Activity type

%% Set up

clear all

columnHeaders = {'Group', 'Mouse ID', 'FOV', 'Data type', 'Block', 'Cell Number', 'Activity', 'Auto-Activity GCamP', 'Auto-Activity Spikes',...
    'GCaMP Peak Amplitude', 'GCaMP P1', 'GCaMP Peak Latency', 'GCaMP Peak Width',...
    'GCaMP Trough Amplitude', 'GCaMP T1', 'GCaMP Trough Latency', 'GCaMP Trough Width',...
    'Spike Peak Amplitude', 'Spike P1', 'Spike Peak Latency', 'Spike Peak Width',...
    'Spike Trough Amplitude', 'Spike T1', 'Spike Trough Latency', 'Spike Trough Width'...
    'Combined GCaMP Amplitude', 'Combined GCaMP L1', 'Combined GCamP L2', 'Combined GCamP Width',...
    'Combined Spike Amplitude', 'Combined Spike L1', 'Combined Spike L2', 'Combined Spike Width'...
    'GCaMP Peak AUC', 'GCaMP Trough AUC', 'Spike Peak AUC', 'Spike Trough AUC', 'Combined GCaMP AUC', 'Combined Spike AUC'};

RF_columnHeaders = {'CF', 'Threshold', 'BW10', 'BW20'};

%% Environment

PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
        cellList_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
        blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
        save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData MET';
        cellList_filename = 'ResponsiveCells';
        
        dataType = 'RF';
        STDlevel = 2;
        AUC_F_level = 0.05;
        AUC_spks_level = 5;
        sort_active = 1;  % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 0; %Plot one graph for each cell
        plot_RF = 1; %Plot receptive fields for each cell
        save_data = 1;
        analyze_by_stim_condition = 1; %determine if cell is active based on individual stim conditions
        
    case 'RD0332' %Carolyn
        cellList_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS';
        blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
        save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS\Locomotor_receptivecell only';
        cellList_filename = 'ResponsiveCellsOnly';
        
        dataType = 'FM'; %To look at one stim type at a time. Leave empty to look at all
        STDlevel = 2;
        AUC_F_level = 0.05;
        AUC_spks_level = 5;
        sort_active = 2; % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials
        plot_graphs = 0;
        save_data = 1;
        analyze_by_stim_condition = 1; %determine if cell is active based on individual stim conditions
        
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
dataTypes = [cellList{:,4}]';
if ~isempty(dataType)
    cellList = cellList(strcmpi(dataTypes, dataType),:);
end

data1 = {}; %Nominal data
autoActivity = cell(size(cellList,1),2); %Auto-determined activity (inhibited/sustained/activated)
data2 = nan(size(cellList,1),30); %Numerical data
raster_F = nan(size(cellList,1),76);
raster_spks = nan(size(cellList,1),76);
RF_F = nan(8,8,size(cellList,1));
RF_isResponsive_F = nan(8,8,size(cellList,1));
RF_spks = nan(8,8,size(cellList,1));
RF_isResponsive_spks = nan(8,8,size(cellList,1));

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
    response_window = block.setup.constant.response_window;
    response_window_in_frames = response_window*framerate;

    
    data1 = [data1; block_cellList(:,1:7)];
    stim_v1 = block.parameters.variable1';
    stim_v2 = block.parameters.variable2';
    
    
    % identify trials to remove (running trials or zero sound trials)
    if strcmpi(dataType,'RF')
        if sort_active ~= 0
            r = find(stim_v2 == 0); %find 0dB trials
            % locomotor or non-loco trials
            if sort_active==1
                rr = find(block.active_trials==1);%find active trials to remove
            elseif sort_active==2
                rr = find(block.active_trials==0);%find inactive trials to remove
            end
            
            ru = union(r,rr);%put the two lists together
            remove = ru; %Making a variable that will be the same across all stim.
            stim_v1(remove,:) = [];
            stim_v2(remove,:) = [];
        else
            remove = find(stim_v2 == 0);
            stim_v1(stim_v2 == 0,:) = [];
            stim_v2(stim_v2 == 0,:) = [];
        end
    else
        error('Stim types other than RF not supported')
    end
    
    store_stim_v1 = stim_v1;
    store_stim_v2 = stim_v2;
    
    for c = 1:size(block_cellList,1)
        cellNumber = block_cellList{c,6};
        
        %when we remove inf below stim might change so refresh it with original stim list
        stim_v1 = store_stim_v1;
        stim_v2 = store_stim_v2;
        
        if isstring(cellNumber)
            %Cell was Nan but it's a string so you can't use isnan function
            autoActivity(count,1:2) = {'NaN'};
            count = count + 1;
            continue; %Go to next cell            
        elseif isnan(cellNumber) || isempty(stim_v1)
            %Cell was not match or there are no active/inactive trials
            autoActivity(count,1:2) = {'NaN'};
            count = count + 1;
            continue; %Go to next cell
        end
        
        cellIndex = find(block.cell_number == cellNumber);
        if isempty(cellIndex)
            warning(['Cell number ' num2str(cellNumber) ' not found.'])
            autoActivity(count,1:2) = {'NaN'};
            count = count + 1;
            continue; %Go to next cell
        end

        %Cell-specific response to sound
        F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:));
        F7_baseline = F7(:,1:nBaselineFrames); %baseline for each trial
        F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
        spks = squeeze(block.aligned_stim.spks_stim(cellIndex,:,:));

        %Eliminate 0dB trials -> If they exist we could potentially use
        %them to confirm above-baseline responses
        %If desired, remove active trials here too. They were identified for the
        %whole block above. Now we're taking them out of each trace.
        F7_df_f(remove,:) = [];
        spks(remove,:) = [];
        
        %Eliminate trials with infinite values
        [inf_rows,~] = find(isinf(F7_df_f));
        remove_inf = unique(inf_rows);
        F7_df_f(remove_inf,:) = [];
        spks(remove_inf,:) = [];
        stim_v1(remove_inf,:) = [];
        stim_v2(remove_inf,:) = [];

        %Get averaged & smoothed response
        if analyze_by_stim_condition %check if each condition is active, then concatenate and keep only active conditions

            F7_by_Stim = [];
            spks_by_Stim = [];

            unique_stim_v1 = unique(stim_v1);
            unique_stim_v2 = flipud(unique(stim_v2));
            stimMat_v1 = [];
            stimMat_v2 = [];
            for v = 1:length(unique_stim_v2) %frequencies
                for vv = 1:length(unique_stim_v1) %intensities
                    stim_rows = intersect(find(stim_v1 == unique_stim_v1(vv)), find(stim_v2 == unique_stim_v2(v)));
                    F7_temp = F7_df_f(stim_rows,:);
                    spks_temp = spks(stim_rows,:);

                    %Average and smooth traces
                    if size(F7_temp,1) == 1 || size(spks_temp,1) == 1 %If only one valid trial
                        F7_temp_smoothed = smooth(F7_temp,3)';
                        spks_temp_smoothed = smooth(spks_temp,3)';
                    else
                        F7_temp_smoothed = smooth(nanmean(F7_temp),3)';
                        spks_temp_smoothed = smooth(nanmean(spks_temp),3)';
                    end

                    %Check if active and add to stim raster
                    RF_isResponsive_F(v,vv,count) = 0;
                    RF_isResponsive_spks(v,vv,count) = 0;
                    
                    [active, tempActivity] = checkIfActive(F7_temp_smoothed, nBaselineFrames, STDlevel, AUC_F_level, 0);
                    if active
                        F7_by_Stim = [F7_by_Stim; F7_temp];
                        if ~isequal(tempActivity,'inhibited')
                            RF_isResponsive_F(v,vv,count) = 1;
                        end
                    end
                    
                    [active, tempActivity] = checkIfActive(spks_temp_smoothed, nBaselineFrames, STDlevel, AUC_F_level, 0);
                    if active
                        spks_by_Stim = [spks_by_Stim; spks_temp];
                        if ~isequal(tempActivity,'inhibited')
                            RF_isResponsive_spks(v,vv,count) = 1;
                        end
                    end
                    
                    %Save averaged activity in RF
                    stimMat_v1(v,vv) = unique_stim_v1(v);
                    stimMat_v2(v,vv) = unique_stim_v2(vv);
                    RF_F(v,vv,count) = nanmean(F7_temp_smoothed(1,nBaselineFrames:response_window_in_frames));
                    RF_spks(v,vv,count) = nanmean(spks_temp_smoothed(1,nBaselineFrames:response_window_in_frames));
                end
            end
            
            if plot_RF
                figure; hold on
                subplot(2,2,1)
                imagesc(RF_F(:,:,count))
                ylabel('dB')
                
                subplot(2,2,2)
                imagesc(RF_spks(:,:,count))
                ylabel('dB')
                
                subplot(2,2,3)
                imagesc(RF_isResponsive_F(:,:,count))
                xlabel('Frequencies')
                
                subplot(2,2,4)
                imagesc(RF_isResponsive_spks(:,:,count))
                xlabel('Frequencies')
                
                suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber)))
            end
            
            if isempty(F7_by_Stim)
                avg_F7_df_f = smooth(nanmean(F7_df_f),3)'; %Take all trials
            elseif size(F7_by_Stim,1) == 1
                avg_F7_df_f = smooth(F7_by_Stim,3)';
            else
                avg_F7_df_f = smooth(nanmean(F7_by_Stim),3)';
            end

            %Spks could have a different number of valid trials than F7
            if isempty(spks_by_Stim)
                avg_spks = smooth(nanmean(spks),3)'; %Take all trials
            elseif size(F7_by_Stim,1) == 1
                avg_spks = smooth(spks_by_Stim,3)';
            else
                avg_spks = smooth(nanmean(spks_by_Stim),3)';
            end

        else
            avg_F7_df_f = smooth(nanmean(F7_df_f),3)';
            avg_spks = smooth(nanmean(spks),3)';
        end

        %Store raster
        if isempty(avg_F7_df_f)
            avg_F7_df_f = nan(1,76);
        end
        if isempty(avg_spks)
            avg_spks = nan(1,76);
        end
        raster_F(count,1:length(avg_F7_df_f)) = avg_F7_df_f;
        raster_spks(count,1:length(avg_spks)) = avg_spks;

        %Compute latencies, width, and amplitude
        if plot_graphs == 1; figure; hold on; end
        for i = 1:2
            if i == 1
                y = avg_F7_df_f;
                AUC_level = AUC_F_level;
                units = 'DF/F';
            elseif i == 2
                y = avg_spks;
                AUC_level = AUC_spks_level;
                units = 'Deconvolved spikes';
            end

            baseline = y(1,1:nBaselineFrames);
            peak_threshold = nanmean(baseline) + STDlevel*std(baseline);
            trough_threshold = nanmean(baseline) - STDlevel*std(baseline);
            response = y(1,nBaselineFrames+1:end);

            %PEAK COMPUTATIONS
            peak_data = nan(1,4);
            [peak, peak_latency] = max(response);
            if peak >= peak_threshold %only store data if peak is above threshold
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
                if isinf(aup)
                    error('Infinite values found')
                end

                %Adjust for baseline
                peak_latency = peak_latency + nBaselineFrames;
                p1_latency = p1_latency + nBaselineFrames;
                p2_latency = p2_latency + peak_latency;

                %Width
                peak_width = p2_latency - p1_latency;
            else
                peak = nan;
                p1 = nan;
                p2 = nan;
                p1_latency = nan;
                p2_latency = nan;
                peak_latency = nan;
                peak_width = nan;
                aup = nan;
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
                t2 = response(trough_latency + t2_latency);

                %AUC
                if isempty(t2_latency)
                    t2_latency_temp = length(response);
                else
                    t2_latency_temp = t2_latency + trough_latency;
                end
                trough_trace = response(1,t1_latency:t2_latency_temp);
                trough_trace(trough_trace > trough_threshold) = trough_threshold;
                trough_trace_no_nan = trough_trace(~isnan(trough_trace));
                aat = trapz(abs(trough_trace_no_nan - trough_threshold)); %Area above trough and below threshold
                if isinf(aat)
                    error('Infinite values found')
                end
                
                %Adjust for baseline
                trough_latency = trough_latency + nBaselineFrames;
                t1_latency = t1_latency + nBaselineFrames;
                t2_latency = t2_latency + trough_latency;

                %Width
                trough_width = t2_latency - t1_latency;
            else
                trough = nan;
                t1 = nan;
                t2 = nan;
                t1_latency = nan;
                t2_latency = nan;
                trough_latency = nan;
                trough_width = nan;
                aat = nan;
            end

            %Store
            if ~isempty(trough);            trough_data(1) = trough;            else;   trough = nan;           end
            if ~isempty(t1_latency);        trough_data(2) = t1_latency;        else;   t1_latency = nan;       end
            if ~isempty(trough_latency);    trough_data(3) = trough_latency;    else;   trough_latency = nan;   end
            if ~isempty(trough_width);      trough_data(4) = trough_width;      else;   trough_width = nan;     end

            %Auto-determined activity (inhibited/sustained/activated)
            if ~isnan(aup)&& aup >= AUC_level; aup_pass = true; else; aup_pass = false; end
            if ~isnan(aat)&& aat >= AUC_level; aat_pass = true; else; aat_pass = false; end

            activity = 'undetermined'; %If it somehow makes it through the conditions without being classified

            if isnan(peak) && isnan(trough)
                activity = 'none';
            elseif ~aat_pass && ~aup_pass
                activity = 'none';
            elseif isnan(peak) && ~isnan(trough) && aat_pass
                activity = 'inhibited';
            elseif ~isnan(peak) && isnan(trough) && aup_pass
                if peak_latency > 40 || isempty(p2_latency)
                    activity = 'sustained';
                else
                    activity = 'activated';
                end
            elseif ~isnan(peak) && ~isnan(trough)
                if (trough_latency < peak_latency) && aat_pass
                    activity = 'inhibited';
                elseif (peak_latency < trough_latency) && aat_pass && ~aup_pass
                    activity = 'inhibited';
                elseif aup_pass
                    if peak_latency > 40 || isempty(p2_latency)
                        activity = 'sustained';
                    else
                        activity = 'activated';
                    end
                else
                    activity = 'none';
                end
            else
                activity = 'none';
            end

            autoActivity{count,i} = activity;

            %Plot
            if plot_graphs == 1
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
                if strcmpi(activity, 'activated') || strcmpi(activity, 'sustained')
                    title([activity ' -  AUC: ' num2str(aup)])
                    plot(p1_latency:(nBaselineFrames + p2_latency_temp), peak_trace, 'g')
                elseif strcmpi(activity, 'inhibited')
                    title([activity ' - AUC: ' num2str(aat)])
                    plot(t1_latency:(nBaselineFrames + t2_latency_temp), trough_trace, 'g')
                else
                    title(activity)
                end

                if i == 2
                    suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber)))
                end
            end

            %Convert to seconds and store data
            peak_data = peak_data./framerate;
            trough_data = trough_data./framerate;

            %For storing combined data, only keep trough OR peak data in one column
            if isequal(activity, 'activated') || isequal(activity, 'sustained')
                combined_data = peak_data;
                combined_auc = aup;
            elseif isequal(activity, 'inhibited')
                combined_data = trough_data;
                combined_auc = aat;
            else
                combined_data = nan(1,4);
                combined_auc = nan;
            end

            if i == 1
                data2(count,1:8) = [peak_data, trough_data]; %GCaMP
                data2(count,17:20) = combined_data; %GCaMP
                data2(count,25:26) = [aup, aat];
                data2(count,29) = combined_auc;
            elseif i == 2
                data2(count,9:16) = [peak_data, trough_data]; %Spikes
                data2(count,21:24) = combined_data; %Spikes
                data2(count,27:28) = [aup, aat];
                data2(count,30) = combined_auc;
            end

        end
        count = count + 1;
    end
end
    
    ExtractedData = struct;
    ExtractedData.DataType = dataType;
    ExtractedData.STDlevel = STDlevel;
    ExtractedData.AUC_F_level = AUC_F_level;
    ExtractedData.AUC_spks_level = AUC_spks_level;
    ExtractedData.Sort_Active = sort_active;
    ExtractedData.Analyze_By_Stim_Condition = analyze_by_stim_condition;
    ExtractedData.ColumnHeaders = columnHeaders;
    ExtractedData.AutoActivity = autoActivity;
    ExtractedData.NominalData = data1;
    ExtractedData.NumericalData = data2;
    ExtractedData.Calcium_Raster = raster_F;
    ExtractedData.Spikes_Raster = raster_spks;
    ExtractedData.RF_F = RF_F;
    ExtractedData.RF_isResponsive_F = RF_isResponsive_F;
    ExtractedData.RF_spks = RF_spks;
    ExtractedData.RF_isResponsive_spks = RF_isResponsive_spks;
    
    %% Save
    if save_data == 1
        cd(save_path)
        d = datestr(now,'yyyymmdd-HHMMSS');
        save(['extractedData_' dataType '_' d '.mat'], 'ExtractedData');
    end
    
%% Plot sorted rasters

%Plot by activity type and peak/trough amplitude

Cells = {'VIP', 'NDNF'};
A = {'activated', 'sustained', 'inhibited'};

cellList = [ExtractedData.NominalData{:,1}]'; %Modify
activityList = [ExtractedData.AutoActivity(:,1)]; %Modify


for c = 1:2 %ndnf vs vip
    currentCells = strcmpi(cellList,Cells{c});

    resorted_raster_F = [];
    resorted_raster_spks = [];
    average_F = [];
    average_spks = [];
    cellOrder = [];

    for i = 1:length(A)
        currentActivity = A{i};
        activeRows = strcmpi(activityList, currentActivity);
        currentRows = and(currentCells, activeRows);
        cellnumbers_A = find(currentRows==1); %cell numbers to use for comparison in locomotor 

        %find rows in current activity type and sort by amplitude and/or latency
        current_raster_F = ExtractedData.Calcium_Raster(currentRows,:);
        current_raster_spks = ExtractedData.Spikes_Raster(currentRows,:);

        %Store average for plots
        if size(current_raster_F,1) == 1
            %Edge case where n = 1 cell
            average_F(i,:) = current_raster_F;
            average_spks(i,:) = current_raster_spks;
        else
            average_F(i,:) = nanmean(current_raster_F);
            average_spks(i,:) = nanmean(current_raster_spks);
        end

        if i == 1 || i == 2 %Peaks for activated and sustained
            current_gcamp_amplitude = ExtractedData.NumericalData(currentRows,1);
            current_spike_amplitude = ExtractedData.NumericalData(currentRows,9);
        elseif i == 3 %Troughs for inhibited
            current_gcamp_amplitude = ExtractedData.NumericalData(currentRows,1);
            current_spike_amplitude = ExtractedData.NumericalData(currentRows,9);
        end

        [~, gcamp_sort_ind] = sort(current_gcamp_amplitude, 'descend');
        [~, spike_sort_ind] = sort(current_spike_amplitude, 'descend');

        %SORT BOTH BY GCAMP (so they match)
        resorted_raster_F = [resorted_raster_F; current_raster_F(gcamp_sort_ind,:)];
        resorted_raster_spks = [resorted_raster_spks; current_raster_spks(gcamp_sort_ind,:)];

        %Save cell order
        cellOrder_not_sorted = find(currentRows);
        cellOrder = [cellOrder; cellOrder_not_sorted(gcamp_sort_ind)];
        CellsInOrder.([Cells{c}]).([A{i}]) = cellnumbers_A ;
        
        end
        if c == 1
        VIP_raster_F = resorted_raster_F;
        VIP_raster_spks = resorted_raster_spks;
        VIP_average_F = average_F;
        VIP_average_spks = average_spks;
        VIP_cell_order = cellOrder;
    elseif c == 2
        NDNF_raster_F = resorted_raster_F;
        NDNF_raster_spks = resorted_raster_spks;
        NDNF_average_F = average_F;
        NDNF_average_spks = average_spks;
        NDNF_cell_order = cellOrder;
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
suptitle(ExtractedData.DataType)

subplot(6,2,1); hold all
plot(VIP_average_F', 'LineWidth', 2)
legend(A)
xlabel('Frames')
ylabel('DF/F')
title('VIP DF/F')

subplot(6,2,2); hold all
plot(VIP_average_spks', 'LineWidth', 2)
legend(A)
xlabel('Frames')
ylabel('Spikes')
title('VIP Spikes')

subplot(6,2,[3,5])
imagesc(VIP_raster_F(:,1:end-1))
ylabel('Cells')
%xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','DF/F');
caxis([-0.5, 1]);

subplot(6,2,[4,6])
imagesc(VIP_raster_spks(:,1:end-1))
ylabel('Cells')
%xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','Spikes');
caxis([0, 100]);

subplot(6,2,7)
plot(NDNF_average_F', 'LineWidth', 2)
legend(A)
xlabel('Frames')
ylabel('DF/F')
title('NDNF DF/F')

subplot(6,2,8)
plot(NDNF_average_spks', 'LineWidth', 2)
legend(A)
xlabel('Frames')
ylabel('Spikes')
title('NDNF Spikes')

subplot(6,2,[9,11])
imagesc(NDNF_raster_F(:,1:end-1))
ylabel('Cells')
xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','DF/F');
caxis([-0.5, 1]);

subplot(6,2,[10,12])
imagesc(NDNF_raster_spks(:,1:end-1))
ylabel('Cells')
xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','Spikes');
caxis([0, 100]);

%% Helper script to MANUALLY remove outliers (look up cell number from raster in VIP/NDNF cellOrder to find out which row to remove)

removeOutlier = 0;

outlier = 268;

if removeOutlier
    ExtractedData.AutoActivity(outlier,:) = [];
    ExtractedData.NominalData(outlier,:) = [];
    ExtractedData.NumericalData(outlier,:) = [];
    ExtractedData.Calcium_Raster(outlier,:) = [];
    ExtractedData.Spikes_Raster(outlier,:) = [];
    ExtractedData.RF_F(:,:,outlier) = [];
    ExtractedData.RF_isResponsive_F(:,:,outlier) = [];
    ExtractedData.RF_spks(:,:,outlier) = [];
    ExtractedData.RF_isResponsive_spks(:,:,outlier) = [];
end