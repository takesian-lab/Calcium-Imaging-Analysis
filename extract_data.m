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

columnHeaders = {'Group', 'Mouse ID', 'FOV', 'Data type', 'Block', 'Cell Number', 'Activity', 'Auto-Activity GCamP', 'Auto-Activity Spikes',...
    'GCaMP Peak Amplitude', 'GCaMP P1', 'GCaMP Peak Latency', 'GCaMP Peak Width',...
    'GCaMP Trough Amplitude', 'GCaMP T1', 'GCaMP Trough Latency', 'GCaMP Trough Width',...
    'Spike Peak Amplitude', 'Spike P1', 'Spike Peak Latency', 'Spike Peak Width',...
    'Spike Trough Amplitude', 'Spike T1', 'Spike Trough Latency', 'Spike Trough Width'};


dataType = 'water'; %To look at one stim type at a time. Leave empty to look at all
STDlevel = 2;
sort_active = 0;
plot_graphs = 0;
save_data = 1;

%% Load data
cellList_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
save_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
cellList_filename = 'Matching cells';
dataTypes = {'FM','RF','SAM','SAMfreq','NoiseITI','water','air'};%order of stim on cellList_filename

cd(cellList_path)
cellList = importfile(cellList_filename);
cellList(1,:) = []; %Remove header

%% Loop through blocks

cd(blocks_path)

%Look at only one stim type at a time
% dataTypes = [cellList{:,4}]';
% if ~isempty(dataType)
%     cellList = cellList(strcmp(dataTypes, dataType),:);
% end

data1 = {}; %Nominal data
autoActivity = cell(size(cellList,1),2); %Auto-determined activity (inhibited/sustained/activated)
data2 = nan(size(cellList,1),16); %Numerical data
raster_F = [];
raster_spks = [];

%Loop through all cells in each block
blocks = [cellList{:,5}]';
uniqueBlocks = unique(blocks);

count = 1;

for b = 1:length(uniqueBlocks)
    currentBlock = uniqueBlocks{b};
    block_cellList = cellList(strcmp(blocks, currentBlock),:);
    load(currentBlock)
    baseline_length = block.setup.constant.baseline_length; %seconds
    framerate = block.setup.framerate;
    nBaselineFrames = baseline_length*framerate; %frames
    
    data1 = [data1; block_cellList(:,1:7)];
    stim_v1 = block.parameters.variable1';
    stim_v2 = block.parameters.variable2';
    
    
    for c = 1:size(block_cellList,1)
        cellNumber = block_cellList{c,6};
        cellIndex = find(block.cell_number == cellNumber);
        if isempty(cellIndex)
            error(['Cell number ' num2str(cellNumber) ' not found.'])
        end
        
        %Cell-specific response to sound
        F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:));
        F7_baseline = F7(:,1:nBaselineFrames); %baseline for each trial
        F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
        spks = squeeze(block.aligned_stim.spks_stim(cellIndex,:,:));
        
        %Eliminate 0dB trials -> If they exist we could potentially use
        %them to confirm above-baseline responses
        %if desired, remove active trials here too
        try
            if dataType == 'FM' | 'RF' | 'NoiseITI';
                if sort_active==1
                    r=find(stim_v2 == 0);
                    rr= find(block.active_trials==1);
                    ru = union(r,rr);
                    F7_df_f(ru,:) = [];
                    spks(ru,:) = [];
                else
                    F7_df_f(stim_v2 == 0,:) = [];
                    spks(stim_v2 == 0,:) = [];
                end
            end
                catch
                    if sort_active==1
                        rr= find(block.active_trials==1);
                        F7_df_f(rr,:) = [];
                        spks(rr,:) = [];
                    end
            end
       
    
        
        
        
        %Get averaged & smoothed response
        avg_F7_df_f = smooth(nanmean(F7_df_f),3)';
        avg_spks = smooth(nanmean(spks),3)';
        
        %Store raster
        %If trial isn't the same length as raster, make them equal size by
        %adding nans at the end
        if count ~= 1
            if length(avg_F7_df_f) ~= size(raster_F,2) %assuming F and spks will be same dim
                dim = max([length(avg_F7_df_f),size(raster_F,2)]);
                while length(avg_F7_df_f) < dim
                    avg_F7_df_f = [avg_F7_df_f, nan];
                end
                while length(avg_spks) < dim
                    avg_spks = [avg_spks, nan];
                end
                while size(raster_F,2) < dim
                    raster_F = [raster_F, nan(size(raster_F,1),1)];
                end
                while size(raster_spks,2) < dim
                    raster_spks = [raster_spks, nan(size(raster_spks,1),1)];
                end
            end
        end
        
        raster_F = [raster_F; avg_F7_df_f];
        raster_spks = [raster_spks; avg_spks];
        
        %Compute latencies, width, and amplitude
        if plot_graphs == 1; figure; hold on; end
        for i = 1:2
            if i == 1
                y = avg_F7_df_f;
                units = 'DF/f';
            elseif i == 2
                y = avg_spks;
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
                [p2_latency] = find(response(1, peak_latency:end) <= peak_threshold, 1, 'first');
                p1 = response(p1_latency);
                p2 = response(peak_latency + p2_latency - 1);
                
                %Adjust for baseline
                peak_latency = peak_latency + nBaselineFrames;
                p1_latency = p1_latency + nBaselineFrames;
                p2_latency = p2_latency + peak_latency - 1;
                
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
                [t2_latency] = find(response(1, trough_latency:end) >= trough_threshold, 1, 'first');
                t1 = response(t1_latency);
                t2 = response(trough_latency + t2_latency - 1);
                
                %Adjust for baseline
                trough_latency = trough_latency + nBaselineFrames;
                t1_latency = t1_latency + nBaselineFrames;
                t2_latency = t2_latency + trough_latency - 1;
                
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
            end
            
            %Store
            if ~isempty(trough);            trough_data(1) = trough;            else;   trough = nan;           end
            if ~isempty(t1_latency);        trough_data(2) = t1_latency;        else;   t1_latency = nan;       end
            if ~isempty(trough_latency);    trough_data(3) = trough_latency;    else;   trough_latency = nan;   end
            if ~isempty(trough_width);      trough_data(4) = trough_width;      else;   trough_width = nan;     end
            
            %Auto-determined activity (inhibited/sustained/activated)
            if isnan(peak) && isnan(trough)
                activity = 'none';
            elseif isnan(peak) && ~isnan(trough)
                activity = 'inhibited';
            elseif ~isnan(peak) && isnan(trough)
                if peak_latency > 40 || isempty(p2_latency)
                    activity = 'sustained';
                else
                    activity = 'activated';
                end
            elseif ~isnan(peak) && ~isnan(trough)
                if trough_latency < peak_latency
                    activity = 'inhibited';
                else
                    if peak_latency > 40 || isempty(p2_latency)
                        activity = 'sustained';
                    else
                        activity = 'activated';
                    end
                end
            else
                activity = 'undetermined';
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
                title(activity)
                if i == 2
                    suptitle(strcat(block.setup.block_supname, ' Cell ', num2str(cellNumber)))
                end
            end
            
            %Convert to seconds and store data
            peak_data = peak_data./framerate;
            trough_data = trough_data./framerate;
            
            if i == 1
                data2(count,1:8) = [peak_data, trough_data]; %GCaMP
            elseif i == 2
                data2(count,9:16) = [peak_data, trough_data]; %Spikes
            end
            
        end
        count = count + 1;
    end
end

ExtractedData = struct;
ExtractedData.DataType = dataType;
ExtractedData.STDlevel = STDlevel;
ExtractedData.Sort_Active = sort_active;
ExtractedData.ColumnHeaders = columnHeaders;
ExtractedData.AutoActivity = autoActivity;
ExtractedData.NominalData = data1;
ExtractedData.NumericalData = data2;
ExtractedData.Calcium_Raster = raster_F;
ExtractedData.Spikes_Raster = raster_spks;

%% Save
if save_data == 1
    cd(save_path)
    save(['extractedData_' dataType '.mat'], 'ExtractedData');
end

%% Plot sorted rasters

%Plot by activity type and peak/trough amplitude

Cells = {'VIP', 'NDNF'};
A = {'activated', 'sustained', 'inhibited'};

cellList = [ExtractedData.NominalData{:,1}]'; %Modify
activityList = [ExtractedData.AutoActivity(:,1)]; %Modify

for c = 1:2
    currentCells = strcmp(cellList,Cells{c});
    
    resorted_raster_F = [];
    resorted_raster_spks = [];
    average_F = [];
    average_spks = [];
    
    for i = 1:length(A)
        currentActivity = A{i};
        activeRows = strcmp(activityList, currentActivity);
        currentRows = and(currentCells, activeRows);

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
        
    end
    if c == 1
        VIP_raster_F = resorted_raster_F;
        VIP_raster_spks = resorted_raster_spks;
        VIP_average_F = average_F;
        VIP_average_spks = average_spks;
    elseif c == 2
        NDNF_raster_F = resorted_raster_F;
        NDNF_raster_spks = resorted_raster_spks;
        NDNF_average_F = average_F;
        NDNF_average_spks = average_spks;
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
imagesc(VIP_raster_F)
ylabel('Cells')
%xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','DF/F');
caxis([0, 1]);

subplot(6,2,[4,6])
imagesc(VIP_raster_spks)
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
imagesc(NDNF_raster_F)
ylabel('Cells')
xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','DF/F');
caxis([0, 1]);

subplot(6,2,[10,12])
imagesc(NDNF_raster_spks)
ylabel('Cells')
xlabel('Frames')
h = colorbar;
set(get(h,'label'),'string','Spikes');
caxis([0, 100]);