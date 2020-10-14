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
    
dataType = 'FM'; %To look at one stim type at a time. Leave empty to look at all
STDlevel = 2;

%% Load data
cellList_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
cellList_filename = 'ResponsiveCells';

cd(cellList_path)
cellList = importfile(cellList_filename);
cellList(1,:) = []; %Remove header

%% Loop through blocks

cd(blocks_path)

%Look at only one stim type at a time
dataTypes = [cellList{:,4}]';
if ~isempty(dataType)
    cellList = cellList(strcmp(dataTypes, dataType),:);
end

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
        
        %Cell-specific response to sound
        F7 = squeeze(block.aligned_stim.F7_stim(cellIndex,:,:));
        F7_baseline = F7(:,1:nBaselineFrames); %baseline for each trial
        F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
        spks = squeeze(block.aligned_stim.spks_stim(cellIndex,:,:));
        
        %Eliminate 0dB trials -> If they exist we could potentially use
        %them to confirm above-baseline responses
        F7_df_f(stim_v2 == 0,:) = [];
        spks(stim_v2 == 0,:) = [];
        
        %Get averaged & smoothed response
        avg_F7_df_f = smooth(mean(F7_df_f),3)';
        avg_spks = smooth(mean(spks),3)';
        
        %Store raster
        if count ~= 1
            if length(avg_F7_df_f) ~= size(raster_F,2)
                dim = max([length(avg_F7_df_f),size(raster_F,2)]);
                while length(avg_F7_df_f) < dim
                    avg_F7_df_f = [avg_F7_df_f, nan];
                end
                while length(avg_spks) < dim
                    avg_spks = [avg_spks, nan];
                end
                while size(raster_F,2) < dim
                    raster_f = [raster_F, nan(size(raster_F,1),1)];
                end
            end
        end
                
        raster_F = [raster_F; avg_F7_df_f];
        raster_spks = [raster_spks; avg_spks];

        %Compute latencies, width, and amplitude
        figure; hold on
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
            subplot(1,2,i); hold on
            plot(y)
            hline(mean(baseline), 'k')
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
ExtractedData.ColumnHeaders = columnHeaders;
ExtractedData.AutoActivity = autoActivity;
ExtractedData.NominalData = data1;
ExtractedData.NumericalData = data2;
ExtractedData.Calcium_Raster = raster_F;
ExtractedData.Spikes_Raster = raster_spks;

%% Save
%save('extractedData.mat', 'ExtractedData');