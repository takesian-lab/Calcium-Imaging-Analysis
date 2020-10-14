function visualize_cell_cgs(block, cellnum) 
% This function allows you to preview the stimulus-evoked responses from
% a single cell (neuron) or selection of cells from a block
%
% cellnum is a 1-D array of cell numbers matching the Suite2p GUI
% If the array is horizontal, the results from all cells will be averaged
% and plotted together. If the array is vertical, each cell will be plotted
% independently
%
% Error bars are SEM when N > 1 (showing inter-cell variability)
% and STD when N == 1 (showing inter-trial variability)
%
% Argument(s): 
%   block (struct)
%   cellnum (1-D array of cell numbers)
% 
% Returns:
% 
% 
% Notes:
%
%
% TODO: 
% Search 'TODO'

%% Magic numbers & Setup

SF = 0.5; %Shrinking factor for traces to appear more spread out (for visualization purposes)
z = 1; %Portion of recording to plot between 0 and 1 e.g. 0.5, 0.33, 1 (for visualization purposes)
ws = 0.15; %Multiplication factor for adding white space to the top of each graph (Above max value)
t1 = 20; %Time of response window onset in frames (for receptive field plots)
t2 = 40; %Time of response window offset in frames (for receptive field plots)

if ~isfield(block,'aligned_stim')
    error('No stim-aligned data to plot');
end    

if size(cellnum,1) > 1 && size(cellnum,2) > 1
    error('cellnum should be a 1-D array')
end

setup = block.setup;
stim_protocol = setup.stim_protocol;
code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'Widefield', 'SAM', 'SAM freq' , 'Behavior', 'Behavior', 'Random H20', 'Noiseburst ITI', 'Random Air', 'Spontaneous', 'Behavior Maryse'};
currentStim = code{stim_protocol};
disp(['Plotting figures for ' currentStim '...'])

%Raw activity
Sound_Time = block.Sound_Time;
all_cell_numbers = block.cell_number;    
F = block.F; %all the cell fluorescence data
Fneu = block.Fneu; %neuropil
F7 = F-setup.constant.neucoeff*Fneu; %neuropil corrected traces
timestamp = block.timestamp; %In seconds
Z = round(length(timestamp)*z);

%Stim-aligned activity
F7_stim = block.aligned_stim.F7_stim;
spks_stim = block.aligned_stim.spks_stim;
baseline_length = block.setup.constant.baseline_length; %seconds
framerate = block.setup.framerate;
nBaselineFrames = baseline_length*framerate; %frames
trial_duration_in_seconds = baseline_length + block.setup.constant.after_stim; %seconds
trial_duration_in_frames = size(F7_stim,3);
x_in_seconds = 0:0.5*(trial_duration_in_frames/trial_duration_in_seconds):trial_duration_in_frames;
x_label_in_seconds = 0:0.5:trial_duration_in_seconds;

%% Plot raw activity of cell(s) for duration of block
%Raw activity of each cell vs. time with locomoter activity beneath
%Each cell is a separate trace

figure('units','normalized','outerposition',[0 0 1 1])
count = 1; %for staggering plot lines

for c = 1:length(cellnum)
    current_cellnum = cellnum(c);

    subplot(3,4,1:8); hold on

    row_num = find(all_cell_numbers == current_cellnum);
    if isempty(row_num)
        error(['Cell ' num2str(current_cellnum) ' was not found']);
    end

    cell_trace = F7(row_num,:);%pull out the full trace for each cell

    mean_gCAMP = nanmean(cell_trace);% average green for each cell
    df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;%(total-mean)/mean
    A = smooth(df_f,10);

    plot(timestamp, A*SF + count,'LineWidth',1);
    count = count + 1;
end

suptitle(block.setup.block_supname)
title('DF/F')
xlim([0 timestamp(Z)])
xlabel('Time (s)')
set(gca, 'YTick', [1:1:count-1])
set(gca, 'YTickLabel', [cellnum(1:count-1)])
ylabel('Cell')

%Vertical lines for sound times
if Z < 15000 && isfield(block, 'Sound_Time') %Don't plot red lines if there is too much data, otherwise its messy
    %plot multicolored lines if less than 8 stim, else plot red lines
    if isfield(block.parameters, 'variable1')
            var1 = unique(block.parameters.variable1);
            variable1 = block.parameters.variable1;
        if length(variable1) > 1 && length(var1) < 8
            for i = 1:length(var1)
                colours = {'r', 'g', 'k', 'b', 'y', 'm', 'c'};
                vline(Sound_Time(variable1 == var1(i)), colours{i})
            end
        else
            vline(Sound_Time, 'r');
        end
    else
        vline(Sound_Time, 'r');
    end
end

%Plot locomotor activity
if ~ismissing(block.setup.Tosca_path)
    loco_time = block.loco_times;
    loco_speed = block.loco_activity;
        
    subplot(3,4,9:12); hold on %loco
    plot(loco_time, loco_speed);
    title('Locomotor activity')
    ylabel('Activity (cm/s)')
    xlim([0 timestamp(Z)])
    xlabel('Time (s)')
end
  
    
%% Plot graphs according to stim presentation

% Average response to all stim
for f = 1:size(cellnum,1) %Individual figures if cellnum is vertical
        current_cells = cellnum(f,:);
        
        row_nums = nan(length(current_cells),1);
        for c = 1:length(current_cells)
            row_nums(c) = find(all_cell_numbers == current_cells(c));
        end
    % individual cells
        if length(current_cells) == 1
            F7_cell = squeeze(F7_stim(row_nums,:,:));
            F7_baseline = F7_cell(:,1:nBaselineFrames); %baseline for each trial
            F7_df_f = (F7_cell-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
            ebar = std(F7_df_f,1)./sqrt(size(F7_df_f,1));
            spks_cell = squeeze(spks_stim(row_nums,:,:));
            ebar_spks = std(spks_cell,1)./sqrt(size(spks_cell,1));
    % average the cells 
        elseif length(current_cells) > 1
            F7_cells = F7_stim(row_nums,:,:);
            F7_baselines = F7_cells(:,:,1:nBaselineFrames);
            F7_df_fs = (F7_cells-nanmean(F7_baselines,3))./nanmean(F7_baselines,3);
            F7_df_f = squeeze(nanmean(F7_df_fs,1));
            ebar = std(F7_df_f,1)/sqrt(length(current_cells)-1);
            spks_cell = squeeze(nanmean(spks_stim(row_nums,:,:),1));
            ebar_spks = std(spks_cell,1)/sqrt(length(current_cells)-1);
        end
        
        figure %hold all
        
        %DF_F average
        subplot(2,2,1)
        total_mean = nanmean(F7_df_f);
%         total_sem = F7_df_f
        shadedErrorBar(1:length(total_mean),smooth(total_mean,3),smooth(ebar,3));
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 trial_duration_in_frames])
        ylabel('DF/F')
        vline(nBaselineFrames,'k') %stim onset
        
        %Spike average
        subplot(2,2,2)
        total_mean = nanmean(spks_cell);
        area(total_mean, 'FaceColor', [0.85 0.85 0.85]);
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 trial_duration_in_frames])
        ylabel('Deconvolved spikes')
        vline(nBaselineFrames,'k') %stim onset
        
        %DF_F raster
        subplot(2,2,3)
        imagesc(F7_df_f)
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 size(F7_df_f,2)])
        ylabel('Trials')
        vline(nBaselineFrames,'k') %stim onset
        
        %Spike raster
        subplot(2,2,4)
        imagesc(spks_cell)
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 size(spks_cell,2)])
        ylabel('Trials')
        vline(nBaselineFrames,'k') %stim onset
        
        if length(current_cells) == 1
            suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
        else
            suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
        end
        
end

%% Average response separated by stim type - air and H20

plotAirOrH20 = 0;

if stim_protocol == 9 %Random H20
    V1 = block.parameters.variable1;
    stim_names = {'H20', 'No H20'};
    plotAirOrH20 = 1;
elseif stim_protocol == 11 %Random Air
    V1 = block.parameters.variable1;
    stim_names = {'Air', 'No Air'};
    plotAirOrH20 = 1;
elseif stim_protocol == 10 %Noiseburst ITI
    V1 = block.parameters.variable1;
    stim_names = {'70dB', '0dB'};
    if length(unique(V1)) == 2 %Some versions of noisburst ITI had multiple dB levels
%         warning('Version of noiseburst ITI with multiple dB levels. Sound vs. No Sound not plotted.')
        plotAirOrH20 = 1;
    end
else %Use same format to plot sound blocks where blank/sham stim were included
    V1 = block.parameters.variable1;
    V2 = block.parameters.variable2;
    
    if any(isnan(V1)) || any(isnan(V2))
        V1 = ~isnan(V1);
        V2 = ~isnan(V2);
        if ~isequal(V1,V2)
            error('Stil need to account for case where V1 and V2 have different Nans')
        end
        stim_names = {'Sound', 'No sound'};
        plotAirOrH20 = 1;
    end
end

if plotAirOrH20
    
    for f = 1:size(cellnum,1) %Individual figures if cellnum is vertical
        current_cells = cellnum(f,:);
        
        row_nums = nan(length(current_cells),1);
        for c = 1:length(current_cells)
            row_nums(c) = find(all_cell_numbers == current_cells(c));
        end
        
        figure; hold all
        
        stimValues = fliplr(unique(V1));  %first is stim, second is sham
        
        for v = 1:length(stimValues)
            stim_name = stim_names(v);
            if length(current_cells) == 1
                F7_cell = squeeze(F7_stim(row_nums,V1 == stimValues(v),:));
                F7_baseline = F7_cell(:,1:nBaselineFrames); %baseline for each trial
                F7_df_f = (F7_cell-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
                ebar = std(F7_df_f,1)./sqrt(size(F7_df_f,1));
                spks_cell = squeeze(spks_stim(row_nums,V1 == stimValues(v),:));
                ebar_spks = std(spks_cell,1);
            elseif length(current_cells) > 1
                F7_cells = F7_stim(row_nums,V1 == stimValues(v),:);
                F7_baselines = F7_cells(:,:,1:nBaselineFrames);
                F7_df_fs = (F7_cells-nanmean(F7_baselines,3))./nanmean(F7_baselines,3);
                F7_df_f = squeeze(nanmean(F7_df_fs,1));
                ebar = std(F7_df_f,1)/sqrt(length(current_cells)-1);
                spks_cell = squeeze(nanmean(spks_stim(row_nums,V1 == stimValues(v),:),1));
                ebar_spks = std(spks_cell,1)/sqrt(length(current_cells)-1);
            end   
            
            %Shift figure position in for loop
            if v == 2; q = 4; else; q = 0; end
        
            %DF_F average
            subplot(4,2,1+q)
            total_mean = nanmean(F7_df_f);
            shadedErrorBar(1:length(total_mean),smooth(total_mean,3),smooth(ebar,3));
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 trial_duration_in_frames])
            ylabel('DF/F')
            vline(nBaselineFrames,'k') %stim onset
            title(stim_name)

            %Spike average
            subplot(4,2,2+q)
            total_mean = nanmean(spks_cell);
            area(total_mean, 'FaceColor', [0.85 0.85 0.85]);
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 trial_duration_in_frames])
            ylabel('Deconvolved spikes')
            vline(nBaselineFrames,'k') %stim onset
            title(stim_name)

            %DF_F raster
            subplot(4,2,3+q)
            imagesc(F7_df_f)
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 trial_duration_in_frames])
            vline(nBaselineFrames,'k') %stim onset
            ylabel('Trials')

            %Spike raster
            subplot(4,2,4+q)
            imagesc(spks_cell)
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 trial_duration_in_frames])
            ylabel('Trials')
            vline(nBaselineFrames,'k') %stim onset
        end

        if length(current_cells) == 1
            suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
        else
            suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
        end            

    end
end

%% Average response separated by stim type - Receptive Field

plotReceptiveField = 0;

if stim_protocol == 2 %Receptive Field
    stim_units = {'kHz', 'dB'};
    V1_label = 'Frequency (kHz)';
    V2_label = 'Intensity (dB)';
    plotReceptiveField = 1;
elseif stim_protocol == 3 %FM sweep
    stim_units = {'', 'dB'};
    V1_label = 'Frequency (kHz)';
    V2_label = 'Speed';
    plotReceptiveField = 1;
elseif stim_protocol == 5 %SAM
    stim_units = {'Hz', ''};
    V1_label = 'Rate (Hz)';
    V2_label = 'Modulation Depth';
    plotReceptiveField = 1; 
elseif stim_protocol == 6 %SAM freq
    stim_units = {'kHz', ''};
    V1_label = 'Frequency (kHz)';
    V2_label = 'Modulation Depth';
    plotReceptiveField = 1; 
end

if plotReceptiveField

    V1 = block.parameters.variable1;
    V2 = block.parameters.variable2;

    V1_stim = unique(V1);
    V2_stim = fliplr(unique(V2));
    
    if any(isnan(V1_stim))
        V1_stim(isnan(V1_stim)) = []; % remove all nans
        %V1_stim(end+1) = NaN; % add the unique one.
    end
    
    if any(isnan(V2_stim))
        V2_stim(isnan(V2_stim)) = []; % remove all nans
        %V2_stim(end+1) = NaN; % add the unique one.
    end
    
    for f = 1:size(cellnum,1) %Individual figures if cellnum is vertical
        current_cells = cellnum(f,:);

        row_nums = nan(length(current_cells),1);
        for c = 1:length(current_cells)
            row_nums(c) = find(all_cell_numbers == current_cells(c));
        end
        
        %Preallocate
        F7_df_f_mat = nan(length(V1_stim)*length(V2_stim),trial_duration_in_frames);
        ebar_mat = nan(length(V1_stim)*length(V2_stim),trial_duration_in_frames);
        spks_mat = nan(length(V1_stim)*length(V2_stim),trial_duration_in_frames);
        ebar_spks_mat = nan(length(V1_stim)*length(V2_stim),trial_duration_in_frames);
        %Record stim as well for sanity check
        V1_list = nan(length(V1_stim)*length(V2_stim),1);
        V2_list = nan(length(V1_stim)*length(V2_stim),1);
        
        count = 1;
        for v = 1:length(V2_stim) %Intensities
            for vv = 1:length(V1_stim) %Frequencies
                
                % this is where I would put active/inactive trials...
                stim_rows = intersect(find(V1 == V1_stim(vv)), find(V2 == V2_stim(v)));
                V1_list(count) = V1_stim(vv);
                V2_list(count) = V2_stim(v);
                
                if length(current_cells) == 1
                    F7_cell = squeeze(F7_stim(row_nums,stim_rows,:));
                    F7_baseline = F7_cell(:,1:nBaselineFrames); %baseline for each trial
                    F7_df_f = (F7_cell-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
                    ebar = std(F7_df_f,1)./sqrt(size(F7_df_f,1));
                    spks_cell = squeeze(spks_stim(row_nums,stim_rows,:));
                    ebar_spks = std(spks_cell,1);
                    F7_df_f_mat(count,:) = nanmean(F7_df_f);
                    ebar_mat(count,:) = ebar;
                    spks_mat(count,:) = nanmean(spks_cell);
                    ebar_spks_mat(count,:) = ebar_spks;
                elseif length(current_cells) > 1
                    F7_cells = F7_stim(row_nums,stim_rows,:);
                    F7_baselines = F7_cells(:,:,1:nBaselineFrames);
                    F7_df_fs = (F7_cells-nanmean(F7_baselines,3))./nanmean(F7_baselines,3);
                    F7_df_f = squeeze(nanmean(F7_df_fs,1));
                    ebar = std(F7_df_f,1)/sqrt(length(current_cells)-1);
                    spks_cell = squeeze(nanmean(spks_stim(row_nums,stim_rows,:),1));
                    ebar_spks = std(spks_cell,1)/sqrt(length(current_cells)-1);
                    F7_df_f_mat(count,:) = nanmean(F7_df_f);
                    ebar_mat(count,:) = ebar;
                    spks_mat(count,:) = nanmean(spks_cell);
                    ebar_spks_mat(count,:) = ebar_spks;
                end
                
                count = count + 1;
            end
        end
    
        %Plot DF/F and Spikes for V1 and V2 separately
        for v = 1:2
            if v == 1
                V_stim = V1_stim;
                V_list = V1_list;
            else
                V_stim = V2_stim;
                V_list = V2_list;
            end
            
            if length(V_stim) == 1
                continue; %Don't plot if there's only 1 stimulus value (e.g. all stim at same dB level)
            end
            
            %Find max values for y limits
            max_df_f = 0;
            max_spks = 0;
            for p = 1:length(V_stim)
                mat_rows = V_list == V_stim(p);
                if sum(mat_rows) == 1
                    df_f_std = F7_df_f_mat(mat_rows,:) + ebar_mat(mat_rows,:);
                    spks_mean = spks_mat(mat_rows,:);
                else
                    df_f_std = nanmean(F7_df_f_mat(mat_rows,:)) + nanstd(F7_df_f_mat(mat_rows,:));
                    spks_mean = nanmean(spks_mat(mat_rows,:));
                end
                if max(df_f_std) > max_df_f
                    max_df_f = max(df_f_std);
                end
                if max(spks_mean) > max_spks
                    max_spks = max(spks_mean);
                end
            end
            
            figure; hold on
            for p = 1:length(V_stim)
                mat_rows = V_list == V_stim(p);
                
                %DF_F average
                subplot(length(V_stim),2,p*2-1)
                if sum(mat_rows) == 1
                    total_mean = F7_df_f_mat(mat_rows,:);
                    ebar = ebar_mat(mat_rows,:);
                else
                    total_mean = nanmean(F7_df_f_mat(mat_rows,:));
                    ebar = std(F7_df_f_mat,1); %recompute error bar for this part - not sure the best way to do this
                end
                shadedErrorBar(1:length(total_mean),smooth(total_mean,3),smooth(ebar,3));
                set(gca, 'XTick', x_in_seconds)
                set(gca, 'XTickLabel', x_label_in_seconds)
                xlim([0 trial_duration_in_frames])
                ylabel([num2str(V_stim(p)) ' ' stim_units{v}])
                ylim([0 max_df_f]);
                vline(nBaselineFrames,'k') %stim onset
                if p == 1
                    title('DF/F')
                elseif p == length(V_stim)
                    xlabel('Time (s)')
                end
               
                %Spike average
                subplot(length(V_stim),2,p*2)
                if sum(mat_rows) == 1
                    total_mean = spks_mat(mat_rows,:);
                else
                    total_mean = nanmean(spks_mat(mat_rows,:));
                end
                area(total_mean, 'FaceColor', [0.85 0.85 0.85]);
                set(gca, 'XTick', x_in_seconds)
                set(gca, 'XTickLabel', x_label_in_seconds)
                xlim([0 trial_duration_in_frames])
                ylabel([num2str(V_stim(p)) ' ' stim_units{v}])
                ylim([0 max_spks+(ws*max_spks)]);
                vline(nBaselineFrames,'k') %stim onset
                if p == 1
                    title('Deconvolved spikes')
                elseif p == length(V_stim)
                    xlabel('Time (s)')
                end
                
            end
            
            if length(current_cells) == 1
                suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
            else
                suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
            end 
        end
                
        %Plot V1 x V2 grid using DF/F or Spikes
        if length(V1_stim) > 1 && length(V2_stim) > 1 %Don't plot if there isn't a grid

            nPlots = length(V1_stim)*length(V2_stim);
            
            %Find max values for y limits
            max_df_f = max(max(F7_df_f_mat,[], 1)) + max(max(ebar_mat,[], 1));
            max_spks = max(max(spks_mat,[], 1));

            %DF_F average
            figure; hold on
            for p = 1:nPlots
                subplot(length(V2_stim),length(V1_stim),p)
                plot(F7_df_f_mat(p,:), 'Linewidth', 2);
                %shadedErrorBar(1:length(F7_df_f_mat(p,:)),F7_df_f_mat(p,:),ebar_mat(p,:));
                set(gca, 'XTick', x_in_seconds)
                set(gca, 'XTickLabel', x_label_in_seconds)
                xlim([0 trial_duration_in_frames])
                ylim([0 max_df_f]);
                vline(nBaselineFrames,'k') %stim onset
                if p <= length(V1_stim)%Top row
                    title([num2str(V1_list(p)) ' ' stim_units{v}])
                end
                if p > nPlots - length(V1_stim) %Bottom row
                    xlabel('Time (s)')
                end
                if ismember(p,[1:length(V1_stim):nPlots]) %First column
                    ylabel([num2str(V2_list(p)) ' ' stim_units{v}])
                end
            end
            
            if length(current_cells) == 1
                suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
            else
                suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
            end 
            
            %Spike average
            figure; hold on
            for p = 1:nPlots
                subplot(length(V2_stim),length(V1_stim),p)
                area(spks_mat(p,:), 'FaceColor',[0, 0.4470, 0.7410]);
                set(gca, 'XTick', x_in_seconds)
                set(gca, 'XTickLabel', x_label_in_seconds)
                xlim([0 trial_duration_in_frames])
                ylim([0 max_spks+(ws*max_spks)]);
                vline(nBaselineFrames,'k') %stim onset
                if p <= length(V1_stim)%Top row
                    title([num2str(V1_list(p)) ' ' stim_units{1}])
                end
                if p > nPlots - length(V1_stim) %Bottom row
                    xlabel('Time (s)')
                end
                if ismember(p,[1:length(V1_stim):nPlots]) %First column
                    ylabel([num2str(V2_list(p)) ' ' stim_units{2}])
                end
            end

            if length(current_cells) == 1
                suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
            else
                suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
            end 
        end
        
                        
        %Plot Receptive field grid using DF/F peak between t1 and t2
        if length(V1_stim) > 1 && length(V2_stim) > 1 %Don't plot if there isn't a grid
        
            RF = nan(length(V2_stim),length(V1_stim));
            count = 1;
            for v = 1:length(V2_stim)
                for vv = 1:length(V1_stim)
                    RF(v,vv) = nanmean(F7_df_f_mat(count,t1:t2));
                    count = count + 1;
                end
            end
    
            figure;
            imagesc(RF)
            set(gca,'XTickLabel',V1_stim)
            set(gca, 'YTickLabel', V2_stim)
            xlabel(V1_label)
            ylabel(V2_label)
            h = colorbar;
            set(get(h,'title'),'string','DF/F');

            if length(current_cells) == 1
                suptitle([block.setup.block_supname...
                    strcat(' Cell #', num2str(cellnum(f)))])
            else
                suptitle([block.setup.block_supname...
                    strcat('Average of ', num2str(length(cellnum)), ' cells')])
            end            
        end   
    end
end

end %function end

