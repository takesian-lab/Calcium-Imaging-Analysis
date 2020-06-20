function visualize_cell(block, cellnum)

% DOCUMENTATION IN PROGRESS
% 
% This function allows you to preview the data from a single cell (neuron) or
% selection of cells from a block
%
% cellnum is a 1-D array of cell numbers or labels, matching the Suite2p GUI
% If the array is horizontal, the results from all cells will be averaged
% and plotted together. If the array is vertical, each cell will be plotted
% independently
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

if ~isfield(block,'aligned_stim')
    error('No stim-aligned data to plot');
end    

if size(cellnum,1) > 1 && size(cellnum,2) > 1
    error('cellnum should be a 1-D array')
end

setup = block.setup;
stim_protocol = setup.stim_protocol;
code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'Widefield', 'SAM', 'SAM freq' , 'Behavior', 'Behavior', 'Random H20', 'Noiseburst ITI', 'Random Air'};
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
trial_duration_in_frames = trial_duration_in_seconds*framerate;
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

    mean_gCAMP = mean(cell_trace);% average green for each cell
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
    if isfield(block, 'locomotion_data')
        loco_data = block.locomotion_data;
    else
        loco_data = block.loco_data;
    end
    
    subplot(3,4,9:12); hold on %loco
    plot(loco_data(:,1), loco_data(:,3));
    title('Locomotor activity')
    ylabel('Activity (cm/s)')
    xlim([0 timestamp(Z)])
    xlabel('Time (s)')
end
  
    
%% Plot graphs according to stim presentation

% Plot 1 - average response to all stim
for f = 1:size(cellnum,1) %Individual figures if cellnum is vertical
        current_cells = cellnum(f,:);
        
        row_nums = nan(length(current_cells),1);
        for c = 1:length(current_cells)
            row_nums(c) = find(all_cell_numbers == current_cells(c));
        end
    
        if length(current_cells) == 1
            F7_cell = squeeze(F7_stim(row_nums,:,:));
            F7_baseline = F7_cell(:,1:nBaselineFrames); %baseline for each trial
            F7_df_f = (F7_cell-mean(F7_baseline,2))./mean(F7_baseline,2); %(total-mean)/mean
            ebar = std(F7_df_f,1);
            spks_cell = squeeze(spks_stim(row_nums,:,:));
            ebar_spks = std(spks_cell,1);
        elseif length(current_cells) > 1
            F7_cells = F7_stim(row_nums,:,:);
            F7_baselines = F7_cells(:,:,1:nBaselineFrames);
            F7_df_fs = (F7_cells-mean(F7_baselines,3))./mean(F7_baselines,3);
            F7_df_f = squeeze(mean(F7_df_fs,1));
            ebar = std(F7_df_f,1)/sqrt(length(current_cells)-1);
            spks_cell = squeeze(mean(spks_stim(row_nums,:,:),1));
            ebar_spks = std(spks_cell,1)/sqrt(length(current_cells)-1);
        end
        
        figure %hold all
        
        %DF_F average
        subplot(2,2,1)
        total_mean = mean(F7_df_f);
        shadedErrorBar(1:length(total_mean),total_mean,ebar);
        vline(nBaselineFrames,'k') %stim onset
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 trial_duration_in_frames])
        ylabel('DF/F')       
        
        %Spike average
        subplot(2,2,2)
        total_mean = mean(spks_cell);
        area(total_mean, 'FaceColor', [0.85 0.85 0.85]);
        vline(nBaselineFrames,'k') %stim onset
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 trial_duration_in_frames])
        ylabel('Deconvolved spikes')
        
        %DF_F raster
        subplot(2,2,3)
        imagesc(F7_df_f)
        vline(nBaselineFrames,'k') %stim onset
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 size(F7_df_f,2)])
        ylabel('Trials')
        
        %Spike raster
        subplot(2,2,4)
        imagesc(spks_cell)
        vline(nBaselineFrames,'k') %stim onset
        set(gca, 'XTick', x_in_seconds)
        set(gca, 'XTickLabel', x_label_in_seconds)
        xlabel('Time (s)')
        xlim([0 size(spks_cell,2)])
        ylabel('Trials')
        
        if length(current_cells) == 1
            suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
        else
            suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
        end
        
end

%% Plot 2 - average response separated by stim type - air and H20

plotAirOrH20 = 0;

if stim_protocol == 9 %Random H20
    stim_names = {'H20', 'No H20'};
    plotAirOrH20 = 1;
elseif stim_protocol == 11 %Random Air
    stim_names = {'Air', 'No Air'};
    plotAirOrH20 = 1;
elseif stim_protocol == 10 %Noiseburst ITI
    stim_names = {'70dB', '0dB'};
    plotAirOrH20 = 1;
end

if plotAirOrH20
    
    V1 = block.parameters.variable1;

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
                F7_df_f = (F7_cell-mean(F7_baseline,2))./mean(F7_baseline,2); %(total-mean)/mean
                ebar = std(F7_df_f,1);
                spks_cell = squeeze(spks_stim(row_nums,V1 == stimValues(v),:));
                ebar_spks = std(spks_cell,1);
            elseif length(current_cells) > 1
                F7_cells = F7_stim(row_nums,V1 == stimValues(v),:);
                F7_baselines = F7_cells(:,:,1:nBaselineFrames);
                F7_df_fs = (F7_cells-mean(F7_baselines,3))./mean(F7_baselines,3);
                F7_df_f = squeeze(mean(F7_df_fs,1));
                ebar = std(F7_df_f,1)/sqrt(length(current_cells)-1);
                spks_cell = squeeze(mean(spks_stim(row_nums,V1 == stimValues(v),:),1));
                ebar_spks = std(spks_cell,1)/sqrt(length(current_cells)-1);
            end   
            
            %Shift figure position in for loop
            if v == 2
                q = 4;
            else
                q = 0;
            end
        
            %DF_F average
            subplot(4,2,1+q)
            total_mean = mean(F7_df_f);
            shadedErrorBar(1:length(total_mean),total_mean,ebar);
            vline(nBaselineFrames,'k') %stim onset
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 trial_duration_in_frames])
            ylabel('DF/F')
            title(stim_name)

            %Spike average
            subplot(4,2,2+q)
            total_mean = mean(spks_cell);
            area(total_mean, 'FaceColor', [0.85 0.85 0.85]);
            vline(nBaselineFrames,'k') %stim onset
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 trial_duration_in_frames])
            ylabel('Deconvolved spikes')
            title(stim_name)

            %DF_F raster
            subplot(4,2,3+q)
            imagesc(F7_df_f)
            vline(nBaselineFrames,'k') %stim onset
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 size(F7_df_f,2)])
            ylabel('Trials')

            %Spike raster
            subplot(4,2,4+q)
            imagesc(spks_cell)
            vline(nBaselineFrames,'k') %stim onset
            set(gca, 'XTick', x_in_seconds)
            set(gca, 'XTickLabel', x_label_in_seconds)
            xlabel('Time (s)')
            xlim([0 size(spks_cell,2)])
            ylabel('Trials')
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

%% Receptive field

V1 = block.parameters.variable1;
V2 = block.parameters.variable2;

freqs = unique(V1);
ints = fliplr(unique(V2));

t1 = 20;
t2 = 40;

for f = 1:length(cellnum)
    
    RF = nan(length(ints),length(freqs));
    
        current_cellnum = cellnum(f);
        row_num = find(all_cell_numbers == current_cellnum);
        F7_cell = squeeze(F7_stim(row_num,:,:));
   
        figure; hold on
        for i = 1:length(ints)
        subplot(length(ints),1,i)
        int_mean = mean(F7_cell(V2 == ints(i),:));
        plot(int_mean);
        title([num2str(ints(i)*100) '%'])
        ylabel('DF/F')
        end
        
        figure; hold on
        q = 1;
        for i = 1:length(ints)
            for k = 1:length(freqs)
                
                rows = intersect(find(V1 == freqs(k)), find(V2 == ints(i)));
        
                subplot(length(ints),length(freqs),q)
                RF_mean = mean(F7_cell(rows,:));
                RF(i,k) = mean(RF_mean(1,t1:t2));
                plot(RF_mean, 'LineWidth', 2);
                ylim([0 3000])
                %ylabel('DF/F')
                %xlabel([num2str(freqs(k))])
                
                q = q+1;
            end
        end
        
        figure;
        imagesc(RF)
        xlabel('Frequency (kHz)')
        set(gca,'XTickLabel',freqs)
        ylabel('Modulation %')
        set(gca, 'YTickLabel', ints)
        h = colorbar;
        set(get(h,'title'),'string','DF/F');
end


end

