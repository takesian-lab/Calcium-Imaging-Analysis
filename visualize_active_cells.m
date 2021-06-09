function visualize_active_cells(block) 
% This function allows you to preview the stimulus-evoked responses from
% neurons deemed to be significantly active within a block
%
%
% Error bars are SEM when N > 1 (showing inter-cell variability)
% and STD when N == 1 (showing inter-trial variability)
%
% Argument(s): 
%   block (struct)
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

%Determine if cell is active
STDlevel = 2;
AUC_F_level = 5;
AUC_S_level = 10;
sort_active = 1;  % 0= dont perform, 1= non-locomotor trials, 2= locomotor trials
analyze_by_stim_condition = 1; %determine if cell is active based on individual stim conditions

disp('===========PARAMETERS===========')
disp(['STD level = ' num2str(STDlevel)])
disp(['AUC level for df/f = ' num2str(AUC_F_level)])
disp(['AUC level for spks = ' num2str(AUC_S_level)])
disp(['Sort active = ' num2str(sort_active)])
disp(['Analyze by stim condition = ' num2str(analyze_by_stim_condition)])
        
%For plotting
SF = 0.5; %Shrinking factor for traces to appear more spread out (for visualization purposes)
z = 1; %Portion of recording to plot between 0 and 1 e.g. 0.5, 0.33, 1 (for visualization purposes)
ws = 0.15; %Multiplication factor for adding white space to the top of each graph (Above max value)
t1 = 20; %Time of response window onset in frames (for receptive field plots)
t2 = 40; %Time of response window offset in frames (for receptive field plots)

if ~isfield(block,'aligned_stim')
    error('No stim-aligned data to plot');
end    

setup = block.setup;
stim_protocol = setup.stim_protocol;
code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'Widefield', 'SAM', 'SAM freq' , 'Behavior', 'Behavior', 'Random H20', 'Noiseburst ITI', 'Random Air', 'Spontaneous', 'Behavior Maryse'};
currentStim = code{stim_protocol};
disp(['Plotting figures for ' currentStim '...'])

%Raw activity
Sound_Time = block.Sound_Time;
all_cell_numbers = block.cell_number;    
F7_full = block.F - (setup.constant.neucoeff*block.Fneu); %neuropil corrected traces
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

%% Find significantly responsive cells

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
switch currentStim
    case {'FM sweep','Receptive Field'}
        stim_v0 = stim_v2;

    case {'Noiseburst ITI', 'Random H20', 'Random Air'}
        %Regular noise is not accounted for yet
        stim_v0 = stim_v1;
        stim_v2 = zeros(size(stim_v1));

    case {'SAM', 'SAM freq'} %NAN trials instead of zeros           
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
    
responsiveCells_F = zeros(size(block.cell_number));
responsiveCells_S = zeros(size(block.cell_number));

for c = 1:size(block.cell_number,1)

    %when we remove inf below stim might change so refresh it with original stim list
    stim_v1 = store_stim_v1;
    stim_v2 = store_stim_v2;
    stimTrials = store_stimTrials;
    blankTrials = store_blankTrials;

    F7 = squeeze(block.aligned_stim.F7_stim(c,:,:));
    F7_baseline = F7(:,1:nBaselineFrames); %baseline for each trial
    F7_df_f = (F7-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %compute df/f: (total-mean)/mean
    spks = squeeze(block.aligned_stim.spks_stim(c,:,:));

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
                [F_active, ~, ~, ~] = checkIfActive_v2(F_temp_smoothed, nBaselineFrames, STDlevel, AUC_F_level, 0);
                [S_active, ~, ~, ~] = checkIfActive_v2(S_temp_smoothed, nBaselineFrames, STDlevel, AUC_F_level, 0);

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

    %CHECK IF RESPONSIVE
    [responsiveCells_F(c), ~, ~, ~] = checkIfActive_v2(avg_F, nBaselineFrames, STDlevel, AUC_F_level, 0);
    [responsiveCells_S(c), ~, ~, ~] = checkIfActive_v2(avg_S, nBaselineFrames, STDlevel, AUC_S_level, 0);
end

disp(['Found ' num2str(sum(responsiveCells_F)) ' responsive F traces and ' num2str(sum(responsiveCells_S)) ' responsive spike traces'])

responsiveCells = find(or(responsiveCells_F, responsiveCells_S));
responsiveCellNums = block.cell_number(responsiveCells);

%% Plot raw activity of cell(s) for duration of block
%Raw activity of each cell vs. time with locomoter activity beneath
%Each cell is a separate trace

figure('units','normalized','outerposition',[0 0 1 1])
count = 1; %for staggering plot lines

for c = 1:length(responsiveCellNums)
    current_cellnum = responsiveCellNums(c);

    subplot(3,4,1:8); hold on

    row_num = find(all_cell_numbers == current_cellnum);
    
    cell_trace = F7_full(row_num,:);%pull out the full trace for each cell

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
set(gca, 'YTickLabel', [responsiveCellNums(1:count-1)])
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
try
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
catch
    disp('Loco data not plotted. Recompile block to get latest loco data.')
end
    
%% Plot graphs according to stim presentation

% Average response to all stim
for f = 1:size(responsiveCellNums,1) %Individual figures if cellnum is vertical
        current_cells = responsiveCellNums(f,:);
        
        row_nums = nan(length(current_cells),1);
        for c = 1:length(current_cells)
            row_nums(c) = find(all_cell_numbers == current_cells(c));
        end
    
        if length(current_cells) == 1
            F7_cell = squeeze(F7_stim(row_nums,:,:));
            F7_baseline = F7_cell(:,1:nBaselineFrames); %baseline for each trial
            F7_df_f = (F7_cell-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
            ebar = std(F7_df_f,1);
            spks_cell = squeeze(spks_stim(row_nums,:,:));
            ebar_spks = std(spks_cell,1);
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
        shadedErrorBar(1:length(total_mean),total_mean,ebar);
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
                strcat(' Cell #', num2str(responsiveCellNums(f)))])
        else
            suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(responsiveCellNums)), ' cells')])
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
    
    for f = 1:size(responsiveCellNums,1) %Individual figures if cellnum is vertical
        current_cells = responsiveCellNums(f,:);
        
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
                ebar = std(F7_df_f,1);
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
            shadedErrorBar(1:length(total_mean),total_mean,ebar);
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
                
                stim_rows = intersect(find(V1 == V1_stim(vv)), find(V2 == V2_stim(v)));
                V1_list(count) = V1_stim(vv);
                V2_list(count) = V2_stim(v);
                
                if length(current_cells) == 1
                    F7_cell = squeeze(F7_stim(row_nums,stim_rows,:));
                    F7_baseline = F7_cell(:,1:nBaselineFrames); %baseline for each trial
                    F7_df_f = (F7_cell-nanmean(F7_baseline,2))./nanmean(F7_baseline,2); %(total-mean)/mean
                    ebar = std(F7_df_f,1);
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
            
            %artificial ylim if no max was found
            if max_spks == 0
                max_spks = 10;
            end
            if max_df_f == 0
                max_df_f = 1;
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
                shadedErrorBar(1:length(total_mean),total_mean,ebar);
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
                strcat(' Cell #', num2str(responsiveCellNums(f)))])
            else
                suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(responsiveCellNums)), ' cells')])
            end 
        end
                
        %Plot V1 x V2 grid using DF/F or Spikes
        if length(V1_stim) > 1 && length(V2_stim) > 1 %Don't plot if there isn't a grid

            nPlots = length(V1_stim)*length(V2_stim);
            
            %Find max values for y limits
            max_df_f = max(max(F7_df_f_mat,[], 1)) + max(max(ebar_mat,[], 1));
            max_spks = max(max(spks_mat,[], 1));

            %artificial ylim if no spikes
            if isnan(max_spks)
                max_spks = 10;
            end
            
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
                strcat(' Cell #', num2str(responsiveCellNums(f)))])
            else
                suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(responsiveCellNums)), ' cells')])
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
                strcat(' Cell #', num2str(responsiveCellNums(f)))])
            else
                suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(responsiveCellNums)), ' cells')])
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
                    strcat(' Cell #', num2str(responsiveCellNums(f)))])
            else
                suptitle([block.setup.block_supname...
                    strcat('Average of ', num2str(length(responsiveCellNums)), ' cells')])
            end            
        end   
    end
end

end %function end

