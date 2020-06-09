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

%% Magic numbers

SF = 0.5; %Shrinking factor for traces to appear more spread out
z = 1; %Portion of recording to plot e.g. 0.5, 0.33, 1


%% Setup

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

%% Plot raw activity of cell(s) for duration of block

Sound_Time = block.Sound_Time;
all_cell_numbers = block.cell_number;    
F = block.F; %all the cell fluorescence data
Fneu = block.Fneu; %neuropil
F7 = F-setup.constant.neucoeff*Fneu; %neuropil corrected traces
timestamp = block.timestamp;
timeUnit = 'Timestamp';
Z = round(length(timestamp)*z);
    
for p = 1:2 %Two plots
    if p == 2 && length(cellnum) == 1
        continue
    end
        
    figure('units','normalized','outerposition',[0 0 1 1])
    count = 1; %for staggering plot lines

    if p == 1
        %PLOT 1 - Raw activity of each cell vs. time with locomoter activity beneath
        %For this graph only, plot one figure with each cell as a separate trace
        for c = 1:length(cellnum) %Cells to average together
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
        
            if length(cellnum) == 1
                ylabel('DF/F')
            else
                set(gca, 'YTick', [1:1:count-1])
                set(gca, 'YTickLabel', [cellnum(1:count-1)])
            end
            suptitle(block.setup.block_supname)

    elseif p == 2
        %PLOT 2 - Averaged raw activity of each cell vs. time with locomoter activity beneath
            subplot(3,4,1:8); hold on
            
            a = nan(length(cellnum),size(F7,2));
            for c = 1:length(cellnum)
                row_num = find(all_cell_numbers == cellnum(c));
                cell_trace = F7(row_num,:);%pull out the full trace for each cell
                mean_gCAMP = mean(cell_trace);% average green for each cell
                df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;%(total-mean)/mean]
                a(c,:) = smooth(df_f,10);
            end
            
            A = mean(a,1);
            plot(timestamp, A*SF,'LineWidth',1);
            ylabel('DF/F')
            suptitle(strcat(block.setup.block_supname, ' - Average of ', num2str(length(cellnum)), ' cells'))
    end
    
    xlim([0 timestamp(Z)])
    xlabel(timeUnit)
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

    subplot(3,4,9:12); hold on %loco
    title('Locomotor activity')
    ylabel('Activity')
    xlabel(timeUnit)
    xlim([0 timestamp(Z)])
    if ~ismissing(block.setup.Tosca_path)
        if isfield(block, 'locomotion_data')
            loco_data = block.locomotion_data;
        else
            loco_data = block.loco_data;
        end
        plot(loco_data(:,1), loco_data(:,3));
    end
end   
    
%% Plot graphs according to stim presentation

F7_stim = block.aligned_stim.F7_stim;

% Plot 1 - average response to all stim
for f = 1:size(cellnum,1) %Individual figures if cellnum is vertical
        current_cells = cellnum(f,:);
        
        row_nums = nan(length(current_cells),1);
        for c = 1:length(current_cells)
            row_nums(c) = find(all_cell_numbers == current_cells(c));
        end
    
        if length(current_cells) == 1
            F7_cell = squeeze(F7_stim(row_nums,:,:));
            ebar = std(F7_cell,1);
        elseif length(current_cells) > 1
            F7_cell = squeeze(mean(F7_stim(row_nums,:,:),1));
            ebar = std(F7_cell,1)/sqrt(length(current_cells)-1);
        end

        figure; hold on
        subplot(1,2,1)
        imagesc(F7_cell)
        ylabel('Trials')
        xlabel('Frames')
        
        subplot(1,2,2)
        total_mean = mean(F7_cell);
        shadedErrorBar(1:length(total_mean),total_mean,ebar);
        ylabel('DF/F')
        xlabel('Frames')
        xlim([0 length(total_mean)])
        
        if length(current_cells) == 1
            suptitle([block.setup.block_supname...
                strcat(' Cell #', num2str(cellnum(f)))])
        else
            suptitle([block.setup.block_supname...
                strcat('Average of ', num2str(length(cellnum)), ' cells')])
        end
        
end

%% Plot 2 - average response separated by stim type

if stim_protocol == 9 %Random H20
    var1 = 'H20';
    var0 = 'No H20';
elseif stim_protocol == 11 %Random Air
    var1 = 'Air';
    var0 = 'No Air';   
else
    continue
end

F7_stim = block.aligned_stim.F7_stim;
V1 = block.parameters.variable1;

for f = 1:length(cellnum)
        current_cellnum = cellnum(f);
        row_num = find(all_cell_numbers == current_cellnum);
        F7_cell = squeeze(F7_stim(row_num,:,:));

        figure; hold on
        subplot(1,3,1)
        imagesc(F7_cell)
        ylabel('Trials')
        xlabel('Frames')
        
        subplot(1,3,2)
        total_mean = mean(F7_cell(V1 == 1,:));
        plot(total_mean);
        ylim([0 1800])
        ylabel('DF/F')
        xlabel('Frames')
        title('H20')
        
        subplot(1,3,3)
        total_mean = mean(F7_cell(V1 == 0,:));
        plot(total_mean);
        ylim([0 1800])
        ylabel('DF/F')
        xlabel('Frames')
        title('No H20')
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

