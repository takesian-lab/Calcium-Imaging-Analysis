function visualize_cell(block, cellnum)

% DOCUMENTATION IN PROGRESS
% 
% This function allows you to preview the data from a single block by
% plotting multiple types of figures
% 
% Argument(s): 
%   block (struct)
%   m (struct)
% 
% Returns:
%   
% 
% Notes:
%
%
% TODO: determine best way to measuer df/F. Currently, using mean trace as
% Fo; however, there are other (better?) ways to do this.
% Search 'TODO'

%% Magic numbers and setup

% neuCorrect = 0.7; %Neuropil correction coefficient
bin = 10; %Number of cells to plot at a time (for visibility)
SF = 0.5; %Shrinking factor for traces to appear more spread out
z = 1; %Portion of recording to plot e.g. 0.5, 0.33, 1

setup = block.setup;

%% Suite2p Section - SKIP if Fall.mat wasn't present

if ismissing(block.setup.suite2p_path)
    disp('Skipping Suite2p data plots...');
else
    
    if isfield(block, 'Sound_Time')
        Sound_Time = block.Sound_Time;
    end
    
    cell_number = block.cell_number;
    redcell = block.redcell;
    F = block.F; %all the cell fluorescence data
    Fneu = block.Fneu; %neuropil
    F7 = F-setup.constant.neucoeff*Fneu; %neuropil corrected traces

    if isfield(block, 'timestamp')
        timestamp = block.timestamp;
        timeUnit = 'Timestamp';
    else
        timestamp = 1:size(F7,2);
        timeUnit = 'Frames';
    end
    
    Z = round(length(timestamp)*z);
        
    %Divide into red and green cells
    %ones variable = row number
    %number variable = suite2p cell labels
    redcell_ones = find(redcell);

    if ~isempty(redcell_ones) 
        redcell_number = cell_number(redcell_ones);
        nonredcell_ones = setdiff(1:length(cell_number), redcell_ones); %what are the active non-red cells
        nonredcell_number = cell_number(nonredcell_ones);
        redcells_exist = 1; %for plotting
    else
        nonredcell_number = cell_number;
        nonredcell_ones = 1:length(cell_number);
        redcells_exist = 0;
    end

     %% Plot DF over F from cells (divided into red and green)
    for f = 1:length(cellnum)
        current_cellnum = cellnum(f);

        figure('units','normalized','outerposition',[0 0 1 1])

        subplot(3,4,1:8) %one cell/row of the graph

        row_num = find(cell_number == current_cellnum);
        count = 1; %for staggering plot lines

        cell_trace = F7(row_num,:);%pull out the total trace for each cell

        mean_gCAMP = mean(cell_trace);% average green for each cell
        df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;%(total-mean)/mean
        A = smooth(df_f,10);

        plot(timestamp, A*SF + count,'LineWidth',1);

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
    xlim([0 timestamp(Z)])
    ylabel('DF/F')
    xlabel(timeUnit)


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
    suptitle(block.setup.block_supname)
    end
        
end

F7_stim = block.aligned_stim.F7_stim;

for f = 1:length(cellnum)
        current_cellnum = cellnum(f);
        row_num = find(cell_number == current_cellnum);
        F7_cell = squeeze(F7_stim(row_num,:,:));

        figure; hold on
        subplot(1,2,1)
        imagesc(F7_cell)
        ylabel('DF/F')
        xlabel('Frames')
        
        subplot(1,2,2)
        total_mean = mean(F7_cell);
        plot(total_mean);
end

%% Air puff


F7_stim = block.aligned_stim.F7_stim;
V1 = block.parameters.variable1;

for f = 1:length(cellnum)
        current_cellnum = cellnum(f);
        row_num = find(cell_number == current_cellnum);
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
        row_num = find(cell_number == current_cellnum);
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

