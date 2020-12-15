function visualize_block(block)
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

%% Behavior Section - SKIP if Tosca data wasn't present
if ismissing(block.setup.Tosca_path)
    disp('Skipping Tosca plots...');
else

    %% Plot locomotor activity

    %THIS IS TEMPORARILY MESSED UP
    if isfield(block, 'loco_activity')
        loco_time = block.loco_times;
        loco_speed = block.loco_activity;
        active_time = loco_speed;
        active_time(active_time < block.setup.constant.locoThresh) = 0;

%     elseif isfield(block, 'locomotion_data')
%         loco_time = block.locomotion_data(:,1);
%         loco_speed = block.locomotion_data(:,3);
%     else
%         loco_time = block.loco_data(:,1);
%         loco_speed = block.loco_data(:,3);
    end
%         
%     if isfield(block, 'active_time')
%         active_time = block.active_time;
%     end

    figure; 

    subplot(2,1,1); hold on
    title('Locomotor activity')
    ylabel('Activity (cm/s)')
    plot(loco_time, loco_speed); hold on
    hline(0.7);

    subplot(2,1,2); hold on
    ylabel('Considered active')
    xlabel('Seconds')
    set(gca, 'ytick', [0 1])
    %if isfield(block, 'active_time') 
        plot(loco_time, active_time > 0); hold on;
    %end
    suptitle(block.setup.block_supname)
   

end %Skip if Tosca data is missing

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

    
    %% Plot raster from spikes (divided into red and green)
   
    spikes = block.spks;
    
    for f = 1:2
        if f == 1
            currentCells = nonredcell_ones;
            currentNumbers = nonredcell_number;
            currentSpikes = spikes(currentCells,:);
            fig_title = 'Green cells raster';
        elseif f == 2 && redcells_exist
            currentCells = redcell_ones;
            currentNumbers = redcell_number;
            currentSpikes = spikes(currentCells,:);
            fig_title = 'Red cells raster';
        else
            continue
        end
        
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2,1,1);
        %colormap('bone')
        imagesc(currentSpikes);
        xlim([0 timestamp(Z)])
        set(gca, 'YTick', [1:1:length(currentNumbers)])
        set(gca, 'YTickLabel', currentNumbers)
        ylabel('Cell number')
        xlabel(timeUnit)
        title(fig_title)
        h = colorbar;
        set(get(h,'label'),'string','Deconvolved');
        
        subplot(2,1,2);
        normSpikes = bsxfun(@rdivide, currentSpikes, max(currentSpikes,[],2));
        imagesc(normSpikes);
        xlim([0 timestamp(Z)])
        set(gca, 'YTick', [1:1:length(currentNumbers)])
        set(gca, 'YTickLabel', currentNumbers)
        ylabel('Cell number')
        xlabel(timeUnit)
        title(fig_title)
        h = colorbar;
        set(get(h,'label'),'string','Deconvolved normalized');
        suptitle(block.setup.block_supname)
    end
    
    
    %% Plot DF over F from cells (divided into red and green)
    for f = 1:2
        if f == 1
            currentCells = nonredcell_ones;
            currentNumbers = nonredcell_number;
            fig_title = 'Green cells - DF over F neuropil corrected';
        elseif f == 2 && redcells_exist
            currentCells = redcell_ones;
            currentNumbers = redcell_number;
            fig_title = 'Red cells - DF over F neuropil corrected';
        else
            continue
        end

        B = floor(length(currentCells)/bin);
        if B == 0
            B = 1;
        end
        extraCells = mod(length(currentCells),bin);
        if extraCells <= 5
            lastBin = bin + extraCells;
        else
            lastBin = extraCells;
        end

        for b = 1:B %Plot one graph for each bin of cells

            c1 = bin*(b-1) + 1; %first cell for this figure
            c2 = bin*(b); %last cell for this figure

            if b == B
                c2 = length(currentCells);
            end

            figure('units','normalized','outerposition',[0 0 1 1])

            subplot(3,4,1:8) %one cell/row of the graph
            
            count = 1; %for staggering plot lines
            
            for a=c1:c2 %for all of the cells in the current bin
                row_num = currentCells(a);
                cell_trace = F7(row_num,:);%pull out the total trace for each cell
                
                mean_gCAMP = mean(cell_trace);% average green for each cell
                df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;%(total-mean)/mean
                A = smooth(df_f,10);

                plot(timestamp, A*SF + count,'LineWidth',1);
                hold on;

                count = count + 1;
            end
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
            ylim([0 (count - 0.5)])
            set(gca, 'YTick', [1:1:count-1])
            set(gca, 'YTickLabel', [currentNumbers(c1:c2)])
            ylabel('Cell number')
            xlabel(timeUnit)
            title(fig_title)

            subplot(3,4,9:12); hold on %loco
            title('Locomotor activity')
            ylabel('Activity')
            xlabel(timeUnit)
            xlim([0 timestamp(Z)])
            if ~ismissing(block.setup.Tosca_path)
                plot(loco_time, loco_speed);
            end
            suptitle(block.setup.block_supname)
        end
        
    end

    %% Plot mean image from suite2p with ROIs outlined and labelled

    stat = block.stat;
    stat = stat(:,2:end); %python to matlab correction

    %Plot maps twice, second will have ROI labels,
    %third will have red cell labels
    for f = 1:3
        if f == 1
            plotROIs = 0;
        elseif f == 2
            plotROIs = 1;
            currentCells = nonredcell_ones;
            currentNumbers = nonredcell_number;
            fig_title = 'Green cells';
        elseif f == 3 && redcells_exist
            plotROIs = 1;
            currentCells = redcell_ones;
            currentNumbers = redcell_number;
            fig_title = 'Red cells';
        else
            continue
        end
 
        figure('units','normalized','outerposition',[0 0 1 1])
        for m = 1:2

            subplot(1,2,m); hold on

            if m == 1
                figtitle = 'Ref Image';
                img = block.img.refImg;
                %img = block.img.meanImg;
                image(img);
                if ismissing(block.setup.VR_name)
                    imagesc(img);
                end
            elseif m == 2 
                figtitle = 'Mean Image Enhanced';
                img = block.img.meanImgE;
                imagesc(img);
            end

            title(figtitle)
            axis square
            xlim([0 512])
            ylim([0 512])
            colormap('bone')
            set(gca,'YDir','reverse')
        end

        if plotROIs == 1

            for a = 1:length(currentCells)
                row_num = currentCells(a);
                xcirc = double(block.stat{1,row_num}.xcirc);
                ycirc = double(block.stat{1,row_num}.ycirc);

                subplot(1,2,1);
                if f == 2
                    plot(xcirc,ycirc,'Linewidth', 1.5);
                elseif f == 3
                    plot(xcirc,ycirc,'Linewidth', 1.5, 'Color', 'r');
                end
                text(max(xcirc),max(ycirc),num2str(currentNumbers(a)), 'Color', 'w');

                subplot(1,2,2);
                if f == 2
                    plot(xcirc,ycirc,'Linewidth', 1.5);
                elseif f == 3
                    plot(xcirc,ycirc,'Linewidth', 1.5, 'Color', 'r');
                end
                text(max(xcirc),max(ycirc),num2str(currentNumbers(a)), 'Color', 'c');
            end

            subplot(1,2,1);
            set(gca,'YDir','reverse')

            subplot(1,2,2);
            set(gca,'YDir','reverse')
            suptitle(block.setup.block_supname)
        end
    end 
    

end %Skip if Suite2p data is missing

end %end function
