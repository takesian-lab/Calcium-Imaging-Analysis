function visualize_block(block)
% Preview the data from a single block

setup = block.setup;

%% Plot locomotor activity

loco_data = block.locomotion_data;
active_time = block.active_time;

figure;

subplot(2,1,1); hold on
title('Locomotor activity')
ylabel('Activity')
plot(loco_data(:,1), loco_data(:,3));

subplot(2,1,2); hold on
ylabel('Considered active')
xlabel('Seconds')
set(gca, 'ytick', [0 1])
plot(loco_data(:,1), active_time > 0); hold on;

%% Plot activity from cells (divided into red and green)

bin = 10; %Number of cells to plot at a time (for visibility)

timestamp = block.timestamp;
Sound_Time = block.Sound_Time;

cell = block.F; %all the cell fluorescence data
Fneu = block.Fneu; %neuropil
iscell = block.iscell;
redcell = block.redcell;
F7 = cell-0.7*Fneu; %neuropil corrected traces
F7 = F7(2:end,:); %python to matlab correction

%Identify cells
cell_id = find(iscell(:,1));
redcell_ones = find(redcell(:,1));

redcells_exist = 0;

if ~isempty(redcell_ones)
    idx = find(ismember(redcell_ones,cell_id)); %identify active red cells within all cells
    redcell_iscell = redcell_ones(idx); %which red cells are active
    nonredcell = setdiff(cell_id, redcell_iscell); %what are the active non-red cells
    
    %redcell_iscell = redcell_iscell - 1; %python to matlab correction
    %nonredcell = nonredcell - 1; %python to matlab correction
    
    redcells_exist = 1; %for plotting
else
    nonredcell = cell_id; % - 1;
end

%Plot one figure for green cells and one for red cells
for f = 1:2
    if f == 1
        currentCells = nonredcell;
        fig_title = 'Green cells';
    elseif f == 2 && redcells_exist
        currentCells = redcell_iscell;
        fig_title = 'Red cells';
    else
        continue
    end
    
    z = 1; %Portion of recording to plot e.g. 0.5, 0.33, 1
    Z = round(length(timestamp)*z);
    SF = 0.5; %Shrinking factor for traces to appear more spread out
    
    B = floor(length(currentCells)/bin);
    extraCells = mod(length(currentCells),bin);
    if extraCells <= 5
        lastBin = bin + extraCells;
    else
        lastBin = extraCells;
    end
    
    for b = 1:B %Plot one graph for each bin of cells
        
        c1 = bin*(b-1) + 1;
        c2 = bin*(b);
        
        if b == B
            c2 = length(currentCells);
        end

        figure('units','normalized','outerposition',[0 0 1 1])

        subplot(3,4,1:8) %one cell/row of the graph
        count = 1;
        
        for a=c1:c2 %for all of the cells
            row_num = currentCells(a);
            cell_trace = F7(row_num,:);%pull out the total trace for each cell
            mean_gCAMP = mean(cell_trace);% average green for each cell
            df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;%(total-mean)/mean
            A = smooth(df_f,10);
            plot(timestamp, A*SF + count,'LineWidth',1);
            hold on;
            
            count = count + 1;
        end
        if Z < 10000 %Don't plot red lines if there is too much data, otherwise its messy
            line = vline(Sound_Time, 'r');
        end
        xlim([0 timestamp(Z)])
        ylim([0 (count - 0.5)])
        set(gca, 'YTick', [1:1:count-1])
        set(gca, 'YTickLabel', [c1:1:c2])
        ylabel('Cell number')
        xlabel('Timestamp')
        title(fig_title)
        
        subplot(3,4,9:12); hold on %loco
        title('Locomotor activity')
        ylabel('Activity')
        xlabel('Timestamp')
        plot(loco_data(:,1), loco_data(:,3));
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
        currentCells = nonredcell;
    elseif f == 3&& redcells_exist
        plotROIs = 1;
        currentCells = redcell_iscell;
    else
        continue
    end
    
    figure('units','normalized','outerposition',[0 0 1 1])
    for m = 1:2

        subplot(1,2,m); hold on

        if m == 1
            figtitle = 'Ref Image';
            img = block.img.refImg;
            image(img);
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
            text(max(xcirc),max(ycirc),num2str(a), 'Color', 'w');

            subplot(1,2,2);
            if f == 2
                plot(xcirc,ycirc,'Linewidth', 1.5);
            elseif f == 3
                plot(xcirc,ycirc,'Linewidth', 1.5, 'Color', 'r');
            end
        end

        subplot(1,2,1);
        set(gca,'YDir','reverse')

        subplot(1,2,2);
        set(gca,'YDir','reverse')
    end
end
