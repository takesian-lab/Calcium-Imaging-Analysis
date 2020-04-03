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
    
    figure;  %one cell/row of the graph
    for a=1:length(currentCells) %for all of the cells

        row_num = currentCells(a);
        cell_trace = F7(row_num,:);%pull out the total trace for each cell
        mean_gCAMP = mean(cell_trace);% average green for each cell
        df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;%(total-mean)/mean
        A = smooth(df_f,10);
        plot(timestamp, A*SF + a,'LineWidth',1);
        hold on;

    end
    line = vline(Sound_Time);
    xlim([0 timestamp(Z)])
    ylim([0 a])
    ylabel('Cell number')
    xlabel('Timestamp')
    title(fig_title)
end

%% Plot mean image from suite2p with ROIs outlined and labelled

figure;

subplot(2,2,1)
imagesc(block.img.refImg)
axis square
colormap('bone')

subplot(2,2,2)
imagesc(block.img.meanImgE)
axis square
colormap('bone')

subplot(2,2,3); hold on
imagesc(block.img.refImg)
axis square
colormap('bone')

for a = 1:length(

