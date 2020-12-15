%plot_ROI_data_from_Fall

Fall = load('Fall.mat'); %Must load like this because iscell is a matlab function and might lead to unexpected errors.

cellsToPlot = [11, 11, 11,11,11,11,11,11,11,11,11];
runIDsToPlot = [6,7,8,9,10,11,5,2,3,4,1];
power = [10,20,30,50,70,80,100,30,50,80,100];
peak = nan(1,3);

for r = 1:length(runIDsToPlot)
    
    currentRunID = runIDsToPlot(r);
    currentCellToPlot = cellsToPlot(r);

    %% Extract info by cycle
    currentRunName = data.Folders{currentRunID,1};
    xMat = data.xMats{currentRunID,1};
    nFrames = size(xMat,2);

    % Cycles = Baseline -> Mark points -> Activation... Repeat
    cycleNumbers = unique(xMat(1,:));
    baseline = 1:3:(cycleNumbers(end)); %t-series always starts with baseline
    activation = 3:3:(cycleNumbers(end));

    %Set new baseline vs. activation parameter
    xMat(3,:) = 0;
    for i = 1:length(activation)
        xMat(3,xMat(1,:) == activation(i)) = 1;
    end

    baselineDuration = max(xMat(2,xMat(3,:) == 0));
    activationDuration = max(xMat(2,xMat(3,:) == 1));
    markPoints = baselineDuration:(baselineDuration + activationDuration):nFrames;

    %% Extract ROI data from Fall

    [Frame_set,~] = get_frames_from_Fall(Fall.ops, currentRunName, 1);
    %Cell and neuropil data
    %Only keep data from 'is cells' within the block's frame set
    block.iscell = Fall.iscell;
    keep_ind = find(block.iscell(:,1));
    block.cell_number = keep_ind-1; %Subtract 1 for the python to matlab correction of cell label
    block.stat = Fall.stat(1,keep_ind);
    block.F = Fall.F(keep_ind,Frame_set);
    block.Fneu = Fall.Fneu(keep_ind,Frame_set);
    block.spks = Fall.spks(keep_ind,Frame_set);
    block.F7 = block.F - 0.7*block.Fneu;

    if isfield(Fall, 'redcell')
        redcell = Fall.redcell; %Not all runs will have red cells
        block.redcell = redcell(keep_ind);
    else
        block.redcell = nan;
    end

%% Plot graph of raw traces

%DF/F
mean_of_full_trace = mean(block.F7,2);
block.F7_df_f = (block.F7-mean_of_full_trace)./mean_of_full_trace; %(total-mean)/mean

count = 0;
% 
% figure; hold on
% for p = 1:size(block.F7_df_f,1) 
%     plot((1:nFrames),block.F7_df_f(p,1:nFrames) + count)
%     count = count + 2;
% end
% title(currentRunName)
% vline(markPoints)
% %ylim([0 21])
% xlim([0,nFrames])

%% Make average mats

nCells = size(block.F7_df_f,1);

baselineMat = nan(nCells,baselineDuration,length(baseline));
activationMat = nan(nCells,activationDuration,length(baseline));
 
for b = 1:length(baseline)
    for c = 1:nCells
    baselineMat(c,:,b) = block.F7_df_f(c, find(xMat(1,:) == baseline(b)));
    end
end

for a = 1:length(activation)
    for c = 1:nCells
    activationMat(c,:,a) = block.F7_df_f(c, find(xMat(1,:) == activation(a)));
    end
end

%average by trial
baselineMat = mean(baselineMat,3);
activationMat = mean(activationMat,3);
avgTrialMat = [baselineMat, activationMat];

activation_peaks = max(activationMat,[],2);
cellInd = find(block.cell_number == currentCellToPlot);
peak(r) = activation_peaks(cellInd);



% 
% %Figure
% nPlots = floor(sqrt(nCells));
% 
% figure; subplot(nPlots,nPlots,1);
% suptitle(currentRunName)
% for c = 1:nCells
%     if c > nPlots^2
%         continue
%     end
%     subplot(nPlots,nPlots,c); hold on
%     title(['Cell #' num2str(block.cell_number(c))])
%     plot(1:size(avgTrialMat,2),avgTrialMat(c,:))
%     line([baselineDuration, baselineDuration], [0 2], 'Color', 'r')
%     xlim([0 size(avgTrialMat,2)])
%     ylabel('F7')
%     xlabel('Frames')
% end

end

figure
bar(peak)
ylabel('Peak dff')
xlabel('Power')
set(gca, 'XTick', 1:length(power))
set(gca, 'XTickLabel', power)
title('YE083020F2 Cell 2')