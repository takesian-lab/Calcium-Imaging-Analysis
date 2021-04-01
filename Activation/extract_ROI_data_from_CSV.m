%extract_ROI_data_from_CSV

%% Make struct with extracted ROI traces
% %All_ROIs = {};
% raw_ROIs = [];
% %copy-paste data from excel into ROIs
% currentRunID = 12; %CHECK data.Folders
% All_ROIs{currentRunID,1} = raw_ROIs;

load('All_ROIs.mat')
load('TSeries_Data_17-Mar-2020-2139.mat')
currentRunID = 4;

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

%% Extract ROI data
%cell 1 is the last cell so flip and rotate ROIs
%cell structure is column 1 = cell 2, column 2 = neuropil for cell 1, etc.

raw_ROIs = All_ROIs{currentRunID,1};
ROIs = raw_ROIs';

%Plot graph of raw traces
figure
plot((1:nFrames),ROIs(:,1:nFrames))
title(currentRunName)
line([markPoints; markPoints],[zeros(1,length(markPoints)); zeros(1,length(markPoints))+21], 'Color', 'r')
ylim([0 21])
xlim([0,nFrames])

%Extract F7
nCells = size(ROIs,1)/2;
cell_ids = 1:2:size(ROIs,1);
neuropil_ids = cell_ids + 1;

cells = ROIs(cell_ids,:);
neuropil = ROIs(neuropil_ids,:);

for i = 1:nCells - 1
    cells(i,:) = cells(i,:) - (nCells - i);
    neuropil(i,:) = neuropil(i,:) - (nCells - i);
end

minCell = min(min(cells));
minNeuropil = min(min(neuropil));

cells = cells + abs(min(minCell, minNeuropil));
neuropil = neuropil + abs(min(minCell, minNeuropil));

F7 = cells - 0.7*neuropil;

%Plot graph of raw traces
figure
plot((1:nFrames),F7(:,1:nFrames))
title(currentRunName)
line([markPoints; markPoints],[zeros(1,length(markPoints)); zeros(1,length(markPoints))+21], 'Color', 'r')
ylim([0 2])

%% Make average mats
baselineMat = nan(nCells,baselineDuration,length(baseline));
activationMat = nan(nCells,activationDuration,length(baseline));
 
for b = 1:length(baseline)
    for c = 1:nCells
    baselineMat(c,:,b) = F7(c, find(xMat(1,:) == baseline(b)));
    end
end

for a = 1:length(activation)
    for c = 1:nCells
    activationMat(c,:,a) = F7(c, find(xMat(1,:) == activation(a)));
    end
end

%average by trial
baselineMat = mean(baselineMat,3);
activationMat = mean(activationMat,3);
avgTrialMat = [baselineMat, activationMat];

%Figure
figure; subplot(4,5,1);
suptitle(currentRunName)
for c = 1:nCells
    if c <= 15
        colour = 'g';
    else
        colour = 'k';
    end
    
    subplot(4,5,c); hold on
    title(['Cell #' num2str(c)])
    plot(1:size(avgTrialMat,2),avgTrialMat(c,:), colour)
    line([baselineDuration, baselineDuration], [0 2], 'Color', 'r')
    xlim([0 size(avgTrialMat,2)])
    ylabel('F7')
    xlabel('Frames')
end