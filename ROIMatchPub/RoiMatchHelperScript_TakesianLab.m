%ROIMatchHelperScript_TakesianLab

%Understanding the roiMatchData struct:

%The first loaded Fall.mat = Experiment 1
%Every other Fall.mat = Experiment 2 to N
%# experiments = N
%# total ROIs = R (sum ROIs in N experiments)
%# found ROIs = r

%roiMatchData = struct with 6 fields:
%rois
%mapping
%comparisonMatrix
%allRois
%refImage
%allSessionMapping

%rois = 1 x N cell containing one struct per experiment
%cellCount = # of ROIs in Fall.mat
%meanFrame = Mean suite2p image from Fall.ops.meanImg
%meanFrameRegistered = meanFrame if experiment 1
%... If not, the meanFrame image is warped to match meanFrame of experiment 1 based on ROI matching
%roiMapRegistered = ROI masks built from Fall.stat. If experiment 1, these are the original ROI locations
%... If not, the ROI masks are warped as above
%committed = A binary (0 or 1) of which ROIs were found to have matches
%... sum(committed) from each experiment is the same

%mapping = R x N double containing the matches for all ROIs
%... R = sum ROIs in N experiments
%... Every ROI is listed for Exp 1 next to its closest match from the other exps
%... e.g. if Exp 1 has 66 Rois, rows 1:66 will be in order of the ROIs from exp 1
%... this is then concatenated vertically with the same thing for Exp 2, etc.
%... If an ROI has no match the row will be replaced with zeros
%... Each column corresponds to the experiments in the order they were loaded
%... The value of the columns corresponds to the ROI number
%... This is NOT the suite2p label

%comparisonMatrix = matrix to make sure every experiment is compared to
%... every other experiment once

%allRois = 1 x N cell containing Fall.mat filepaths

%refImage = meanFrame from experiment 1

%allSessionMapping = r x N double of the final ROIs chosen 
%... this is equivalent to 'mapping' except including only the final ROIs

%% Combine data across multiple roiMatchPub files

clear all

batchProcessing = 1;
batchFiles = {};
batchLabels = {};
b = 1;

while batchProcessing
    
%Load roiMatchData file
disp('Load roiMatchData file')
[match_filename,match_filepath] = uigetfile('.mat');
load([match_filepath '/' match_filename])
disp(match_filename)

%Load Fall.mat files from roiMatchData
%We need a minimum of 2 files
falls = {};
for f = 1:size(roiMatchData.allRois,2)
    falls{f} = load(roiMatchData.allRois{1,f});
    batchFiles{b,f} = roiMatchData.allRois{1,f};
end


%% Check that loaded fall.mats match the roiMatchPub (as best as we can)
% If the user modified the fall.mat ROIs anytime after matching, the suite2p cell
% labels might have changed. We don't have a great way to catch this except
% to check that total number of ROIs is the same (but total N could stay the same
% even if labels change so it's not perfect)

for f = 1:length(falls)
    nROIs = sum(falls{1,f}.iscell(:,1));
    cellCount = roiMatchData.rois{1,f}.cellCount;

    if nROIs ~= cellCount
        error('Suite2p labels have changed since matching was performed')
    end
end

%% Convert ROI mapping to suite2p labels

allSessionMapping = roiMatchData.allSessionMapping;
all_ROIs_to_label = []; %Convert ROI to suite2p labels
unmatchedCells = {};

for f = 1:length(falls)

    current_suite2p_labels = find(falls{1,f}.iscell(:,1)) - 1;
    ROI_to_label = [];
        
    %Use manual labelling if available
    if isfield(roiMatchData, 'manualMatching')
        ROI_to_label = roiMatchData.manualMatching.matchedCells{f,1} - 1; %Corrected auto-matches, -1 is for python correction
        
        %New manual matches
        if ~isempty(roiMatchData.manualMatching.newMatches)
            newMatches = roiMatchData.manualMatching.newMatches(:,f) - 1; %-1 is for python correction
            ROI_to_label = [ROI_to_label; newMatches];
        end
            
    else
        current_allSessionMapping = roiMatchData.allSessionMapping(:,f); %original matches
        for r = 1:size(current_allSessionMapping,1)
            ROI_to_label(r,1) = current_suite2p_labels(current_allSessionMapping(r));
        end

    end
    
    %Add any cells that didn't have a match
    unmatchedCells{f} = setdiff(current_suite2p_labels, ROI_to_label);
    
    all_ROIs_to_label = [all_ROIs_to_label, ROI_to_label];
end

%Append unmatched cells to ROI labels by padding with nans
for f = 1:length(falls)
    current_unmatchedCells = unmatchedCells{f};
    nanmat = nan(length(current_unmatchedCells),length(falls));
    nanmat(:,f) = current_unmatchedCells;
    all_ROIs_to_label = [all_ROIs_to_label; nanmat];
end

batchLabels{b,1} = all_ROIs_to_label;
b = b + 1; %batch loop count
batchProcessing = input('Continue batch processing? 1 for yes, 0 for no: ');

end %while loop

%% Combine multiple ROIs_to_label using batchLabels

for b = 1:size(batchLabels,1)
    
    if b == 1
        combinedLabels = batchLabels{1}; 
        columns = batchFiles(1,:);
        continue
    end
    
    if isempty(batchLabels{b})
        batchLabels{b} = [nan, nan]; %fake data so that column numbers will stay the same
    end
    
    if size(batchLabels{b},2) > 2
        error('This part of the script only works with 2-file matches right now')
    end
    
    column_index1 = find(strcmp(batchFiles{b,1}, columns));
    column_index2 = find(strcmp(batchFiles{b,2}, columns));

    if isempty(column_index1) && isempty(column_index2) %Files do not match previous files
        %Pad combinedLabels with NaNs
        c1 = size(combinedLabels,2) + 1;
        c2 = c1 + 1;
        combinedLabels(:,c1:c2) = nan;
        %Pad new labels with NaNs
        nanmat = nan(size(batchLabels{b},1),c2);
        nanmat(:,c1:c2) = batchLabels{b};
        %Append new labels to the end of combinedLabels
        combinedLabels = [combinedLabels; nanmat];
    
    elseif ~isempty(column_index1) || isempty(column_index2) %One file matches previous file
        existingList = batchLabels{b}(:,1);
        newList = labels{b}(:,2);
        combinedLabels(:,column_count) = nan;
        for j = 1:length(existingList)
            if ismember(existingList(j),combinedLabels(:,column_index1))
                list_ind = find(combinedLabels(:,column_index1) == existingList(j));
                combinedLabels(list_ind,column_count) = newList(j);
            else
                listLength = size(combinedLabels,1);
                combinedLabels(listLength + 1, :) = nan;
                combinedLabels(listLength + 1, column_index1) = existingList(j);
                combinedLabels(listLength + 1, column_count) = newList(j);
            end
        end
        columns(column_count) = file2;
        column_count = column_count + 1;

    elseif ~isempty(column_index1) && ~isempty(column_index2) %Both File 1 and 2 match previous files
        existingList1 = labels{b}(:,1);
        existingList2 = labels{b}(:,2);
        for j = 1:length(existingList1)
            if ismember(existingList1(j),combinedLabels(:,column_index1))
                list_ind = find(combinedLabels(:,column_index1) == existingList1(j));
                label_to_check = combinedLabels(list_ind, column_index2);
                if isnan(label_to_check)
                    combinedLabels(list_ind, column_index2) = existingList2(j);
                elseif label_to_check == existingList2(j) %labels match
                    continue
                else %labels don't match
                    listLength = size(combinedLabels,1);
                    combinedLabels(listLength + 1, :) = nan;
                    combinedLabels(listLength + 1, column_index1) = existingList1(j);
                    combinedLabels(listLength + 1, column_index2) = existingList2(j);
                end
            else
                listLength = size(combinedLabels,1);
                combinedLabels(listLength + 1, :) = nan;
                combinedLabels(listLength + 1, column_index1) = existingList1(j);
                combinedLabels(listLength + 1, column_index2) = existingList2(j);                
            end
        end
    end
end


% Save combined matching data
save('CombinedMatchingData', 'batchLabels', 'combinedLlabel');


%If you want to check the labels against each other quickly,
%look at the variable all_ROIs_to_label