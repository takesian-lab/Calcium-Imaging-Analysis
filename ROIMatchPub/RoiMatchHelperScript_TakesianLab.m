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

%% Load files
clear all

batchProcessing = 1;
batchLabels = [];
b = 1;

while batchProcessing
    
%Load roiMatchData file
disp('Load roiMatchData file')
[match_filename,match_filepath] = uigetfile('.mat');
load(match_filename)
disp(match_filename)
fall1Number = input('Input first fall number: ');
fall2Number = input('Input second fall number: ');
batchLabels(b,1) = fall1Number;
batchLabels(b,2) = fall2Number;

%Load Fall.mat files in the same order they were loaded for roiMatchPub
%We need a minimum of 2 files
falls = {};
fall_filenames = {};
fall_filepaths = {};
nFilesLoaded = 1;
allFilesLoaded = 0;
while ~allFilesLoaded
    if nFilesLoaded == 1
        disp('Load first fall.mat file')
    elseif nFilesLoaded == 2
        disp('Load second fall.mat file')
    else
        disp('Load additional fall.mat file')
    end
    [fall_filenames{nFilesLoaded},fall_filepaths{nFilesLoaded}] = uigetfile('.mat');
    falls{nFilesLoaded} = load(fall_filenames{nFilesLoaded});
    disp(fall_filenames{nFilesLoaded});
    nFilesLoaded = nFilesLoaded + 1;
    if nFilesLoaded > 2
        allFilesLoaded = input('Done loading files? 1 for yes, 0 for no: ');
    end
end

%% Check that loaded files match the roiMatchPub (as best as we can)
% I can think of two ways to do this:
% 1) check that the correct number of fall.mats were loaded
% 2) check that the fall.mat has the same number of ROIs
% 3) check that the mean image is identical

% 1) check that the correct number of fall.mats were loaded
if length(falls) ~= size(roiMatchData.allSessionMapping,2)
    error('Number of fall.mat files does not match the roiMatchPub file')
end

for f = 1:length(falls)
    % 2) check that the fall.mat has the same number of ROIs
    nROIs = sum(falls{1,f}.iscell(:,1));
    cellCount = roiMatchData.rois{1,f}.cellCount;
    if nROIs ~= cellCount
        error('One of the loaded fall.mats does not match the roiMatchPub file')
    end
    
    % 3) check that the mean image is identical
    meanImg = falls{1,f}.ops.meanImg;
    meanFrame = roiMatchData.rois{1,f}.meanFrame;
    if ~isequal(meanImg,meanFrame)
        error('One of the loaded fall.mats does not match the roiMatchPub file')
    end   
end

%% Convert ROI mapping to suite2p labels

allSessionMapping = roiMatchData.allSessionMapping;
all_ROIs_to_label = []; %Convert ROI to suite2p labels
all_manualMatches = []; %Manual matches added by user during manual matching (Be more careful with these)

for f = 1:length(falls)
    
    current_suite2p_labels = find(falls{1,f}.iscell(:,1)) - 1;
    ROI_to_label = [];
    ROI_to_label_new = [];
        
    %Use manual labelling if available
    if isfield(roiMatchData, 'manualMatching')
        ROI_to_label = roiMatchData.manualMatching.matchedCells{f,1} - 1; %Corrected auto-matches, -1 is for python correction
        
        %If new manual matches exist, save to new variable to not mix with auto-matches
        if ~isempty(roiMatchData.manualMatching.newMatches) %Manual matches
            ROI_to_label_new = roiMatchData.manualMatching.newMatches(:,f) - 1; %-1 is for python correction
            all_manualMatches = [all_manualMatches, ROI_to_label_new];
        else
            ROI_to_label_new = nan;
            all_manualMatches = [all_manualMatches, ROI_to_label_new];
        end
            
    else
        current_allSessionMapping = roiMatchData.allSessionMapping(:,f); %original matches
        for r = 1:size(current_allSessionMapping,1)
            ROI_to_label(r,1) = current_suite2p_labels(current_allSessionMapping(r));
        end
        
        ROI_to_label_new = nan;
        all_manualMatches = [all_manualMatches, ROI_to_label_new];
    end
    
    all_ROIs_to_label = [all_ROIs_to_label, ROI_to_label];
    
    falls{1,f}.allSessionMapping = ROI_to_label;
    falls{1,f}.roiMatchData = roiMatchData;
    falls{1,f}.roiMatchData.fall_filenames = fall_filenames;
    falls{1,f}.roiMatchData.fall_filepaths = fall_filepaths;
end

ROIs_to_label{1,b} = all_ROIs_to_label;
ManualROIs_to_label{1,b} = all_manualMatches;
b = b + 1; %batch loop count
batchProcessing = input('Continue batch processing? 1 for yes, 0 for no: ');

end %while loop

%% Combine multiple ROIs_to_label and ManualROIs_to_label using batchLabels

for r = 1:2
    if r == 1
        labels = ROIs_to_label;
    elseif r == 2
        labels = ManualROIs_to_label;
    end
    
    combinedLabels = []; 
    columns = []; 
    firstSetOfLabels = 1;
    
    for i = 1:size(batchLabels,1)
        
        if isempty(labels{i})
            labels{i} = [nan, nan]; %fake data so that column numbers will stay the same
        end
           
        file1 = batchLabels(i,1);
        file2 = batchLabels(i,2);
        
        if firstSetOfLabels
            combinedLabels(:,file1) = labels{i}(:,1);
            combinedLabels(:,file2) = labels{i}(:,2);
            columns(1) = file1;
            columns(2) = file2;
            column_count = 3;
            firstSetOfLabels = 0;
            continue
        end

        column_index1 = find(columns == file1);
        column_index2 = find(columns == file2);
        if ~isempty(column_index1) && isempty(column_index2) %File 1 matches previous file
            existingList = labels{i}(:,1);
            newList = labels{i}(:,2);
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
            existingList1 = labels{i}(:,1);
            existingList2 = labels{i}(:,2);
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
 
    if r == 1
        combined_ROIs_to_label = combinedLabels;
    elseif r == 2
        combined_ManualROIs_to_label = combinedLabels;
    end
end

% Save combined matching data

save('CombinedMatchingData', 'batchLabels', 'combined_ROIs_to_label', 'combined_ManualROIs_to_label');


%% Save fall.mats

% for f = 1:length(falls)
%     fullpath = [fall_filepaths{f} '/' fall_filenames{f}];
%     fall = falls{1,f};
%     save(fullpath, '-struct', 'fall');
% end
% 
% disp('Finished processing!')

%If you want to check the labels against each other quickly,
%look at the variable all_ROIs_to_label