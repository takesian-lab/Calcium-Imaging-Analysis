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

%Load roiMatchData file
disp('Load roiMatchData file')
[match_filename,match_filepath] = uigetfile('.mat');
load(match_filename)

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
all_ROIs_to_label = []; %For verification purposes

for f = 1:length(falls)
    current_allSessionMapping = roiMatchData.allSessionMapping(:,f);
    current_suite2p_labels = find(falls{1,f}.iscell(:,1)) - 1;
    ROI_to_label = [];
    for r = 1:length(current_allSessionMapping)
        ROI_to_label(r,1) = current_suite2p_labels(current_allSessionMapping(r));
    end
    all_ROIs_to_label = [all_ROIs_to_label, ROI_to_label];
    falls{1,f}.allSessionMapping = ROI_to_label;
    falls{1,f}.roiMatchData = roiMatchData;
    falls{1,f}.roiMatchData.fall_filenames = fall_filenames;
    falls{1,f}.roiMatchData.fall_filepaths = fall_filepaths;
end

%% Save fall.mats

for f = 1:length(falls)
    fullpath = [fall_filepaths{f} '/' fall_filenames{f}];
    fall = falls{1,f};
    save(fullpath, '-struct', 'fall');
end

disp('Finished processing!')

%If you want to check the labels against each other quickly,
%look at the variable all_ROIs_to_label