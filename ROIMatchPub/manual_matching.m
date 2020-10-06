%manual_matching
%Load roiMatchData file
disp('Load roiMatchData file')
[match_filename,match_filepath] = uigetfile('.mat');
load(match_filename)

filepaths = roiMatchData.allRois;
mapping = roiMatchData.allSessionMapping;

noMatchingROIs = 0;
if isempty(mapping)
    noMatchingROIs = 1;
end

%%
figure;

matchedCells = {};
unmatchedCells = {};
disp('Loading figure...')

for f = 1:length(filepaths)
    
    %Get ROI boundaries from fall.mat
    Fall = load(filepaths{1,f});
    cellValid = Fall.iscell(:,1);
    cellIDMap = zeros(size(Fall.ops.meanImg));
    validCellList = find(cellValid(:,1)==1);
    if ~noMatchingROIs
        matchedCells{f,1} = validCellList(mapping(:,f));
        unmatchedCells{f,1} = setdiff(validCellList,matchedCells{f,1});
    else
        unmatchedCells{f,1} = validCellList;
    end
    stat{f,1} = Fall.stat;

    subplot(1,2,f)
    meanFrame = Fall.ops.meanImg;
    imagesc(imadjust(int16(meanFrame)));
    colormap gray;
    hold on
    
end

%%

if ~noMatchingROIs
    
for i = 1:size(mapping,1)
    
    cellID_A = matchedCells{1,1}(i);
    cellID_B = matchedCells{2,1}(i);
    
    subplot(1,2,1)

    roiPix = sub2ind(size(cellIDMap),stat{1,1}{cellID_A}.ypix+1,stat{1,1}{cellID_A}.xpix+1);
    cellIDMap(roiPix) = 1;
    B = bwboundaries(cellIDMap);
    h1 = visboundaries(B, 'Color', 'g');
    cellIDMap = zeros(size(Fall.ops.meanImg));
    
    subplot(1,2,2)
    
    roiPix = sub2ind(size(cellIDMap),stat{2,1}{cellID_B}.ypix+1,stat{2,1}{cellID_B}.xpix+1);
    cellIDMap(roiPix) = 1;
    B = bwboundaries(cellIDMap);
    h2 = visboundaries(B, 'Color', 'g');
    cellIDMap = zeros(size(Fall.ops.meanImg));
    
    userInput = [];
    while(isempty(userInput))
        userInput = input('Do these ROIs match? 1 for yes, 0 for no: ');
    end
    if userInput == 0
        unmatchedCells{1,1} = [unmatchedCells{1,1}; cellID_A];
        unmatchedCells{2,1} = [unmatchedCells{2,1}; cellID_B];
    elseif userInput ~= 1
        warning('Answer was not 0 or 1. If you meant to type 0, please type ctrl+c to start over.')
    end
    
    delete([h1 h2])
end

end

%% Match A that are not on B

%which Fall has fewer unmatched cells
size_unmatchedCells1 = size(unmatchedCells{1,1},1);
size_unmatchedCells2 = size(unmatchedCells{2,1},1);
[~,min_ind] = min([size_unmatchedCells1, size_unmatchedCells2]);
[~,max_ind] = max([size_unmatchedCells1, size_unmatchedCells2]);

if min_ind == max_ind
    min_ind = 1;
    max_ind = 2;
end

subplot(1,2,max_ind)

if ~noMatchingROIs
    for i = 1:size(matchedCells{max_ind,1},1)
        cellID = matchedCells{max_ind,1}(i);
        roiPix = sub2ind(size(cellIDMap),stat{max_ind,1}{cellID}.ypix+1,stat{max_ind,1}{cellID}.xpix+1);
        cellIDMap(roiPix) = cellID;
    end
B = bwboundaries(cellIDMap);
h1 = visboundaries(B, 'Color', 'g');
cellIDMap = zeros(size(Fall.ops.meanImg));
end

for i = 1:size(unmatchedCells{max_ind,1},1)
    cellID = unmatchedCells{max_ind,1}(i);
    roiPix = sub2ind(size(cellIDMap),stat{max_ind,1}{cellID}.ypix+1,stat{max_ind,1}{cellID}.xpix+1);
    cellIDMap(roiPix) = cellID;
end
B = bwboundaries(cellIDMap);
h2 = visboundaries(B, 'Color', 'r');
%cellIDMap = zeros(size(Fall.ops.meanImg));

subplot(1,2,min_ind)
cellIDMap2 = zeros(size(Fall.ops.meanImg));
 
if ~noMatchingROIs
    for i = 1:size(matchedCells{min_ind,1},1)
        cellID = matchedCells{min_ind,1}(i);
        roiPix = sub2ind(size(cellIDMap2),stat{min_ind,1}{cellID}.ypix+1,stat{min_ind,1}{cellID}.xpix+1);
        cellIDMap2(roiPix) = cellID;
    end

    B = bwboundaries(cellIDMap2);
    h3 = visboundaries(B, 'Color', 'g');
end

cellMatch = [];
matchCount = 1;

for i = 1:size(unmatchedCells{min_ind,1},1)
    cellID = unmatchedCells{min_ind,1}(i);
    cellIDMap2 = zeros(size(Fall.ops.meanImg));
    roiPix = sub2ind(size(cellIDMap2),stat{min_ind,1}{cellID}.ypix+1,stat{min_ind,1}{cellID}.xpix+1);
    cellIDMap2(roiPix) = cellID;
    subplot(1,2,min_ind)
    B = bwboundaries(cellIDMap2);
    h4 = visboundaries(B, 'Color', 'r');
    
    userInput = [];
    while(isempty(userInput))
        userInput = input('Does this ROI have a match? 1 for yes, 0 for no: ');
    end
    if userInput == 0
    elseif userInput == 1
        [x,y] = ginput(1);
        while cellIDMap(round(y),round(x)) == 0
            disp('Try again')
            [x,y] = ginput(1);
        end
        cellMatch(matchCount,max_ind) = cellID;
        cellMatch(matchCount,min_ind) = cellIDMap(round(y),round(x));
        matchCount = matchCount + 1;
        
        subplot(1,2,max_ind)
        h5 = visboundaries(B, 'Color', 'c');

    elseif userInput ~= 1
        warning('Answer was not 0 or 1. Please type ctrl+c to start over.')
    end
    
        delete([h4])

end


disp('Saving data')
manualMatching.matchedCells = matchedCells;
manualMatching.unmatchedCells = unmatchedCells;
manualMatching.newMatches = cellMatch;
roiMatchData.manualMatching = manualMatching;
save(match_filename,'roiMatchData');
disp('All done')
