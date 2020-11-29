%Script to extract y from t-series

mouse_dir = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\NxDE072420M2\2020-11-19';
cd(mouse_dir)
currentFolders = dir;

%% Get list of TSeries folders

isDir = [currentFolders.isdir]';
isTseries = zeros(size(currentFolders,1),1);
for i = 1:size(currentFolders,1)
    folderNames = {currentFolders.name}';
    currentFolderName = folderNames{i,1};
    if length(currentFolderName) < 7; continue; end
    isTseries(i,1) = input(['Include ' currentFolderName ' in TSeries? 1 for yes, 0 for no: ']);
     %isTseries(i,1) = strcmp(currentFolderName(1:7),'TSeries') || strcmp(currentFolderName(1:7),'Tseries');
end
isTseries = logical(isTseries);

tseriesFolders = folderNames;
tseriesFolders((isDir + isTseries) < 2,:) = [];

%% Go into each folder and extract TSeries

data = struct;
data.Directory = mouse_dir;
data.Folders = tseriesFolders;
data.xMats = {};

for t = 1:length(tseriesFolders)
    currentFolder = tseriesFolders{t,1};
    cd([mouse_dir '\' currentFolder])
    tseriesFiles = dir('*.tif');
    tseriesFileNames = {tseriesFiles.name}';
    fileNameSplit = {};
    for c = 1:length(tseriesFileNames)
        fileNameSplit = [fileNameSplit; strsplit(tseriesFileNames{c,1},'_')];
    end
   
    %Make xMat for nFrames x [Cycle, Image]
    nFrames = size(fileNameSplit,1)/2; %Divide by 2 because there are 2 ch
    nParameters = size(fileNameSplit,2);
    ind_image = nParameters;
    ind_ch = nParameters - 1;
    ind_cycle = nParameters - 2;
    
    %extract parameters
    cell_image = [fileNameSplit(:,ind_image)];
    cell_ch = [fileNameSplit(:,ind_ch)];
    cell_cycle = [fileNameSplit(:,ind_cycle)];
    
    %transform parameters to doubles
    for c = 1:length(cell_image)
        cell_image{c,1} = cell_image{c,1}(1:6);
        cell_ch{c,1} = cell_ch{c,1}(end);
        cell_cycle{c,1} = cell_cycle{c,1}(end-4:end);
    end
    dbl_image = str2double(cell_image)';
    dbl_ch = str2double(cell_ch)';
    dbl_cycle = str2double(cell_cycle)';
    
    %only keep values from channel 1
    dbl_image1 = dbl_image(dbl_ch == 1);
    dbl_cycle1 = dbl_cycle(dbl_ch == 1);
    
    xMat = [dbl_cycle1; dbl_image1];
    
    data.xMats{t,1} = xMat;
end

cd(mouse_dir)
currentTime = clock;
save(['TSeries_Data_' date '-' num2str(currentTime(4)) num2str(currentTime(5))], 'data')