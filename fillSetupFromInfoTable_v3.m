function [data] = fillSetupFromInfoTable_v3(Info, compiled_blocks_path, stim_protocol, stim_name)
% This function creates a data.mat file with all of the block data for 
% a given experiment specified by Info
% 
% Argument(s): 
%   Info - Info table loaded from Info excel spreadsheet
%   compiled_blocks_path (string) - filepath where the compiled blocks are saved
%   stim_protocol - number corresponding to stim type that will be analyzed
%   stim_name - if this is included, use a specified stim name instead of protocol to choose which files to analyze
%               For example, there could be multiple types of behavior stim compiled with the same stim_protocol
%   
% Returns:
%   data(struct)
% 
% Notes:
%   v3 January 2021 uses blocks saved with compile_blocks_from_info_v2
%
% TODO: allow stim_protocol to take a list of stims instead of just one
% Search 'TODO'

%% Info.mat Table is a variable that stores all recording file information

%  Currently, the column order of Info is:
I = 1; %ignore (0 or 1 allowing user to "hide" some files from analysis)
P = 2; %pathname
M = 3; %mouse name
D = 4; %experiment date
B = 5; %block name
F = 6; %FOV
TS = 7; %Tosca session #
TR = 8; %Tosca run #
AP = 9; %part of the path, folder where fall.mats are stored
RR = 10; %do you have red cells? 0 or 1 %Not used right now
VN = 11; %full voltage recording name (if widefield only)
SN = 12; %type of stim presentation in plain text
SP = 13; %stim protocol (number corresponding to stim)
GT = 14; %f, m, or s depending on GCaMP type
EG = 15; %name of experimental group or condition

%% We will be looking for files that match stim_protocol
%Later we could update this to also look for other parameters

%stim protocol code is:
%noiseburst = 1
%ReceptiveField = 2
%FM sweep = 3
%SAM = 6
%widefield = 4
%SAM freq = 6
%Behavior = 7 and 8
%Random H20 = 9
%Noiseburst_ITI = 10
%Random air puff = 11
%Spontaneous = 12
%Maryse's Behavior = 13

if nargin < 5
    matchStimName = 0;
else
    matchStimName = 1;
end

code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'Widefield', 'SAM', 'SAM freq' , 'Behavior', 'Behavior',...
    'Random H20', 'Noiseburst ITI', 'Random Air', 'Spontaneous', 'Maryse Behavior'};
disp(['Analyzing ' code{stim_protocol} ' files'])

%% Create data structure

%Add Imaging_set (block #) to Info
IS = size(Info,2) + 1;
Info{1,IS} = 'imaging_set';
for i = 2:size(Info,1)
    if ~ismissing(Info{i,B})
        block_name = Info{i,B}{1};
        Info{i,IS} = str2double(block_name(1,end-2:end));
    end
end

Info(1,:) = []; %Remove header from Info

%Set up data
data = struct;
data.setup = struct;
data.setup.stim_protocol = stim_protocol;
data.setup.Info = Info;
data.setup.compiled_blocks_path = compiled_blocks_path;

%Find table rows that match stim_protocol
setup = data.setup;
stims = [Info{:,SP}]';
if matchStimName
    stimNames = [Info{:,SN}]';
    matching_stims = intersect(find(stims == setup.stim_protocol), find(stimNames == stim_name));
    currentInfo = Info(matching_stims,:);
else
    matching_stims = stims == setup.stim_protocol;
    currentInfo = Info(matching_stims,:);
end

%Remove rows that are set to "Ignore"
ignore = [currentInfo{:,I}]';
currentInfo = currentInfo(ignore == 0,:);

if isempty(currentInfo)
    error('No data found to compile. Check Info sheet.')
end

%% Fill setup and data

%Make data struct with the following format:
%data.([mouseID]).(['ImagingBlock' Imaging_Num]).VARIABLE
%And: data.([mouseID]).parameters

cd(compiled_blocks_path)

%Per each mouse, combine blocks that come from the same FOV
mice = cellstr([currentInfo{:,M}])';
uniqueMice = unique(mice);

for i = 1:length(uniqueMice)
    currentMouse = uniqueMice{i};
    matching_mice = strcmp(currentMouse, mice);
    
    currentInfo_M = currentInfo(matching_mice,:);
    FOVs = [currentInfo_M{:,F}]'; %The FOV corresponds to which data was run together in Suite2p.
    %If left empty it will be NAN and all NANs will be treated as separate FOVs
    %THIS DOES NOT CURRENTLY WORK IF YOU HAVE MIXED NAN AND NUMERICAL FOVs FOR THE SAME ANIMAL
    uniqueFOVs = unique(FOVs);
    
    for j = 1:length(uniqueFOVs)
        currentFOV = uniqueFOVs(j);
        if ismissing(currentFOV)
            currentInfo_F = currentInfo_M(j,:);
        else
            currentInfo_F = currentInfo_M(FOVs == currentFOV,:); 
        end
        
        %Fill setup with structure {Mouse, FOV}
        setup.mousename{i,j}         =   currentMouse;
        setup.FOVs{i,j}              =   currentFOV;
        setup.stim_name{i,j}         =   [currentInfo_F{:,SN}];
        setup.expt_date{i,j}         =   [currentInfo_F{:,D}];
        setup.block_name{i,j}        =   [currentInfo_F{:,B}];
        setup.Imaging_sets{i,j}      =   [currentInfo_F{:,IS}];
        setup.Session{i,j}           =   [currentInfo_F{:,TS}];
        setup.Tosca_Runs{i,j}        =   [currentInfo_F{:,TR}];
                
        %Preallocate variables to record filenames
        block_filenames = cell(1,size(currentInfo_F,1));
        unique_block_names = cell(1,size(currentInfo_F,1));
        
        %Build and record filenames
        for f = 1:size(currentInfo_F,1)

            if ~ismissing(setup.FOVs{i,j})
                FOVtag = strcat('_FOV', setup.FOVs{i,j});
            else
                FOVtag = '';
            end
            
            if ~ismissing(setup.block_name{i,j})
                block_name = setup.block_name{i,j}{f};
                blockTag = ['_Block' block_name(1,end-2:end)];
                uniqueBlockTag = ['Block' num2str([currentInfo_F{f,IS}]) '_'];
            else
                blockTag = '';
                uniqueBlockTag = '';
            end
   
            if ~ismissing(currentInfo_F{f,VN})
                widefieldTag = 'widefield-';
            else
                widefieldTag = '';
            end

            block_filenames{1,f} = strcat('Compiled_', setup.mousename{i,j}, FOVtag, '_', setup.expt_date{i,j}{f}, ...
                '_Session', sprintf('%02d',setup.Session{i,j}(f)), '_Run', sprintf('%02d',setup.Tosca_Runs{i,j}(f)),...
                 blockTag, '_', widefieldTag, setup.stim_name{i,j}{f});
            
            load(block_filenames{1,f})
            disp(block_filenames{1,f})
        
            unique_block_names{1,f} = strcat(uniqueBlockTag,...
                'Session', num2str([currentInfo_F{f,TS}]), '_Run', num2str([currentInfo_F{f,TR}]));

            if isfield(block, 'parameters')
                data.([currentMouse]).parameters = block.parameters; %This will be written over every time
                %This assumes that the block parameters are the same for every stim paradigm, but they might not be
                %For example if some trials are lost. We'll have to fix this at some point.
            else
                block.parameters = nan;
            end
            data.([currentMouse]).([unique_block_names{1,f}]) = block;
            clear('block');
        end
        setup.block_filename{i,j} = [block_filenames{:,:}];
        setup.unique_block_names{i,j} = [string(unique_block_names(:,:))];
    end
end

data.setup = setup;
end