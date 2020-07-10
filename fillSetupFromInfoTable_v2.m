function [data] = fillSetupFromInfoTable_v2(Info, compiled_blocks_path, stim_protocol)
% This function creates a data.mat file with all of the block data for 
% a given experiment specified by Info
% 
% Argument(s): 
%   Info - Info table loaded from Info excel spreadsheet
%   compiled_blocks_path (string) - filepath where the compiled blocks are saved
%   stim_protocol - number corresponding to stim type that will be analyzed
%   
% Returns:
%   data(struct)
% 
% Notes:
%
%
% TODO: allow stim_protocol to take a list of stims instead of just one
% Search 'TODO'

%% Info.mat Table is a variable that stores all recording file information

%  Currently, the column order of Info is:
I = 1; %ignore (0 or 1 allowing user to "hide" some files from analysis)
P = 2; %first part of the path
U = 3; %part of the path, not every user will have this, okay to leave empty
M = 4; %part of the path, no underscores
D = 5; %part of the path, YYYY-MM-DD
B = 6; %part of the path - full block name used for BOT
F = 7; %which data to consider as coming from the same FOV, per mouse
IS = 8; %block or BOT numbers
TS = 9; %Tosca session #
TR = 10; %Tosca run #
AP = 11; %part of the path, folder where fall.mats are stored
FR = 12; %15 or 30, eventually we can detect this automatically
RR = 13; %do you have red cells? 0 or 1 %Not used right now
VR = 14; %voltage recording (0 for widefield, 1 for 2p)
VN = 15; %full voltage recording name (if widefield only)
SN = 16; %type of stim presentation in plain text
SP = 17; %stim protocol (number corresponding to stim)
GT = 18; %f, m, or s depending on GCaMP type
EG = 19; %name of experimental group or condition

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

code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'Widefield', 'SAM', 'SAM freq' , 'Behavior', 'Behavior', 'Random H20', 'Noiseburst ITI', 'Random Air'};
disp(['Analyzing ' code{stim_protocol} ' files'])

%% Create data structure

Info(1,:) = []; %Remove header from Info

data = struct;
data.setup = struct;
data.setup.stim_protocol = stim_protocol;
data.setup.Info = Info;
data.setup.compiled_blocks_path = compiled_blocks_path;

%Find table rows that match stim_protocol
setup = data.setup;
stims = [Info{:,SP}]';
matching_stims = stims == setup.stim_protocol;
currentInfo = Info(matching_stims,:);

%Remove rows that are set to "Ignore"
ignore = [currentInfo{:,I}]';
currentInfo = currentInfo(ignore == 0,:);

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
    FOVs = [currentInfo_M{:,F}]'; %This requires every block to have an FOV number, which corresponds
    %to which data was run together in Suite2p. If left empty there will be an error
    uniqueFOVs = unique(FOVs);
    
    for j = 1:length(uniqueFOVs)
        currentFOV = uniqueFOVs(j);
        currentInfo_R = currentInfo_M(FOVs == currentFOV,:); 
        
        try
            gcamp_type = [currentInfo_R{:,GT}];
        catch
            gcamp_type = nan;
        end
        
        try
            expt_group = [currentInfo_R{:,EG}];
        catch
            expt_group = nan;
        end

        %Fill setup with structure {Mouse, FOV}
        setup.mousename{i,j}         =   currentMouse;
        setup.FOVs{i,j}              =   currentFOV;
        setup.gcamp_type{i,j}        =   gcamp_type;
        setup.expt_group{i,j}        =   expt_group;
        setup.expt_date{i,j}         =   [currentInfo_R{:,D}];
        setup.Imaging_sets{i,j}      =   [currentInfo_R{:,IS}];
        setup.Session{i,j}           =   [currentInfo_R{:,TS}];
        setup.Tosca_Runs{i,j}        =   [currentInfo_R{:,TR}];
        setup.FrameRate{i,j}         =   [currentInfo_R{:,FR}] 
                
        block_filenames = cell(1,size(currentInfo_R,1));
        unique_block_names = cell(1,size(currentInfo_R,1));
        path_name =currentInfo_R{:,P};
        
        for r = 1:size(currentInfo_R,1)
            
            Block_number = sprintf('%03d',currentInfo_R{r,IS});
  
            if currentInfo_R{r,VR} == 0
                widefieldTag = 'widefield-';
                setup.VRname = [currentInfo_R{:,VN}];
                setup.BOTname = [currentInfo_R{:,B}];
                setup.VRpath = strcat(path_name, '/', currentMouse, '/', setup.expt_date{i,j}, '/', setup.VRname); 
                setup.BOTpath = strcat(path_name, '/', currentMouse, '/', setup.expt_date{i,j}, '/', setup.BOTname); 
                
                
            else
                widefieldTag = '';
            end
    
            block_filenames{1,r} = strcat('Compiled_', currentMouse, '_', [currentInfo_R{r,D}], ...
            '_Block_', Block_number, '_Session_', num2str([currentInfo_R{r,TS}]), ...
            '_Run_', num2str([currentInfo_R{r,TR}]), '_', widefieldTag, [currentInfo_R{r,SN}]);
            
            load(block_filenames{1,r})
        
            unique_block_names{1,r} = strcat('Block', num2str([currentInfo_R{r,IS}]),...
                '_Session', num2str([currentInfo_R{r,TS}]), '_Run', num2str([currentInfo_R{r,TR}]));

            if isfield(block, 'parameters')
                data.([currentMouse]).parameters = block.parameters; %This will be written over every time
                %This assumes that the block parameters are the same for every stim paradigm, but they might not be
                %For example if some trials are lost. We'll have to fix this at some point.
            else
                block.parameters = nan;
            end
            data.([currentMouse]).([unique_block_names{1,r}]) = block;
            clear('block');
        end
        setup.block_filename{i,j} = [block_filenames{:,:}];
        setup.unique_block_names{i,j} = [string(unique_block_names(:,:))];
    end
end

data.setup = setup;
end