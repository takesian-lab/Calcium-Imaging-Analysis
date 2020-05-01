function [data, setup] = fillSetupFromInfoTable_wisam(setup, Info, compiled_blocks_path)
% [data, setup] = fillSetupFromInfoTable_wisam(setup, Info, compiled_blocks_path)
%
% THIS DOCUMENTATION IS A WORK IN PROGRESS:
%
%   This function...
%
% Details:
%   More info...
% 
% Arguments: 
%   setup (struct)
%   Info (excel sheet)
%   compiled_blocks_path (string)
%
% Returns:
%   data (struct)
%   setup (struct)
%
% Search 'TODO'

%% Info.mat Table is a variable that stores all recording file information
% TODO: Is all of this strong typing necessary?
% TODO: Do we need to translate from numeric values?
% Currently, the column order of Info is:
I = 1;      % ignore (0 or 1 allowing user to "hide" some files from analysis)
P = 2;      % first part of the path
U = 3;      % part of the path, not every user will have this, okay to leave empty
M = 4;      % part of the path, no underscores
D = 5;      % part of the path, YYYY-MM-DD
B = 6;      % part of the path - full block name used for BOT
F = 7;      % which data to consider as coming from the same FOV, per mouse
IS = 8;     % block or BOT numbers
TS = 9;     % Tosca session #
TR = 10;    % Tosca run #
AP = 11;    % part of the path, folder where fall.mats are stored
FR = 12;    % 15 or 30, eventually we can detect this automatically
RR = 13;    % do you have red cells? 0 or 1 %Not used right now
VR = 14;    % voltage recording (0 for widefield, 1 for 2p)
VN = 15;    % full voltage recording name (if widefield only)
SN = 16;    % type of stim presentation in plain text
SP = 17;    % stim protocol (number corresponding to stim)

%% Look for files that match stim_protocol
% Later we could update this to also look for other parameters

% stim protocol codes are:
% noiseburst = 1
% ReceptiveField = 2
% FM sweep = 3
% SAM = 5
% widefield = 4
% SAM freq = 6
code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'SAM', 'Widefield', 'SAM freq'};
disp(['Analyzing ' code{setup.stim_protocol} ' files'])

% Remove header from Info
Info(1,:) = [];

% Find table rows that match stim_protocol
stims = [Info{:,SP}]';
matching_stims = stims == setup.stim_protocol;
currentInfo = Info(matching_stims,:);

% Remove rows that are set to "Ignore"
ignore = [currentInfo{:,I}]';
currentInfo = currentInfo(ignore == 0,:);

%% Fill setup and data

% Make data struct with the following format:
% data.([mouseID]).(['ImagingBlock' Imaging_Num]).VARIABLE
% And: data.([mouseID]).parameters

cd(compiled_blocks_path)
data = struct;

% Per each mouse, combine blocks that come from the same FOV
mice = cellstr([currentInfo{:,M}])';
uniqueMice = unique(mice);

% Loop over all unique mice in the excel sheet
for i = 1:length(uniqueMice)
    
    % How many unique mice are in the excel sheet
    currentMouse = uniqueMice{i};
    % TODO: Check how this works with a sheet with multiple mice / blocks, etc
    matching_mice = strcmp(currentMouse, mice);
    
    % Collect data across excel sheet rows cooresponding to a particular
    % mouse (by ID)
    % TODO: Check how this works with a sheet with multiple mice / blocks, etc
    currentInfo_M = currentInfo(matching_mice,:);
    FOVs = [currentInfo_M{:,F}]'; 
    % This requires every block to have an FOV number, which corresponds
    % to which data was run together in Suite2p. If left empty there will be an error
    uniqueFOVs = unique(FOVs);
    
    % For a given unique mouse in the excel sheet:
    % Loop over all unique fields of view (FOV)
    for j = 1:length(uniqueFOVs)
        
        currentFOV = uniqueFOVs(j);
        % Grab all the rows in the excel sheet cooresponding to the current
        % mouse and the current FOV
        currentInfo_R = currentInfo_M(FOVs == currentFOV,:); 
        
        % Fill setup with structure {Mouse, FOV}
        % Current Mouse ID
        setup.mousename{i,j}         =   currentMouse;
        % Current FOV for a given mouse ID
        setup.FOVs{i,j}              =   currentFOV;
        % Experiment Date
        setup.expt_date{i,j}         =   [currentInfo_R{:,D}];
        % block or BOT numbers
        setup.Imaging_sets{i,j}      =   [currentInfo_R{:,IS}];
        % Tosca Session(s)
        setup.Session{i,j}           =   [currentInfo_R{:,TS}];
        % Tosca Run(s)
        setup.Tosca_Runs{i,j}        =   [currentInfo_R{:,TR}];
        
        % Pre-allocate cell array for block filenames
        block_filenames = cell(1,size(currentInfo_R,1));
        % Pre-allocate cell array for block filenames
        unique_block_names = cell(1,size(currentInfo_R,1));
        
        % For a given unique mouse and a unique field of view (FOV):
        % Loop over the blocks
        for r = 1:size(currentInfo_R,1)
            
            % String for current block number
            Block_number = sprintf('%03d',currentInfo_R{r,IS});
  
            % voltage recording (0 for widefield, 1 for 2p)
            if currentInfo_R{r,VR} == 0
                widefieldTag = 'widefield-';
            else
                widefieldTag = '';
            end
    
            % Create the full block file name
            block_filenames{1,r} = strcat('Compiled_', currentMouse, '_', [currentInfo_R{r,D}], ...
            '_Block_', Block_number, '_Session_', num2str([currentInfo_R{r,TS}]), ...
            '_Run_', num2str([currentInfo_R{r,TR}]), '_', widefieldTag, [currentInfo_R{r,SN}]);
            
            % TODO: Why is this being loaded?
            % I thought the point of this function is to generate this
            % compiled file
            load(block_filenames{1,r})
        
            % TODO: Why are we doing this?
            % This appears to be generating a unique block name string 
            % 'Block#_Session#_Run#'
            unique_block_names{1,r} = strcat('Block', num2str([currentInfo_R{r,IS}]),...
                '_Session', num2str([currentInfo_R{r,TS}]), '_Run', num2str([currentInfo_R{r,TR}]));

            % TODO: 'block' is undefined? Why?
            if isfield(block, 'parameters')
                % TODO: FIX THIS
                data.([currentMouse]).parameters = block.parameters; % This will be written over every time
                % This assumes that the block parameters are the same for every stim paradigm, but they might not be
                % For example if some trials are lost. We'll have to fix this at some point.
            else block.parameters = nan;
            end
            
            data.([currentMouse]).([unique_block_names{1,r}]) = block;
            % TODO: Why is this happening?
            clear('block');
        end
        
        setup.block_filename{i,j} = [block_filenames{:,:}];
        setup.unique_block_names{i,j} = [string(unique_block_names(:,:))];
    end
end

end