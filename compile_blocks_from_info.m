%% compile_blocks_from_info
%
%  This script saves a compiled 'block.mat' file for each row in Info.
%  A compiled block contains Suite2p, Bruker, and Tosca data.
%
%  The Info spreadsheet contains the information required to access the
%  data by specifying the filepaths to each.
%
%  If data is missing, the script will still run in order to allow for the
%  user to look at components of the data individually.
%
%  Uses the functions:
%  - importfile
%  - define_behavior_singleblock
%  - FreqDisc_Behavior_singleblock
%  - define_sound_singleblock
%  - define_loco_singleblock
%  - define_suite2p_singleblock
%  - align_to_stim
%  - visualize_block
%    
%  Use visualize_block to preview the block contents once they've been compiled.
%
%  TAKESIAN LAB - March 2020

%% Load Info.mat and change user-specific options

visualize = 0; %1 to plot figures of the block immediately, 0 to skip
recompile = 1; %1 to save over previously compiled blocks, 0 to skip
checkOps = 0; %1 to check Fall.ops against user-specified ops.mat file

%% set up values for 'align to stim'

% How many seconds of baseline?
constant.baseline_length = 0.5;

% How many seconds after stim should we look at?
constant.after_stim = 3;

% Define (in seconds) where to look for the response peak?
constant.response_window = 2;

% define where to look for locomotor responses, in sec?
constant.locowindow = 2.5;

%minimum amout of time (sec) that mouse is moving to be considered active
constant.locoThresh = 0.8;

% Define the neuropil coefficient
% TODO: automatically grab this from Suite2p
constant.neucoeff = 0.7;

%% 
PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
        info_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study';
        save_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\CompiledBlocks';
        %info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
        %save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\Compiled Blocks';
        info_filename = 'Info_NxDC030220F2';
        ops_filename = 'Maryse_ops2.mat';
        
    case 'TAKESIANLAB2P' %2P computer
        info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
        save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\Compiled Blocks';
        info_filename = 'Info';    
        
    case 'RD0332' %Carolyn
        info_path = 'D:\2P analysis\2P local data\Carolyn';
        save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Compiled Blocks';
        info_filename = 'Info_widefield_YD111219M1';
        
         case 'RD-6-TAK2' %Esther's computer
        info_path = 'Z:\Carolyn\2P Imaging data\SSRI study with Jacob';
        save_path = 'Z:\Carolyn\2P Imaging data\SSRI study with Jacob\Compiled Blocks';
        info_filename = 'Info_widefield_YD111219M1';
        
    case 'RD0386' %Wisam
        % INSERT PATHS HERE
        info_filename = 'Info';
    otherwise
        disp('Computer does not match known users')
        return
end

cd(info_path)
Info = importfile(info_filename);

if checkOps
    load(ops_filename);
    user_ops = ops;
    clear('ops');
    user_ops.checkOps = 1;
else
    user_ops.checkOps = 0;
end

%% Compile all blocks unless they are set to "Ignore"
%  No need to change any variables below this point

%Add extra empty columns when updates might be ahead of people's Info sheets
lastColumn = 19; %WARNING: Magic number
if size(Info,2) < lastColumn
    for i = (size(Info,2) + 1):lastColumn
        Info{1,i} = 'Expand Columns';
    end
    warning('gcamp_type and expt_group columns missing from Info')
end
        
%Remove header from Info
Info(1,:) = [];

%Remove rows that are set to "Ignore"
ignore = [Info{:,1}]';
currentInfo = Info(ignore == 0,:);

%Loop through all remaining rows
for i = 1:size(currentInfo,1)

    %Create setup variable that will contain all the necessary information about the block
    setup = struct;
    setup.constant = constant;
    setup.Info              =   Info;                   %Record Info for records only
    setup.pathname          =   [currentInfo{i,2}];     %first part of the path
    setup.username          =   [currentInfo{i,3}];     %part of the path, not every user will have this, okay to leave empty
    setup.mousename         =   [currentInfo{i,4}];     %part of the path, no underscores
    setup.expt_date         =   [currentInfo{i,5}];     %part of the path, YYYY-MM-DD
    setup.block_name        =   [currentInfo{i,6}];     %part of the path - full block name used for BOT
    setup.FOV               =   [currentInfo{i,7}];     %which data to consider as coming from the same field of view, per mouse
    setup.imaging_set       =   [currentInfo{i,8}];     %block or BOT numbers
    setup.Tosca_session     =   [currentInfo{i,9}];     %Tosca session
    setup.Tosca_run         =   [currentInfo{i,10}];    %Tosca run
    setup.analysis_name     =   [currentInfo{i,11}];    %part of the path, folder where fall.mats are stored
    setup.framerate         =   [currentInfo{i,12}];    %15 or 30, eventually we can detect this automatically
    setup.run_redcell       =   [currentInfo{i,13}];    %do you have red cells? 0 or 1
    setup.voltage_recording =   [currentInfo{i,14}];    %0 for widefield, 1 for 2p
    setup.VR_name           =   [currentInfo{i,15}];    %full voltage recording name (if widefield only)
    setup.stim_name         =   [currentInfo{i,16}];    %type of stim presentation in plain text
    setup.stim_protocol     =   [currentInfo{i,17}];    %number corresponding to stim protocol
    setup.gcamp_type        =   [currentInfo{i,18}];    %f, m, or s depending on GCaMP type
    setup.expt_group        =   [currentInfo{i,19}];    %name of experimental group or condition
    
    Block_number = sprintf('%03d',setup.imaging_set);
    
    if setup.voltage_recording == 0
        widefieldTag = 'widefield-';
    else
        widefieldTag = '';
    end
    
    setup.block_filename = strcat('Compiled_', setup.mousename, '_', setup.expt_date, ...
        '_Block_', Block_number, '_Session_', num2str(setup.Tosca_session), ...
        '_Run_', num2str(setup.Tosca_run), '_', widefieldTag, setup.stim_name);
    setup.block_supname = strcat(setup.mousename, ' ', setup.expt_date, ...
        ' ', setup.stim_name, ' ', Block_number);
    
    %Skip files that have previously been compiled
    if ~recompile
        cd(save_path)
        if isfile(strcat(setup.block_filename, '.mat'))
            disp(setup.block_filename);
            disp('Skipping (already compiled)');
            continue
        end
    end
    
    %Not every user has a username folder, allow for this column to be empty
    if ~ismissing(setup.username)
        usernameSlash = strcat(setup.username, '/');
    else
        usernameSlash = '';
    end
    
    %Establish and test paths, allowing for paths to be missing
    %TOSCA PATH
    if ismissing(setup.Tosca_session)
        setup.Tosca_path = nan;
    else
        setup.Tosca_path = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/Tosca_', setup.mousename, {'/Session '}, num2str(setup.Tosca_session));
        if ~isfolder(setup.Tosca_path)
            disp(setup.block_filename)
            error('Your Tosca path is incorrect.')
        end
    end
    
    %BLOCK PATH
    if ismissing(setup.block_name)
        setup.block_path = nan;
    else
        setup.block_path   = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.expt_date, '/', setup.block_name);
        if ~isfolder(setup.block_path)
            disp(setup.block_filename)
            error('Your block path is incorrect.')
        end
    end
    
    %VR PATH
    if ismissing(setup.VR_name)
        setup.VR_path = nan;
    else
        setup.VR_path  = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.expt_date, '/', setup.VR_name);
        if ~isfolder(setup.VR_path)
            disp(setup.block_filename)
            error('Your voltage recording path is incorrect.')
        end
    end
    
    %SUITE2P PATH
    if ismissing(setup.analysis_name)
        setup.suite2p_path = nan;
    else
        setup.suite2p_path = strcat(setup.pathname, '/', usernameSlash, '/', setup.analysis_name);
        if ~isfolder(setup.suite2p_path)
            disp(setup.block_filename)
            error('Your Suite2p analysis path is incorrect.')
        end
    end
   
    
    %% COMPILE BLOCK
    disp('Processing...');
    disp(setup.block_filename);

    block = struct;
    block.setup = setup;

    % Behavior, locomotion, sound, and Suite2p data

    %pull out the Tosca-derived, behaviorally relevant data
    [block] = define_behavior_singleblock(block);
    [block] = FreqDisc_Behavior_singleblock(block);

    %pull out the Bruker-derived timestamps from BOTs and Voltage Recordings
    [block] = define_sound_singleblock(block);

    %determine which trials are considered "active (locomotor)"
    % This might not be necessary to do here, but leaving in for now.
    [block] = define_loco_singleblock(block);

    %pull out block-specific data from Fall.mat
    [block] = define_suite2p_singleblock(block, user_ops);
    
    %find the stim-aligned traces
    [block] = align_to_stim(block);
    
    %Optionally visually check block
    if visualize == 1
        visualize_block(block);
    end
    
    %% Save block
    disp('Saving...');
    cd(save_path)
    save(setup.block_filename, 'block');
end       

disp('Finished compiling all blocks.');