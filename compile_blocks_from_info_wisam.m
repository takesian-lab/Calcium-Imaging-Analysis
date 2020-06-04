%% compile_blocks_from_info
%
%  Save one compiled 'block.mat' file for each block listed in info.mat
%  A compiled block contains the 2p and Tosca data for that recording session
%
%  Info.mat is a variable that stores all recording file information
%  where each row of Info corresponds to a single recording block
%
%  Goal: If data is missing, the script will still run in order to allow for the
%  user to look at Tosca data separate from 2p data or vice versa.
%
%  Uses the function compile_block and verify_block
%
%  Use visualize_block to preview the block contents prior to analysis
%
%  Maryse Thomas - March 2020

%% Wisam's Notes

% This script: 
% compiles the behavior, locomotion, sound, and Suite2p data for every
% block
% 
% Search 'TODO'

%% CLEAN AND CLEAR

close all
clear 
clc

%% Load Info.mat

% 1 to plot figures of the block immediately, 0 to skip
visualize = 0; 
% 1 to save over previously compiled blocks, 0 to skip
recompile = 1; 

% % % % % % PC_name = getenv('computername');
% % % % % % 
% % % % % % switch PC_name
% % % % % %     case 'RD0366' %Maryse
% % % % % %         info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
% % % % % %         save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocks';
% % % % % %         info_filename = 'Info';
% % % % % %     case 'RD0332' %Carolyn
% % % % % %         info_path = 'D:\2P analysis\2P local data\Carolyn';
% % % % % %         save_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Daily Imaging';
% % % % % %         info_filename = 'Info';
% % % % % %     case 'RD0386' %Wisam
% % % % % %         % INSERT PATHS HERE
% % % % % %         info_filename = 'Info';
% % % % % %     otherwise
% % % % % %         disp('Computer does not match known users')
% % % % % % end

% Path to Info file (an excel sheet)
info_path = '/Users/wisamreid/Documents/School/Research (Harvard)/Takesian/2P/Thy1 Experiments/YD111219F3-2P-noisebursts';
% Path to Compiled file path (a .mat file) 
save_path = '/Users/wisamreid/Documents/School/Research (Harvard)/Takesian/2P/Thy1 Experiments/YD111219F3-2P-noisebursts/Compiled';
info_filename = 'Info';

% TODO: Do we need to cd around?
cd(info_path)
% Import info file (an excel sheet)
Info = importfile(info_filename);

%% Compile all blocks unless they are set to "Ignore"

% TODO: What if we keep these headers and use them to create our own
% namespaces?
% Remove header from Info
Info(1,:) = [];

% Remove rows that are set to "Ignore"
ignore = [Info{:,1}]';
currentInfo = Info(ignore == 0,:);

% Loop through all remaining rows in the excel sheet
% For every row: 
% 1) create a 'setup' struct containing all the necessary information about the block
% 2) setup the file paths 
% 3) and compile the block data
for i = 1:size(currentInfo,1)

    % Create setup variable that will contain all the necessary information about the block
    setup = struct;
    setup.Info              =   Info;                   % Record Info for records only
    setup.pathname          =   [currentInfo{i,2}];     % First part of the path
    setup.username          =   [currentInfo{i,3}];     % Part of the path, not every user will have this, okay to leave empty
    setup.mousename         =   [currentInfo{i,4}];     % Part of the path, no underscores
    setup.expt_date         =   [currentInfo{i,5}];     % Part of the path, YYYY-MM-DD
    setup.block_name        =   [currentInfo{i,6}];     % Part of the path - full block name used for BOT
    setup.FOV               =   [currentInfo{i,7}];     % Which data to consider as coming from the same field of view, per mouse
    setup.imaging_set       =   [currentInfo{i,8}];     % Block or BOT numbers
    setup.Tosca_session     =   [currentInfo{i,9}];     % Tosca session
    setup.Tosca_run         =   [currentInfo{i,10}];    % Tosca run
    setup.analysis_name     =   [currentInfo{i,11}];    % Part of the path, folder where fall.mats are stored
    setup.framerate         =   [currentInfo{i,12}];    % 15 or 30, eventually we can detect this automatically
    setup.run_redcell       =   [currentInfo{i,13}];    % Do you have red cells? 0 or 1
    setup.voltage_recording =   [currentInfo{i,14}];    % 0 for widefield, 1 for 2p
    setup.VR_name           =   [currentInfo{i,15}];    % Full voltage recording name (if widefield only)
    setup.stim_name         =   [currentInfo{i,16}];    % Type of stim presentation in plain text
    setup.stim_protocol     =   [currentInfo{i,17}];    % Number corresponding to stim protocol

    Block_number = sprintf('%03d',setup.imaging_set);
    
    % Widefield recordings do not require a voltage recording (VR)
    if setup.voltage_recording == 0
        widefieldTag = 'widefield-';
    else
        widefieldTag = '';
    end
    
    % Build the full filename for the block
    setup.block_filename = strcat('Compiled_', setup.mousename, '_', setup.expt_date, ...
        '_Block_', Block_number, '_Session_', num2str(setup.Tosca_session), ...
        '_Run_', num2str(setup.Tosca_run), '_', widefieldTag, setup.stim_name);
    
    % Skip files that have previously been compiled
    % TODO: check this logic
    if ~recompile
        cd(save_path)
        if isfile(strcat(setup.block_filename, '.mat'))
            continue
        end
    end
    
    % Not every user has a username folder, allow for this column to be empty
    if ~ismissing(setup.username)
        usernameSlash = strcat(setup.username, '/');
    else
        usernameSlash = '';
    end
    
    % TODO: Clarify these notes.  Is this the Tosca path?
    % Establish and test paths, allowing for paths to be missing
    % TOSCA PATH
    if ismissing(setup.Tosca_session)
        setup.Tosca_path = nan;
    else
        setup.Tosca_path = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/Tosca_', setup.mousename, {'/Session '}, num2str(setup.Tosca_session));
        try
            cd(setup.Tosca_path)
        catch
            disp(setup.block_filename)
            error('Your Tosca path is incorrect.')
        end
    end
    
    % BLOCK PATH
    % These are the raw tiff files from prairie view
    if ismissing(setup.block_name)
        setup.block_path = nan;
    else
        setup.block_path   = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.expt_date, '/', setup.block_name);
        try
            cd(setup.block_path)
        catch
            disp(setup.block_filename)
            error('Your block path is incorrect.')
        end
    end
    
    % Voltage recording (VR) PATH
    if ismissing(setup.VR_name)
        setup.VR_path = nan;
    else
        setup.VR_path  = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.expt_date, '/', setup.VR_name);
        try
            cd(setup.VR_path)
        catch
            disp(setup.block_filename)
            error('Your voltage recording path is incorrect.')
        end
    end
    
    % SUITE2P PATH 
    % This is looking for the 'Fall.mat' file generated by Suite2P
    if ismissing(setup.analysis_name)
        setup.suite2p_path = nan;
    else
        setup.suite2p_path = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.analysis_name);
        try
            cd(setup.suite2p_path)
        catch
            disp(setup.block_filename)
            error('Your Suite2p analysis path is incorrect.')
        end
    end
   
    
    %% COMPILE BLOCK
    % Behavior, locomotion, sound, and Suite2p data
    
    disp('Processing...');
    disp(setup.block_filename);

    % Setup a block struct
    block = struct;
    % Store the setup for the given block
    block.setup = setup;

    % Pull out the Tosca-derived, behaviorally relevant data
    % [block] = define_behavior_singleblock(block);
    [block] = define_behavior_singleblock_wisam(block);
    % TODO: What is this?
    % [block] = FreqDisc_Behavior_singleblock(block);
    [block] = FreqDisc_Behavior_singleblock_wisam(block);

    % pull out the Bruker-derived timestamps from BOTs and Voltage Recordings
    % [block] = define_sound_singleblock(block);
    [block] = define_sound_singleblock_wisam(block);

    % determine which trials are considered "active (locomotor)"
    % This might not be necessary to do here, but leaving in for now.
    % TODO: Is this needed?
    [block] = define_loco_singleblock(block);

    % pull out block-specific data from Fall.mat
    [block] = define_suite2p_singleblock(block);
    
    %find the stim-aligned traces
%     [block] = align_to_stim(block);
    [block] = align_to_stim_wisam(block);
    
    % Optionally visually check block
    if visualize == 1
        visualize_block(block);
    end
    
    %% Save block
    disp('Saving...');
    cd(save_path)
    save(setup.block_filename, 'block');
    
end       

disp('Finished compiling all blocks.');