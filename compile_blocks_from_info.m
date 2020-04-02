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

%% Load Info.mat

visualize = 0; %1 to plot figures of the block immediately, 0 to skip
recompile = 1; %1 to save over previously compiled blocks, 0 to skip

info_path = 'D:/2P analysis/2P local data/Carolyn';
save_path = 'D:/2P analysis/2P local data/Carolyn/analyzed/Daily Imaging';
cd(info_path)
Info = importfile('Info');

%% Compile all blocks unless they are set to "Ignore"

%Remove header from Info
Info(1,:) = [];

%Remove rows that are set to "Ignore"
ignore = [Info{:,1}]';
currentInfo = Info(ignore == 0,:);

%Loop through all remaining rows
for i = 1:size(currentInfo,1)

    %Create setup variable that will contain all the necessary information about the block
    setup = struct;
    setup.Info              =   Info;                   %Record Info for records only
    setup.pathname          =   [currentInfo{i,2}];     %first part of the path
    setup.username          =   [currentInfo{i,3}];     %part of the path, not every user will have this, okay to leave empty
    setup.mousename         =   [currentInfo{i,4}];     %part of the path, no underscores
    setup.expt_date         =   [currentInfo{i,5}];     %part of the path, YYYY-MM-DD
    setup.block_name        =   [currentInfo{i,6}];     %part of the path - full block name used for BOT
    setup.ROI               =   [currentInfo{i,7}];     %which data to consider as coming from the same ROI, per mouse
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

    Block_number = sprintf('%03d',setup.imaging_set);
    
%     setup.block_filename = strcat('Compiled_', setup.mousename, '_', setup.expt_date, ...
%         '_Block_', num2str('%03d',setup.imaging_set), '_Session_', num2str(setup.Tosca_session), ...
%         '_Run_', num2str(setup.Tosca_run), '_', setup.stim_name);
    
     setup.block_filename = strcat('Compiled_', setup.mousename, '_', setup.expt_date, ...
        '_Block_', Block_number, '_Session_', num2str(setup.Tosca_session), ...
        '_Run_', num2str(setup.Tosca_run), '_', setup.stim_name);
    
    %Skip files that have previously been compiled
    if ~recompile
        cd(save_path)
        if isfile(strcat(setup.block_filename, '.mat'))
            continue
        end
    end
    
    %Not every user has a username folder, allow for this column to be empty
    if ~isnan(setup.username)
        usernameSlash = strcat(setup.username, '/');
    else
        usernameSlash = '';
    end
    
    setup.block_path   = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.expt_date, '/', setup.block_name);
    setup.suite2p_path = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.analysis_name);
    setup.Tosca_path   = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/Tosca_', setup.mousename, {'/Session '}, num2str(setup.Tosca_session));
    if setup.voltage_recording == 0
        setup.VR_path  = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.expt_date, '/', setup.VR_name);
    end
    
    %Test paths prior to starting to compile
    try
        cd(setup.block_path)
        cd(setup.suite2p_path)
        cd(setup.Tosca_path)
        if setup.voltage_recording == 0
            cd(setup.VR_path)
        end
    catch
        disp(setup.block_filename)
        error('One of your paths is incorrect.')
    end
    
    %% COMPILE BLOCK
    disp('Processing...');
    disp(setup.block_name);

    block = struct;
    block.setup = setup;

    % Behavior, locomotion, sound, and Suite2p data

    %pull out the Tosca-derived, behaviorally relevant data
    [block] = define_behavior_singleblock(block);

    %pull out the Bruker-derived timestamps from BOTs and Voltage Recordings
    [block] = define_sound_singleblock(block);

    %determine which trials are considered "active (locomotor)"
    % This might not be necessary to do here, but leaving in for now.
    [block] = define_loco_singleblock(block);

    %pull out block-specific data from Fall.mat
    [block] = define_suite2p_singleblock(block);
    
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