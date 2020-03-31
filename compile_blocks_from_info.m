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

info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
cd(info_path)
load('Info.mat')

%% Compile all blocks unless they are set to "Ignore"

%Remove rows that are set to "Ignore"
ignore = [Info{:,1}]';
currentInfo = Info(ignore == 0,:);

%Loop through all remaining rows
for i = 1:size(currentInfo,1)

    %Create setup variable that will contain all the necessary information about the block
    setup = struct;
    setup.Info              =   Info;                   %Record Info for records only
    setup.username          =   [currentInfo{i,2}];     %part of the path, not every user will have this, okay to leave empty
    setup.mousename         =   [currentInfo{i,3}];     %part of the path, no underscores
    setup.pathname          =   [currentInfo{i,4}];     %first part of the path
    setup.expt_date         =   [currentInfo{i,5}];     %part of the path, YYYY-MM-DD
    setup.ROI               =   [currentInfo{i,6}];     %which data to consider as coming from the same ROI, per mouse
    setup.imaging_set       =   [currentInfo{i,7}];     %block or BOT numbers
    setup.Tosca_session     =   [currentInfo{i,8}];     %Tosca session
    setup.Tosca_run         =   [currentInfo{i,9}];     %Tosca run
    setup.analysis_name     =   [currentInfo{i,10}];    %part of the path, folder where fall.mats are stored
    setup.framerate         =   [currentInfo{i,11}];    %15 or 30, eventually we can detect this automatically
    setup.run_redcell       =   [currentInfo{i,12}];    %do you have red cells? 0 or 1
    setup.voltage_recording =   [currentInfo{i,13}];    %0 for widefield, 1 for 2p
    setup.stim_name         =   [currentInfo{i,14}];    %type of stim presentation in plain text
    setup.stim_protocol     =   [currentInfo{i,15}];    %number corresponding to stim protocol

    setup.block_name = strcat('Compiled_', setup.mousename, '_', setup.expt_date, ...
        '_Block_', num2str(setup.imaging_set), '_Session_', num2str(setup.Tosca_session), ...
        '_Run_', num2str(setup.Tosca_run), '_', setup.stim_name);
    
    %Not every user has a username folder, allow for this column to be empty
    if ~isempty(setup.username)
        usernameSlash = strcat(setup.username, '/');
    else
        usernameSlash = '';
    end
    
    setup.suite2p_path = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/', setup.analysis_name);
    setup.Tosca_path   = strcat(setup.pathname, '/', usernameSlash, setup.mousename, '/Tosca_', setup.mousename);
    
    block = compile_block(setup);

    %Optionally visually check block
    if visualize == 1
        visualize_block(block);
    end
    
    cd(save_path)
    save(setup.block_name, 'block');
    
end        