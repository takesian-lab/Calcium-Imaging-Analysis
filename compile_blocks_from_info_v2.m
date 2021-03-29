%% compile_blocks_from_info_v2
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
%    
%  Use visualize_block and visualize_cell to preview the block contents once they've been compiled.
%
%  TAKESIAN LAB - March 2020
%  v2 December 2020 saves blocks with new naming system

%% Load Info.mat and change user-specific options

recompile = 0; %1 to save over previously compiled blocks, 0 to skip
checkOps = 0; %1 to check Fall.ops against user-specified ops.mat file

%% set up values for 'align to stim'
% Ndnf vs. Vip project: 0.5, 2.5, 1.5, 1.5, 0.8, 0.7

% How many seconds of baseline?
constant.baseline_length = 0.5;

% How many seconds after stim should we look at?
% Can be overwritten by column in info
constant.after_stim = 2.5;

% Define (in seconds) where to look for the response peak?
constant.response_window = 1.5;

% define where to look for locomotor responses, in sec?
constant.locowindow = 1.5;

%minimum amout of time (sec) that mouse is moving to be considered active
constant.locoThresh = 0.8;

% Define the neuropil coefficient
% This will be checked against Suite2p value and a warning will appear if they do not match
constant.neucoeff = 0.7;

%% 
PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
%         info_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study';
%         save_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\CompiledBlocks_V2';
        %info_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\CompiledBlocks_BehaviorStim';
        %save_path = 'D:\Data\2p\VIPvsNDNF_response_stimuli_study\CompiledBlocks_BehaviorStim';
        info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
        save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\Compiled Blocks v2';
        info_filename = 'Info v2';
        ops_filename = 'Maryse_ops2.mat';
         
    case 'TAKESIANLAB2P' %2P computer
        info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots';
        save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p data\Behavior Pilots\Compiled Blocks';
        info_filename = 'Info';    
        
    case 'RD0332' %Carolyn
       info_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Info Sheets';
        save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Compiled Blocks';
        info_filename = 'Info_rerunapan';
        
    case 'RD-6-TAK2' %Esther's computer
        info_path = 'Z:\Carolyn\2P Imaging data\SSRI study with Jacob';
        save_path = 'D:\test';
        info_filename = 'Info_VxDD053120M2';
        
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
 
%Remove header from Info
Info(1,:) = [];

%Remove rows that are set to "Ignore"
ignore = [Info{:,1}]';
currentInfo = Info(ignore == 0,:);

if isempty(currentInfo)
    error('No data found to compile. Check Info sheet.')
end

%Loop through all remaining rows
for i = 1:size(currentInfo,1)

    %Create setup variable that will contain all the necessary information about the block
    setup = struct;
    setup.constant          =   constant;
    setup.Info              =   Info;                   %Record Info for records only
    setup.pathname          =   [currentInfo{i,2}];     %first part of the path
    setup.mousename         =   [currentInfo{i,3}];     %part of the path, no underscores
    setup.expt_date         =   [currentInfo{i,4}];     %part of the path, YYYY-MM-DD
    setup.block_name        =   [currentInfo{i,5}];     %part of the path, full block name used for BOT
    setup.FOV               =   [currentInfo{i,6}];     %which data to consider as coming from the same field of view, per mouse
    setup.Tosca_session     =   [currentInfo{i,7}];     %Tosca session
    setup.Tosca_run         =   [currentInfo{i,8}];     %Tosca run
    setup.analysis_name     =   [currentInfo{i,9}];     %part of the path, folder where fall.mats are stored
    setup.run_redcell       =   [currentInfo{i,10}];    %do you have red cells? 0 or 1
    setup.VR_name           =   [currentInfo{i,11}];    %full voltage recording name (if widefield only)
    setup.stim_name         =   [currentInfo{i,12}];    %type of stim presentation in plain text
    setup.stim_protocol     =   [currentInfo{i,13}];    %number corresponding to stim protocol
    setup.gcamp_type        =   [currentInfo{i,14}];    %f, m, or s depending on GCaMP type
    setup.expt_group        =   [currentInfo{i,15}];    %name of experimental group or condition
    setup.imaging_depth     =   [currentInfo{i,16}];    %imaging depth in microns
    setup.after_stim        =   [currentInfo{i,17}];    %overwrite constant.after_stim
 
    if ~ismissing(setup.after_stim)
        setup.constant.after_stim = setup.after_stim;
    end
    
    if ~ismissing(setup.FOV)
        try
            FOVtag = ['_FOV' setup.FOV{1}];
        catch
            FOVtag = ['_FOV' num2str(setup.FOV(1))];
        end
    else
        FOVtag = '';
    end

    if ~ismissing(setup.block_name)
        block_number = setup.block_name{1}(1,end-2:end);
        setup.imaging_set = str2double(block_number);
        blockTag = ['_Block' block_number];
    else
        blockTag = '';
    end
    
    if ~ismissing(setup.VR_name)
        widefieldTag = 'widefield-';
    else
        widefieldTag = '';
    end
    
    setup.block_filename = strcat('Compiled_', setup.mousename, FOVtag, '_', setup.expt_date, ...
        '_Session', sprintf('%02d',setup.Tosca_session), '_Run', sprintf('%02d',setup.Tosca_run),...
        blockTag, '_', widefieldTag, setup.stim_name);
    setup.block_supname = strcat(setup.mousename, '-', FOVtag(2:end), '-', setup.expt_date, ...
        '-', setup.stim_name, '-Session', sprintf('%02d',setup.Tosca_session), '-Run', sprintf('%02d',setup.Tosca_run),...
        blockTag(2:end));
    
    %Skip files that have previously been compiled
    if ~recompile
        cd(save_path)
        if isfile(strcat(setup.block_filename, '.mat'))
            disp(setup.block_filename);
            disp('Skipping (already compiled)');
            continue
        end
    end
    
    %Establish and test paths, allowing for paths to be missing
    %TOSCA PATH
    if ismissing(setup.Tosca_session)
        setup.Tosca_path = nan;
    else
        setup.Tosca_path = strcat(setup.pathname, '/', setup.mousename, '/Tosca_', setup.mousename, {'/Session '}, num2str(setup.Tosca_session));
        if ~isfolder(setup.Tosca_path)
            disp(setup.block_filename)
            error('Your Tosca path is incorrect.')
        end
    end
    
    %BLOCK PATH
    if ismissing(setup.block_name)
        setup.block_path = nan;
    else
        setup.block_path   = strcat(setup.pathname, '/', setup.mousename, '/', setup.expt_date, '/', setup.block_name);
        if ~isfolder(setup.block_path)
            disp(setup.block_filename)
            error('Your block path is incorrect.')
        end
    end
    
    %VR PATH
    if ismissing(setup.VR_name)
        setup.VR_path = nan;
    else
        setup.VR_path  = strcat(setup.pathname, '/', setup.mousename, '/', setup.expt_date, '/', setup.VR_name);
        if ~isfolder(setup.VR_path)
            disp(setup.block_filename)
            error('Your voltage recording path is incorrect.')
        end
    end
    
    %SUITE2P PATH
    if ismissing(setup.analysis_name)
        setup.suite2p_path = nan;
    else
        setup.suite2p_path = strcat(setup.pathname, '/', setup.analysis_name);
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
    %determine which trials are considered "active (locomotor)"
    [block] = define_sound_singleblock(block);
    [block] = define_loco_singleblock(block);

    %pull out block-specific data from Fall.mat
    [block] = define_suite2p_singleblock(block, user_ops);
    
    %find the stim-aligned traces
    [block] = align_to_stim(block);
    
    %% Save block
    disp('Saving...');
    cd(save_path)
    save(setup.block_filename, 'block');
end       

disp('Finished compiling all blocks.');