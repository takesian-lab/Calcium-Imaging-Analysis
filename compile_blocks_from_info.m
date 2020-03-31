%% compile_blocks_from_info
%  Save one compiled 'block.mat' file for each block listed in info.mat
%  A compiled block contains the 2p and Tosca data for that recording session
%
%  If data is missing, the script will still run in order to allow for the
%  user to look at Tosca data separate from 2p data or vice versa.
%
%  Uses the functions compile_block, fillSetupFromInfoTable
%
%  Use the script verify_block to preview the block contents prior to analysis
%
%  Maryse Thomas - March 2020

%% Load Info

info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
cd(info_path)
load('Info.mat')

%Info.mat Table is a variable that stores all recording file information
%where each row of Info corresponds to a single recording block

%  Currently, the column order of Info is:
I = 1; %ignore (0 or 1 allowing user to "hide" some files from analysis)
U = 2; %user name %Not used right now
M = 3; %mouse name
P = 4; %path name (path to main folder of mouse's data) %Not used right now
D = 5; %experiment date
R = 6; %ROI (which data to consider as coming from the same ROI)
IS = 7; %imaging set (BOT #s)
TS = 8; %Tosca session #
TR = 9; %Tosca run #
AP = 10; %analysis path name (path to suite2p analysis)
FR = 11; %frame rate (15 or 30)
RR = 12; %do you have red cells? 0 or 1 %Not used right now
VR = 13; %voltage recording (0 for widefield, 1 for 2p)
SN = 14; %stim name (type of stim presentation in plain text)
SP = 15; %stim protocol (number corresponding to stim)

%% Compile all blocks unless they are set to "Ignore"

%Remove rows that are set to "Ignore"
ignore = [Info{:,I}]';
currentInfo = Info(ignore == 0,:);

%Loop through all remaining rows
for i = 1:size(currentInfo,1)

    %Create setup variable that will contain all the necessary information about the block
    setup = struct;
    setup.Info              =   Info;
    setup.username          =   [currentInfo_R{i,U}];
    setup.mousename         =   [currentInfo_R{i,M}];
    setup.pathname          =   [currentInfo_R{i,P}];
    setup.expt_date         =   [currentInfo_R{i,D}];
    setup.ROI               =   [currentInfo_R{i,R}];
    setup.imaging_set       =   [currentInfo_R{i,IS}];
    setup.Tosca_session     =   [currentInfo_R{i,TS}];
    setup.Tosca_run         =   [currentInfo_R{i,TR}];
    setup.analysis_name     =   [currentInfo_R{i,AP}];
    setup.framerate         =   [currentInfo_R{i,FR}];
    setup.run_redcell       =   [currentInfo_R{i,RR}];
    setup.voltage_recording =   [currentInfo_R{i,VR}];
    setup.stim_name         =   [currentInfo_R{i,SN}];
    setup.stim_protocol     =   [currentInfo_R{i,SP}];

    setup.block_name = strcat('Compiled_', setup.mousename, '_', setup.expt_date, ...
        '_Block', setup.imaging_set, '_Session', setup.Tosca_session, '_Run', setup.Tosca_run, ...
        setup.stim_name);
    
    suite2p_path = strcat(setup.pathname, '/', setup.mousename, '/', setup.analysis_name);
    Tosca_path = strcat(setup.pathname, '/', setup.mousename, '/Tosca_', setup.mousename);
    
    block = compile_block(setup, suite2p_path, Tosca_path);

    cd(save_path)
    save(setup.block_name, 'block.mat');
    
end        