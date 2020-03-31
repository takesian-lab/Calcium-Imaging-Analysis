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

%% Create setup variable for files in Info.mat

%% Info.mat Table is a variable that stores all recording file information

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

%% Look for files that match stim_protocol
%Later we could update this to also look for other parameters

%stim protocol code is:
%noiseburst=1
%ReceptiveField=2
%FM sweep=3
%SAM = 5
%widefield=4
%SAM freq = 6
code = {'Noiseburst', 'Receptive Field', 'FM sweep', 'SAM', 'Widefield', 'SAM freq'};
display(['Analyzing ' code{setup.stim_protocol} ' files'])

%Find table rows that match stim_protocol
stims = [Info{:,SP}]';
matching_stims = stims == setup.stim_protocol;
currentInfo = Info(matching_stims,:);

%Remove rows that are set to "Ignore"
ignore = [currentInfo{:,I}]';
currentInfo = currentInfo(ignore == 0,:);

setup = struct;
setup.username = ''; %'Carolyn'
setup.path_name = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study'; %Include end slash if you're Carolyn
cd(setup.path_name)
load('Info.mat')
setup = fillSetupFromInfoTable(setup, Info);
setup.Info = Info;


        tempFrame_set = {};
        for k = 1:size(currentInfo_R,1)
            path = strcat(currentInfo_R{k,P}, '/', currentInfo_R{k,M}, '/', currentInfo_R{k,AP});
            cd(path);
            load('Fall.mat', 'ops');
            Imaging_Block = currentInfo_R{k,IS};
            showTable = 0;
            if k == 1
                display(currentMouse);
                showTable = 1;
            end
            tempFrame_set{1,k} = get_frames_from_Fall(ops,Imaging_Block,showTable);
        end
        
        %Fill setup with structure {Mouse, ROI}
        setup.mousename{i,j}         =   currentMouse;
        setup.expt_date{i,j}         =   [currentInfo_R{:,D}];
        setup.Imaging_sets{i,j}      =   [currentInfo_R{:,IS}];
        setup.Session{i,j}           =   [currentInfo_R{:,TS}];
        setup.Tosca_Runs{i,j}        =   [currentInfo_R{:,TR}];
        setup.analysis_name{i,j}     =   [currentInfo_R{:,AP}];
        setup.framerate{i,j}         =   [currentInfo_R{:,FR}];
        %setup.run_redcell{i,j}       =   [currentInfo_R{:,RR}];
        setup.voltage_recording{i,j} =   [currentInfo_R{:,VR}];
        setup.Frame_set{i,j}         =   tempFrame_set;
        
%%