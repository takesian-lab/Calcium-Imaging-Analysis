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

setup = struct;
setup.username = ''; %'Carolyn'
setup.path_name = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study'; %Include end slash if you're Carolyn
cd(setup.path_name)
load('Info.mat')
setup = fillSetupFromInfoTable(setup, Info);
setup.Info = Info;

%%