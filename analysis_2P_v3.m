clear all;

%% Noiseburst all cells from suite2p
%Anne Takesian - 2/22/2019
%updated Carolyn, compatible with Python version of Suite2p. Also does Red vs Green cell 7/23/19
%Updated Feb 2020, CGS - put most of the analysis into functions.
%Updated April 2020, MET - V3 created to load compiled blocks

%% Load Data if it already exists, otherwise create new Data struct

loadPreviousData = 0;

if loadPreviousData
    %Load data
    [FileName,PathName] = uigetfile('.mat');
    load([PathName '/' FileName])
else

    %% Magic numbers and define what type of analysis you are doing
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


    stim_protocol = 1;
    run_redcell = 1;
    std_level = 3;
    std_level_byStim = 3;

    %% Load Info.mat
    % Make setup and data structure out of all blocks that correspond to stim_protocol
    % Later we can also add other things like groups

    PC_name = getenv('computername');

    switch PC_name
        case 'RD0366' %Maryse
            info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            %compiled_blocks_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocks';
            compiled_blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\analyzed\Daily_Imaging';
            save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            info_filename = 'Info_VxAC031419M1';
        case 'RD0332' %Carolyn
            info_path = 'D:\2P analysis\2P local data\Carolyn';
%             compiled_blocks_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Daily Imaging';
            compiled_blocks_path = 'Z:\Maryse\2p analysis\Compiled blocks';
            save_path = 'D:\2P analysis\2P local data\Carolyn';
            info_filename = 'Info';
        case 'RD0386' %Wisam
            % INSERT PATHS HERE
            info_filename = 'Info';
        otherwise
            disp('Computer does not match known users')
    end

    cd(info_path)
    Info = importfile(info_filename);

    %Create data structure for files corresponding to stim_protocol
    [data] = fillSetupFromInfoTable_v2(Info, compiled_blocks_path, stim_protocol);
    data.setup.run_redcell = run_redcell;
end

%% Now find processed suite2P data
    [data]=concatBlocks_aligned(data);
    [data]=df_F(data);

%% sound responsive cells - all sounds averaged together, plot means
% in figure1, you get a grid of mean activity per cell.
% Black=nonresponsive, blue=responsive, and cyan=negatively responsive
% in figure2, you get 3 images - magenta = mean of all cells, blue = mean
% of responsive cells, and cyan = mean of negatively responsive cells

%this only works for green data currently
[data] = isresponsive_all(data,std_level);

%% pull out responsive cells by stim type, plot


[data] = isresponsive_byStim(data,std_level_byStim);
[data] = plotbystim(data,run_redcell); 


%% plot +/- locomotion
[data] = plot_loco_by_stim(data);

%% Save data

if loadPreviousData
    cd(PathName) %Save in the same place you loaded data from
    save([FileName(1:end-4) '_reload'])
else
    cd(save_path)
    d = datestr(now,'yyyymmdd-HHMMSS');
    save(['Data_' d '.mat'], 'data');
end