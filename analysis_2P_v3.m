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


    stim_protocol = 2;
    run_redcell = 1;

    stim_protocol = 11;
    run_redcell = 0;
    std_level = 1.5;
    std_level_byStim = 1.5;

    %% Load Info.mat
    % Make setup and data structure out of all blocks that correspond to stim_protocol
    % Later we can also add other things like groups

    PC_name = getenv('computername');

    switch PC_name
        case 'RD0366' %Maryse
            info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            compiled_blocks_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocks';
            save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            info_filename = 'Info_NxDB092719M2';
        case 'RD0332' %Carolyn
            info_path = 'D:\2P analysis\2P local data\Carolyn';
            compiled_blocks_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Daily Imaging';
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
[data] = isresponsive_all(data,std_level,run_redcell);

%% pull out responsive cells by stim type, plot
%magenta = all cells
%blue = positively responsive cells
%cyan=negatively responsive cells - update this 

[data] = isresponsive_byStim(data,std_level_byStim);
[data] = plotbystim(data,run_redcell); 


%% plot +/- locomotion

setup = data.setup;

for a=1:size(setup.mousename,1)
    for b=1:size(setup.mousename,2)
        
        if isempty(setup.mousename{a,b})
            continue;
        end
        
        mouseID=setup.mousename{a,b};
        FOV=setup.FOVs{a,b};
        
         Loc_trial=find(data.([mouseID]).loco);
         noLoc_trial=find(data.([mouseID]).loco==0);
         a_green = squeeze(mean(mean(data.([mouseID]).traces_G(:,Loc_trial,:), 2),1)); %Active trials
         for i=1:size(data.([mouseID]).traces_G,1)
             a_sem = std(data.([mouseID]).traces_G(i,Loc_trial,:))./sqrt(size(data.([mouseID]).traces_G(i,Loc_trial,:),2));
         end
         
         
         b_green = squeeze(mean(mean(data.([mouseID]).traces_G(:,noLoc_trial,:), 2),1)); %Non-active trials
        for i=1:size(data.([mouseID]).traces_G,1)
             b_sem = std(data.([mouseID]).traces_G(i,noLoc_trial,:))./sqrt(size(data.([mouseID]).traces_G(i,noLoc_trial,:),2));
         end
         x=1:length(a_green);
         
         figure
         
         subplot(2,1,1); hold on
         title([mouseID ' FOV ' num2str(FOV)])
         shadedErrorBar(x,smooth(a_green,5),smooth(a_sem,5),'lineprops','-b','transparent',1);
         legend({'Active trials'})
         xlabel('Frames')
         ylabel('Delta F/F')
         
         subplot(2,1,2); hold on
         shadedErrorBar(x,smooth(b_green,5),smooth(b_sem,5),'lineprops','-k','transparent',1);
         legend({'Inactive trials'})
         xlabel('Frames')
         ylabel('Delta F/F')
end
end

%% Save data

if loadPreviousData
    cd(PathName) %Save in the same place you loaded data from
    save([FileName(1:end-4) '_reload'])
else
    cd(save_path)
    d = datestr(now,'yyyymmdd-HHMMSS');
    save(['Data_' d '.mat'], 'data');
end