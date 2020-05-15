%% CLEAN AND CLEAR

clear 
close all 
clc

%% Noiseburst all cells from suite2p
%
% Anne Takesian - 2/22/2019
% Updated Carolyn, compatible with Python version of Suite2p. Also does Red vs Green cell 7/23/19
% Updated Feb 2020, CGS - put most of the analysis into functions.
% Updated April 2020, MET - V3 created to load compiled blocks
% Updated April 2020, by Wisam Reid

%% Define what type of analysis you are doing

% stim protocol code is:
% noiseburst=1
% ReceptiveField=2
% FM sweep=3
% SAM = 6
% widefield=4
% SAM freq = 6

% TODO: Why is this not directly extracted from the table?
% I think this meant to selective on a larger table with potentially
% multiple stim protocols to choose from.
stim_protocol = 1;

%% Set Data Paths and Load Info.mat
% Make setup and data structure out of all blocks that correspond to stim_protocol

info_path = '/Users/wisamreid/Documents/School/Research (Harvard)/Takesian/2P/Thy1 Experiments/YD111219F3-2P-noisebursts';
compiled_blocks_path = '/Users/wisamreid/Documents/School/Research (Harvard)/Takesian/2P/Thy1 Experiments/YD111219F3-2P-noisebursts/Compiled';
info_filename = 'Info';

% TODO: Do we need to cd into a different directory?
cd(info_path)
Info = importfile(info_filename);

%%

% Create setup variable for files corresponding to stim_protocol
setup = struct;
setup.Info = Info;
setup.stim_protocol = stim_protocol;
setup.run_redcell = 0;
% TODO: This is broken
% [data, setup] = fillSetupFromInfoTable_v2(setup, Info, compiled_blocks_path);


%%%%%% FINISH NOTES IN THIS FILE
[data, setup] = fillSetupFromInfoTable_wisam(setup, Info, compiled_blocks_path);
data.setup = setup; %Save this info

%% Now find processed suite2P data

% TODO: why do we need to do this with the red cells?
if setup.run_redcell == 0
    
    %%%%%% FINISH NOTES IN THIS FILE
    % [data] = Noiseburst_analysis_greenonly_v2(data,setup);
    [data] = Noiseburst_analysis_greenonly_wisam(data,setup);
    % Red cells need to be updated and checked to make sure that they work.
elseif setup.run_redcell==1
% % % % % [data,traces_R,traces_G] = Noiseburst_analysis(a,Frames,Frame_rate,Imaging_Block_String,Imaging_Num,mouseID,date,Sound_Time,...
% % % % %     timestamp,i,analysis_folder,path_name,length_sound_trial_first,username,data);
    disp('Tried running red cells...')
end

%% sound responsive cells - all sounds averaged together, plot means
% Figure 1, you get a grid of mean activity per cell.
% Black = nonresponsive, blue = responsive, and cyan = negatively responsive
% in figure2, you get 3 images - magenta = mean of all cells, blue = mean
% of responsive cells, and cyan = mean of negatively responsive cells

% TODO: 
% This only works for green data currently
% TODO: Is there a better way to do this?
% Perhaps we could calculate the noise floor and measure from there.
% i.e. a threshold on the probability of noise
std_level = 0; % 1.5; % set this here to change std
% TODO: do we need both data and setup?
% [data] = isresponsive_all(data,setup,std_level)
[data] = isresponsive_all_wisam(data,setup,std_level)
clear std_level

%% pull out responsive cells by stim type, plot

% magenta = all cells
% blue = positively responsive cells
% cyan = negatively responsive cells - update this 

% % % % % % % std_level = 0; %1.5; % set this here to change std 
% % % % % % % [data] = isresponsive_byStim(data,setup,std_level)
% % % % % % % clear std_level
% % % % % % % [data]= plotbystim(setup,data)

%% plot +/- locomotion

% % % % % % % for a=1:length(setup.mousename)
% % % % % % %     
% % % % % % %     mouseID=setup.mousename{(a)};
% % % % % % %     Loc_trial=find(data.([mouseID]).parameters.loco);
% % % % % % %     noLoc_trial=find(data.([mouseID]).parameters.loco==0);
% % % % % % %     a_green = squeeze(mean(mean(data.([mouseID]).traces_G(:,Loc_trial,:), 2),1)); % Active trials
% % % % % % %     
% % % % % % %     for i=1:size(data.([mouseID]).traces_G,1)
% % % % % % %         a_sem = std(data.([mouseID]).traces_G(i,Loc_trial,:))./sqrt(size(data.([mouseID]).traces_G(i,Loc_trial,:),2));
% % % % % % %     end
% % % % % % % 
% % % % % % %     b_green = squeeze(mean(mean(data.([mouseID]).traces_G(:,noLoc_trial,:), 2),1)); %Non-active trials
% % % % % % %     
% % % % % % %     for i=1:size(data.([mouseID]).traces_G,1)
% % % % % % %         b_sem = std(data.([mouseID]).traces_G(i,noLoc_trial,:))./sqrt(size(data.([mouseID]).traces_G(i,noLoc_trial,:),2));
% % % % % % %     end
% % % % % % %     
% % % % % % %     x=1:length(a_green);
% % % % % % % 
% % % % % % %     figure
% % % % % % % 
% % % % % % %     subplot(2,1,1); hold on
% % % % % % %     title(mouseID)
% % % % % % %     shadedErrorBar(x,smooth(a_green,5),smooth(a_sem,5),'lineprops','-b','transparent',1);
% % % % % % %     legend({'Active trials'})
% % % % % % %     xlabel('Frames')
% % % % % % %     ylabel('Delta F/F')
% % % % % % % 
% % % % % % %     subplot(2,1,2); hold on
% % % % % % %     shadedErrorBar(x,smooth(b_green,5),smooth(b_sem,5),'lineprops','-k','transparent',1);
% % % % % % %     legend({'Inactive trials'})
% % % % % % %     xlabel('Frames')
% % % % % % %     ylabel('Delta F/F')
% % % % % % % end

%% save data

% folder = 'Z:\Carolyn\Figures for SfN 2019';
% addpath(folder);
% save('noiseburst_analysis.mat', 'data');

