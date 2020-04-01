clear all;
%% Noiseburst all cells from suite2p
%Anne Takesian - 2/22/2019
%updated Carolyn, compatible with Python version of Suite2p. Also does Red vs Green cell 7/23/19
%Updated Feb 2020, CGS - put most of the analysis into functions.
%Updated April 2020, MET - V3 created to load compiled blocks

%% define what type of analysis you are doing
%stim protocol code is:
%noiseburst=1
%ReceptiveField=2
%FM sweep=3
%SAM = 6
%widefield=4
%SAM freq = 6

stim_protocol=1;

%% Load Info.mat
% Make setup and data structure out of all blocks that correspond to stim_protocol
% Later we can also add other things like groups

info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
compiled_blocks_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocks';
%save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study'; %To save figures later
cd(info_path)
Info = importfile('Info');

%Create setup variable for files corresponding to stim_protocol
setup = struct;
setup.Info = Info;
setup.stim_protocol = stim_protocol;
[data, setup] = fillSetupFromInfoTable_v2(setup, Info, compiled_blocks_path);

%% Now find processed suite2P data
if setup.run_redcell==0
    [data]=Noiseburst_analysis_greenonly_v2(data,setup);
    %red cells need to be updated and checked to make sure that they work.
elseif setup.run_redcell==1
    [data,traces_R,traces_G]=Noiseburst_analysis(a,Frames,Frame_rate,Imaging_Block_String,Imaging_Num,mouseID,date,Sound_Time,...
        timestamp,i,analysis_folder,path_name,length_sound_trial_first,username,data);
end
%% sound responsive cells - all sounds averaged together, plot means
% in figure1, you get a grid of mean activity per cell.
% Black=nonresponsive, blue=responsive, and cyan=negatively responsive
% in figure2, you get 3 images - magenta = mean of all cells, blue = mean
% of responsive cells, and cyan = mean of negatively responsive cells

%this only works for green data currently
std_level = 1.5;%set this here to change std
[data] = isresponsive_all(data,setup,std_level)
%
clear std_level

%% pull out responsive cells by stim type, plot
%magenta = all cells
%blue = positively responsive cells
%cyan=negatively responsive cells - update this 

std_level = 1.5;%set this here to change std 
[data] = isresponsive_byStim(data,setup,std_level)
clear std_level
[data]= plotbystim(setup,data)



%% plot +/- locomotion
for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)};
         Loc_trial=find(data.([mouseID]).parameters.loco);
         noLoc_trial=find(data.([mouseID]).parameters.loco==0);
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
         title(mouseID)
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



%% save data
% folder = 'Z:\Carolyn\Figures for SfN 2019';
% addpath(folder);
% save('noiseburst_analysis.mat', 'data');

