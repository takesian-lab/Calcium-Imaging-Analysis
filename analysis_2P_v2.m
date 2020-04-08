clear all;
%% Noiseburst all cells from suite2p
%Anne Takesian - 2/22/2019
%updated Carolyn, compatible with Python version of Suite2p. Also does Red vs Green cell 7/23/19
%Updated Feb 2020, CGS - put most of the analysis into functions.

% Maryse is awesome!

%define what type of analysis you are doing
%stim protocol code is:
%noiseburst=1
%ReceptiveField=2
%FM sweep=3
%SAM = 5
%widefield=4
%SAM freq = 6
%Go/No-go behavior = 7

stim_protocol=1;

%Create setup variable for files corresponding to stim_protocol using Info.mat
setup = struct;
setup.username = ''; %'Carolyn'
%setup.path_name = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
setup.path_name = 'D:/Data/ReceptiveField example';
setup.stim_protocol = stim_protocol;
setup.run_redcell = 0;
cd(setup.path_name)
load('Info.mat')
setup = fillSetupFromInfoTable(setup, Info);
setup.Info = Info;

%% Behavior, locomotion, and sound

%pull out the Tosca-derived, behaviorally relevant data
[data] = behavior_RF(setup);

%pull out the Bruker-derived timestamps
[data] = define_sound(data,setup);

%determine which trials are considered "active (locomotor)"
% [loco_activity,isLocoSound,data] = isLoco(setup,data);
[data] = define_loco(setup,data);

clear adjusted_times timestamp Var1 Var2 New_Sound_Times isLocoSound...
    loco_activity loco_times locTime2 New_sound_times Sound_Time sound_v1 sound_v2 start_time;

%% Now find processed suite2P data
if setup.run_redcell==0
    [data]=Noiseburst_analysis_greenonly(data,setup);
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
[data]=isresponsive_all(data,setup,std_level)
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

