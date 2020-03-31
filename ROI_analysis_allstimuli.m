clear all;
%% Noiseburst all cells from suite2p 
    %Anne Takesian - 2/22/2019
    %updated Carolyn, compatible with Python version of Suite2p. Also does Red vs Green cell 7/23/19
    
%Define analysis-specific info here:        
        setup.username='Carolyn';
        setup.mousename={'VxDB100819F2'};%,'VxDB070919M3','VxDB082019M1'};
        setup.path_name = '/Volumes/Seagate Backup Plus Drive/VIP_NdnF_project/2P local data/';
        %when was the experiment performed?
        setup.expt_date={'2020-01-30'};%;'2019-09-24';'2019-12-04'};
        setup.Imaging_sets = [1]; %BOT numbers
setup.Session = {'2'};%tosca
setup.Tosca_Runs = [1];%;2;1];%tosca
setup.analysis_name = {'ROI_1'};%;'noiseburst2';'noiseburst1'}; %what did you call the analyzed data folder?
setup.Frame_set = {(1:3676)};%(1:2411);(1:1969)};
setup.framerate= {'30'};%,'15','15'};
setup.run_redcell=0;%do you have red cells? y/n
 setup.voltage_recording = [1];% this will always be 1 for 2p, this number changes for widefield
% cd('D:\');

%ddfine what type of analysis you are doing
setup.stim_protocol=1;
%stim protocol code is:
%noiseburst=1
%ReceptiveField=2
%FM sweep=3
%widefield=4

%% Behavior, locomotion, and sound
%this works right now as long as it is one animal and one BOT file - it may
%break with more.

for a=1:length(mousename)
    mouseID=mousename{(a)}
    Tosca_Session=Session{(a)}
    date=expt_date{(a)};
    Imaging_Block=Imaging_sets(a,:)
    Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder
    Tosca_Run_number = num2str(Tosca_Runs(a));
    folder = sprintf([path_name username '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]);
    cd(folder)
    
    Var1=[]; Var2=[]; isLocoSound=[];sound_v1=[]; sound_v2=[];
    
    for i=1:length(Imaging_Block(a,:))
        Imaging_Block_String = num2str(Imaging_Block(i))
        Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
        %pull out the Tosca-derived, behaviorally relevant data
        [Var1,Var2,New_sound_times,start_time,loco_times,loco_activity] = behavior_RF(stim_protocol,a,username,mouseID,Tosca_Session,Tosca_Runs,path_name,Tosca_Run_number,Var1,Var2);
 
        %pull out the Bruker-derived timestamps
        [Sound_Time,timestamp,sound_v1,sound_v2,locTime2]=define_sound(sound_v1,sound_v2,stim_protocol,i,path_name,username,mouseID,date,Var1,Var2,New_sound_times,voltage_recording,start_time,Imaging_Block_String,Imaging_Num,loco_times,isLocoSound);
        %determine which trials are considered "active (locomotor)"
        [loco_activity,isLocoSound] = isLoco(isLocoSound,loco_activity,Sound_Time,locTime2);
        
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp=timestamp;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time=Sound_Time;
    end
        data.([mouseID]).parameters(:,1)=isLocoSound;%index of locomotor times
        data.([mouseID]).parameters(:,2)=sound_v1;%index of variable1 (frequency)
        data.([mouseID]).parameters(:,3)=sound_v2;%index of variable 2 (level)
        clear adjusted_times timestamp Var1 Var2 New_Sound_Times isLocoSound...
            loco_activity loco_times locTime2 New_sound_times Sound_Time;
        
    end

%% Now find processed suite2P data 


length_sound_trial_first =[];
for a=1:length(mousename)
    mouseID=mousename{(a)}
    Frame_rate=framerate{(a)};
    date=expt_date{(a)};
    Frames=cell2mat(Frame_set(a));
    analysis_folder=analysis_name{(a)}
    Imaging_Block=Imaging_sets(a,:)
    traces_G=[];
    if run_redcell==1
        traces_R=[];
    end

    for i=1:length(Imaging_Block(a))
        Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
        Imaging_Block_String = num2str(Imaging_Block(i));
        timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
        Sound_Time = data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;
        folder = sprintf([path_name username '/analyzed/' mouseID '/' date '/' analysis_folder '/suite2p/plane0'])
        cd(folder);
       
        if run_redcell==1
            [data,traces_R,traces_G]=Noiseburst_analysis(a,Frames,Frame_rate,Imaging_Block_String,Imaging_Num,mouseID,date,Sound_Time,...
                timestamp,i,analysis_folder,path_name,length_sound_trial_first,username,data);
   end
    
    if run_redcell==0

    [data,traces_G]=Noiseburst_analysis_greenonly(traces_G,a,Frames,Frame_rate,Imaging_Block_String,Imaging_Num,mouseID,date,Sound_Time,...
        timestamp,i,analysis_folder,path,length_sound_trial_first,username,data);
       data.([mouseID]).traces.green.raw=traces_G;
%        peak_G=
%         peak_maxG=
%         peak_minG=
%         baseline=
%         raw_trace_G=
       

    end
    end
end
%% sound responsive cells - all sounds
all_cells_trace_green = data.combined.trace.green_all;
all_cells_response_green = data.combined.response.green_all;
all_cell_peak_green = data.combined.peak.green_all;
all_cells_avg_around_peak_green = data.combined.avg_around_peak_green;
all_cell_STD_green = data.combined.STDbase.green_all;

if run_redcell==1
all_cells_trace_red = data.combined.trace.red_all;
all_cells_response_red = data.combined.response.red_all;
all_cell_peak_red = data.combined.peak.red_all;
all_cells_avg_around_peak_red = data.combined.avg_around_peak_red;
all_cell_STD_red = data.combined.STDbase.red_all;
end

std_level = 1.25;%set this here to change std

%find green responsive cells and plot results
for i=1:size(all_cells_response_green,1);
    a_green(i,:) = mean(all_cells_trace_green(i,:,:), 2); %mean of all trials for each cell
    SEM_a_green(i,:) = std(all_cells_trace_green(i,:,:))./sqrt(size(all_cells_trace_green(i,:,:),2));
    
    b_green(i,:) = mean(all_cells_response_green(i,:), 2); %average means responses (sound to 2s post sound) across trials for each cell
    c_green(i,:) = mean(all_cell_peak_green(i,:), 2); %average peak response
    d_green(i,:) = mean(all_cells_avg_around_peak_green(i,:), 2);
    e_green(i,:) = mean(all_cells_neg_avg_around_peak_green(i,:), 2); %average around negative peak
    
    %determine whether each cell is responsive (defined as average response
    %more than 2 STDS above baseline) 
    f = mean(all_cell_STD_green(i,:), 2); %average baseline STD across trials for each cell
    isResponsive_green(i) = d_green(i,:) > std_level*mean(f) & b_green(i,:)>0; %will be 0 or 1
    isNegResponsive_green(i) = e_green(i,:) < -std_level*mean(f) & b_green(i,:)<0;     
end

first_cell = 1;
last_cell = size(all_cells_response_green,1);
num_cells = last_cell-first_cell+1;
figure;

for i=first_cell:last_cell  
    %plot mean traces across cells with response means - responsive cells are
    %green, non-responsive are black 
    subplot_num = i-first_cell+1;
    subplot(ceil(sqrt(num_cells)),ceil(sqrt(num_cells)),subplot_num);
    x_green = 1:length(a_green(i,:));
    if isResponsive_green(i) == 1
        shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-b','transparent',1); hold on;
    else
        if isNegResponsive_green(i)== 1
            shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-c','transparent',1); hold on;
        else
            shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-k','transparent',1); hold on;
        end
    end
    
    plot(b_green(i),'o'); %plot average response
    plot(d_green(i),'o'); %plot average around peak
    plot(e_green(i),'o'); 
end  

clear first_cell last_cell num_cells d

%find red responsive cells and plot results
if run_redcell==1
for i=1:size(all_cells_response_red,1);
    a_red(i,:) = mean(all_cells_trace_red(i,:,:), 2); %mean of all trials for each cell
    SEM_a_red(i,:) = std(all_cells_trace_red(i,:,:))./sqrt(size(all_cells_trace_red(i,:,:),2));
    
    b_red(i,:) = mean(all_cells_response_red(i,:), 2); %average means responses (sound to 2s post sound) across trials for each cell
    c_red(i,:) = mean(all_cell_peak_red(i,:), 2); %average peak response
    d_red(i,:) = mean(all_cells_avg_around_peak_red(i,:), 2);
    e_red(i,:) = mean(all_cells_neg_avg_around_peak_red(i,:), 2); %average around negative peak
    
    %determine whether each cell is responsive (defined as average response
    %more than 2 STDS above baseline) 
    f = mean(all_cell_STD_red(i,:), 2); %average baseline STD across trials for each cell
    isResponsive_red(i) = d_red(i,:) > std_level*mean(f) & b_red(i,:)>0; %will be 0 or 1
    isNegResponsive_red(i) = e_red(i,:) < -std_level*mean(f)  & b_red(i,:)<0;     
end


first_cell = 1;
last_cell = size(all_cells_response_red,1);
num_cells = last_cell-first_cell+1;
figure;

for i=first_cell:last_cell  
    %plot mean traces across cells with response means - responsive cells are
    %red, non-responsive are black 
    subplot_num = i-first_cell+1;
    subplot(ceil(sqrt(num_cells)),ceil(sqrt(num_cells)),subplot_num);
    x_red = 1:length(a_red(i,:));
    if isResponsive_red(i) == 1
        shadedErrorBar(x_red,smooth((a_red(i,:)),10),smooth((SEM_a_red(i,:)),10),'lineprops','-r','transparent',1); hold on;
    else
        if isNegResponsive_red(i)== 1
            shadedErrorBar(x_red,smooth((a_red(i,:)),10),smooth((SEM_a_red(i,:)),10),'lineprops','-m','transparent',1); hold on;
        else
            shadedErrorBar(x_red,smooth((a_red(i,:)),10),smooth((SEM_a_red(i,:)),10),'lineprops','-k','transparent',1); hold on;
        end
    end
    
    plot(b_red(i),'o'); %plot average response
    plot(d_red(i),'o'); %plot average around peak
    plot(e_red(i),'o'); 
end  
end


data.combined.isResponsiveGreen=isResponsive_green';
data.combined.isNegResponsive_green = isNegResponsive_green';

if run_redcell==1
data.combined.isResponsiveRed=isResponsive_red';
data.combined.isNegResponsive_red = isNegResponsive_red';
end
% Plot traces of all Positive Responsive Cells  w SEM (for both red and green
% cells)
timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
responsive_cells_green = a_green(isResponsive_green,:,:);
mean_responsive_green = mean(responsive_cells_green, 1);
SEM_responsive_green = std(responsive_cells_green)./sqrt(size(responsive_cells_green,1));
x_green = timestamp(1:length(mean_responsive_green));

figure; hold on;
shadedErrorBar(x_green,smooth((mean_responsive_green),10),smooth((SEM_responsive_green),10),'lineprops','-b','transparent',1); hold on;

if run_redcell==1
responsive_cells_red = a_red(isResponsive_red,:,:);
mean_responsive_red = mean(responsive_cells_red, 1);
SEM_responsive_red = std(responsive_cells_red)./sqrt(size(responsive_cells_red,1));
 x_red = timestamp(1:length(mean_responsive_red));
 shadedErrorBar(x_red,smooth((mean_responsive_red),10),smooth((SEM_responsive_red),10),'lineprops','-r','transparent',1);
end

 



title('Mean of All Positively Responsive Cells');

% Plot traces of all Negative Responsive Cells  w SEM (for both red and green
% cells)
timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
neg_responsive_cells_green = a_green(isNegResponsive_green,:,:);
neg_mean_responsive_green = mean(neg_responsive_cells_green, 1);
neg_SEM_responsive_green = std(neg_responsive_cells_green)./sqrt(size(neg_responsive_cells_green,1));
 x_green = timestamp(1:length(neg_mean_responsive_green));
figure; hold on;
shadedErrorBar(x_green,smooth((neg_mean_responsive_green),10),smooth((neg_SEM_responsive_green),10),'lineprops','-b','transparent',1); hold on;

if run_redcell==1
neg_responsive_cells_red = a_red(isNegResponsive_red,:,:);
neg_mean_responsive_red = mean(neg_responsive_cells_red, 1);
neg_SEM_responsive_red = std(neg_responsive_cells_red)./sqrt(size(neg_responsive_cells_red,1));
 x_red = timestamp(1:length(neg_mean_responsive_red));
 shadedErrorBar(x_red,smooth((neg_mean_responsive_red),10),smooth((neg_SEM_responsive_red),10),'lineprops','-r','transparent',1);

end
title('Mean of All Negatively Responsive Cells');

%% sound responsive cells - by stim


%% mean across all cells

%% mean across stim parameters

%%
%% Plot traces of all cells with locomotion during sounds versus without
%green
mean_loco_green = mean(all_cells_isLocoSound_green,1);
SEM_loco_green = std(all_cells_isLocoSound_green,1)./sqrt(size(all_cells_isLocoSound_green,1));

mean_noloco_green = mean(all_cells_noLocoSound_green,1);
SEM_noloco_green = std(all_cells_noLocoSound_green,1)./sqrt(size(all_cells_noLocoSound_green,1));

%red
mean_loco_red = mean(all_cells_isLocoSound_red,1);
SEM_loco_red = std(all_cells_isLocoSound_red,1)./sqrt(size(all_cells_isLocoSound_red,1));

mean_noloco_red = mean(all_cells_noLocoSound_red,1);
SEM_noloco_red = std(all_cells_noLocoSound_red,1)./sqrt(size(all_cells_noLocoSound_red,1));

timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
x = squeeze(timestamp(1:length(mean_loco_green)));

  figure; 
       shadedErrorBar(x,smooth((mean_loco_green),5),smooth((SEM_loco_green),5),'lineprops','-b','transparent',1); hold on;
       shadedErrorBar(x,smooth((mean_noloco_green),5),smooth((SEM_noloco_green),5),'lineprops','-k','transparent',1); hold on;
       title('Effect of Locomotion - non VIP neurons, green is loco')
     %  legend('locomotion','no locomotion')
       
  figure;     
      
       shadedErrorBar(x,smooth((mean_loco_red),5),smooth((SEM_loco_red),5),'lineprops','-r','transparent',1); hold on;
       shadedErrorBar(x,smooth((mean_noloco_red),5),smooth((SEM_noloco_red),5),'lineprops','-k','transparent',1); 
       title('Effect of Locomotion - VIP neurons, red is loco')
       %legend('locomotion','no locomotion'); hold on;


        
%% Find means across all cells
all_cells_trace_green =[];% this matrix will have trace/green cell - will concatenate imaging blocks in this section
all_cells_response_green =[];%response to sound, no baseline - measured from frame before sound
all_cell_peak_green = [];%peak responses/green cell
all_cells_avg_around_peak_green = [];
all_cells_neg_avg_around_peak_green = [];
all_cell_STD_green = [];

if run_redcell==1
all_cells_trace_red = []; %this matrix will have trace/red cell - will concatenate imaging blocks in this section
all_cells_response_red= [];
all_cell_peak_red = [];
all_cells_avg_around_peak_red = [];
all_cells_neg_avg_around_peak_red = [];
all_cell_STD_red = [];
end


for a=1:length(mousename)
    mouseID=mousename{(a)};
    date=expt_date{(a)};
    Imaging_Block=Imaging_sets(a,:)
for i=1:length(Imaging_Block(i))
    Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
    trace_around_sound_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_green;%trace(dF/F) around sound for each non-VIP (cellxsoundxtime matrix)
    
    avg_sound_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_sound_green;
    avg_around_peak_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_around_peak_green;
    neg_avg_around_peak_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).neg_avg_around_peak_green;
    peak_sound_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).peak_sound_green;
    std_baseline_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).std_baseline_green; 
 %concatenate the imaging blocks
    all_cells_trace_green = [all_cells_trace_green; trace_around_sound_green]; 
    all_cells_response_green = [all_cells_response_green; avg_sound_green];
    all_cell_peak_green = [all_cell_peak_green; peak_sound_green];
    all_cells_avg_around_peak_green = [all_cells_avg_around_peak_green; avg_around_peak_green];
    all_cells_neg_avg_around_peak_green = [all_cells_neg_avg_around_peak_green; neg_avg_around_peak_green];
    all_cell_STD_green = [all_cell_STD_green; std_baseline_green];
end
end
data.combined.trace.green_all=all_cells_trace_green;
data.combined.response.green_all=all_cells_response_green;
data.combined.peak.green_all=all_cell_peak_green;
data.combined.avg_around_peak_green = all_cells_avg_around_peak_green;
data.combined.neg_avg_around_peak_green = all_cells_neg_avg_around_peak_green;
data.combined.STDbase.green_all=all_cell_STD_green;

mean_cell_green = squeeze(mean(mean(all_cells_trace_green(:,:,:),1),2));
mean_response_green = mean(all_cells_response_green(:,:),2);
mean_peak_green = mean(all_cell_peak_green(:,:),2);
data.combined.trace.green_mean=mean_cell_green;
data.combined.response.green_mean=mean_response_green;
data.combined.peak.green_mean=mean_peak_green;

%red
if run_redcell==1
for a=1:length(mousename)
    mouseID=mousename{(a)};
    date=expt_date{(a)};
for i=1:length(Imaging_Block(a))
    Imaging_Num =  sprintf( '%03d', Imaging_Block(a));
    trace_around_sound_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_red;%trace(dF/F) around sound for each non-VIP (cellxsoundxtime matrix)
   if size(trace_around_sound_red,3)==52
        trace_around_sound_red=cat(3,trace_around_sound_red,zeros(size(trace_around_sound_red,1),size(trace_around_sound_red,2),1));
        trace_around_sound_red(trace_around_sound_red==0)=NaN;
    end
    avg_sound_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_sound_red;
    avg_around_peak_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_around_peak_red;
    neg_avg_around_peak_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).neg_avg_around_peak_red;
    peak_sound_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).peak_sound_red;
    std_baseline_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).std_baseline_red; 
 %concatenate the imaging blocks
    all_cells_trace_red = [all_cells_trace_red; trace_around_sound_red]; 
    all_cells_response_red = [all_cells_response_red; avg_sound_red];
    all_cell_peak_red = [all_cell_peak_red; peak_sound_red];
    all_cells_avg_around_peak_red = [all_cells_avg_around_peak_red; avg_around_peak_red];
    all_cells_neg_avg_around_peak_red = [all_cells_neg_avg_around_peak_red; neg_avg_around_peak_red];
    all_cell_STD_red = [all_cell_STD_red; std_baseline_red];
end
end
data.combined.trace.red_all=all_cells_trace_red;
data.combined.response.red_all=all_cells_response_red;
data.combined.peak.red_all=all_cell_peak_red;
data.combined.avg_around_peak_red = all_cells_avg_around_peak_red;
data.combined.neg_avg_around_peak_red = all_cells_neg_avg_around_peak_red;
data.combined.STDbase.red_all=all_cell_STD_red;

mean_cell_red = squeeze(mean(mean(all_cells_trace_red(:,:,:),1),2));
mean_response_red = mean(all_cells_response_red(:,:),2);
mean_peak_red = mean(all_cell_peak_red(:,:),2);
data.combined.trace.red_mean=mean_cell_red;    
data.combined.response.red_mean=mean_response_red;
data.combined.peak.red_mean=mean_peak_red;
end

%% Plot traces of all average interneurons w SEM
timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
y_green = mean_cell_green;
SEM_greens = std(mean(all_cells_trace_green,1))./sqrt(size(all_cells_trace_green,1));
 x_green = timestamp(1:length(y_green));
 
%x_green = 1:length(y_green);
%x_red = 1:length(y_red);
 
figure; hold on;
% shadedErrorBar(x_green,y_red,SEM_red,'lineprops','-r','transparent',1);hold on;
% shadedErrorBar(x_red,y_green,SEM_greens,'lineprops','-b','transparent',1);

shadedErrorBar(x_green,smooth((y_green),5),smooth((SEM_greens),5),'lineprops','-b','transparent',1);hold on;

if run_redcell==1
    y_red = mean_cell_red;
    SEM_red = std(mean(all_cells_trace_red,1))./sqrt(size(all_cells_trace_red,1));
    x_red = timestamp(1:length(y_red));
    shadedErrorBar(x_red,smooth((y_red),5),smooth((SEM_red),5),'lineprops','-r','transparent',1);
end
%% Plot mean traces of individual cells and determine if responsive

all_cells_trace_green = data.combined.trace.green_all;
all_cells_response_green = data.combined.response.green_all;
all_cell_peak_green = data.combined.peak.green_all;
all_cells_avg_around_peak_green = data.combined.avg_around_peak_green;
all_cell_STD_green = data.combined.STDbase.green_all;

if run_redcell==1
all_cells_trace_red = data.combined.trace.red_all;
all_cells_response_red = data.combined.response.red_all;
all_cell_peak_red = data.combined.peak.red_all;
all_cells_avg_around_peak_red = data.combined.avg_around_peak_red;
all_cell_STD_red = data.combined.STDbase.red_all;
end

std_level = 1.25;%set this here to change std

%find green responsive cells and plot results
for i=1:size(all_cells_response_green,1);
    a_green(i,:) = mean(all_cells_trace_green(i,:,:), 2); %mean of all trials for each cell
    SEM_a_green(i,:) = std(all_cells_trace_green(i,:,:))./sqrt(size(all_cells_trace_green(i,:,:),2));
    
    b_green(i,:) = mean(all_cells_response_green(i,:), 2); %average means responses (sound to 2s post sound) across trials for each cell
    c_green(i,:) = mean(all_cell_peak_green(i,:), 2); %average peak response
    d_green(i,:) = mean(all_cells_avg_around_peak_green(i,:), 2);
    e_green(i,:) = mean(all_cells_neg_avg_around_peak_green(i,:), 2); %average around negative peak
    
    %determine whether each cell is responsive (defined as average response
    %more than 2 STDS above baseline) 
    f = mean(all_cell_STD_green(i,:), 2); %average baseline STD across trials for each cell
    isResponsive_green(i) = d_green(i,:) > std_level*mean(f) & b_green(i,:)>0; %will be 0 or 1
    isNegResponsive_green(i) = e_green(i,:) < -std_level*mean(f) & b_green(i,:)<0;     
end

first_cell = 1;
last_cell = size(all_cells_response_green,1);
num_cells = last_cell-first_cell+1;
figure;

for i=first_cell:last_cell  
    %plot mean traces across cells with response means - responsive cells are
    %green, non-responsive are black 
    subplot_num = i-first_cell+1;
    subplot(ceil(sqrt(num_cells)),ceil(sqrt(num_cells)),subplot_num);
    x_green = 1:length(a_green(i,:));
    if isResponsive_green(i) == 1
        shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-b','transparent',1); hold on;
    else
        if isNegResponsive_green(i)== 1
            shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-c','transparent',1); hold on;
        else
            shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-k','transparent',1); hold on;
        end
    end
    
    plot(b_green(i),'o'); %plot average response
    plot(d_green(i),'o'); %plot average around peak
    plot(e_green(i),'o'); 
end  

clear first_cell last_cell num_cells d

%find red responsive cells and plot results
if run_redcell==1
for i=1:size(all_cells_response_red,1);
    a_red(i,:) = mean(all_cells_trace_red(i,:,:), 2); %mean of all trials for each cell
    SEM_a_red(i,:) = std(all_cells_trace_red(i,:,:))./sqrt(size(all_cells_trace_red(i,:,:),2));
    
    b_red(i,:) = mean(all_cells_response_red(i,:), 2); %average means responses (sound to 2s post sound) across trials for each cell
    c_red(i,:) = mean(all_cell_peak_red(i,:), 2); %average peak response
    d_red(i,:) = mean(all_cells_avg_around_peak_red(i,:), 2);
    e_red(i,:) = mean(all_cells_neg_avg_around_peak_red(i,:), 2); %average around negative peak
    
    %determine whether each cell is responsive (defined as average response
    %more than 2 STDS above baseline) 
    f = mean(all_cell_STD_red(i,:), 2); %average baseline STD across trials for each cell
    isResponsive_red(i) = d_red(i,:) > std_level*mean(f) & b_red(i,:)>0; %will be 0 or 1
    isNegResponsive_red(i) = e_red(i,:) < -std_level*mean(f)  & b_red(i,:)<0;     
end


first_cell = 1;
last_cell = size(all_cells_response_red,1);
num_cells = last_cell-first_cell+1;
figure;

for i=first_cell:last_cell  
    %plot mean traces across cells with response means - responsive cells are
    %red, non-responsive are black 
    subplot_num = i-first_cell+1;
    subplot(ceil(sqrt(num_cells)),ceil(sqrt(num_cells)),subplot_num);
    x_red = 1:length(a_red(i,:));
    if isResponsive_red(i) == 1
        shadedErrorBar(x_red,smooth((a_red(i,:)),10),smooth((SEM_a_red(i,:)),10),'lineprops','-r','transparent',1); hold on;
    else
        if isNegResponsive_red(i)== 1
            shadedErrorBar(x_red,smooth((a_red(i,:)),10),smooth((SEM_a_red(i,:)),10),'lineprops','-m','transparent',1); hold on;
        else
            shadedErrorBar(x_red,smooth((a_red(i,:)),10),smooth((SEM_a_red(i,:)),10),'lineprops','-k','transparent',1); hold on;
        end
    end
    
    plot(b_red(i),'o'); %plot average response
    plot(d_red(i),'o'); %plot average around peak
    plot(e_red(i),'o'); 
end  
end


data.combined.isResponsiveGreen=isResponsive_green';
data.combined.isNegResponsive_green = isNegResponsive_green';

if run_redcell==1
data.combined.isResponsiveRed=isResponsive_red';
data.combined.isNegResponsive_red = isNegResponsive_red';
end
% Plot traces of all Positive Responsive Cells  w SEM (for both red and green
% cells)
timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
responsive_cells_green = a_green(isResponsive_green,:,:);
mean_responsive_green = mean(responsive_cells_green, 1);
SEM_responsive_green = std(responsive_cells_green)./sqrt(size(responsive_cells_green,1));
x_green = timestamp(1:length(mean_responsive_green));

figure; hold on;
shadedErrorBar(x_green,smooth((mean_responsive_green),10),smooth((SEM_responsive_green),10),'lineprops','-b','transparent',1); hold on;

if run_redcell==1
responsive_cells_red = a_red(isResponsive_red,:,:);
mean_responsive_red = mean(responsive_cells_red, 1);
SEM_responsive_red = std(responsive_cells_red)./sqrt(size(responsive_cells_red,1));
 x_red = timestamp(1:length(mean_responsive_red));
 shadedErrorBar(x_red,smooth((mean_responsive_red),10),smooth((SEM_responsive_red),10),'lineprops','-r','transparent',1);
end

 



title('Mean of All Positively Responsive Cells');

% Plot traces of all Negative Responsive Cells  w SEM (for both red and green
% cells)
timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
neg_responsive_cells_green = a_green(isNegResponsive_green,:,:);
neg_mean_responsive_green = mean(neg_responsive_cells_green, 1);
neg_SEM_responsive_green = std(neg_responsive_cells_green)./sqrt(size(neg_responsive_cells_green,1));
 x_green = timestamp(1:length(neg_mean_responsive_green));
figure; hold on;
shadedErrorBar(x_green,smooth((neg_mean_responsive_green),10),smooth((neg_SEM_responsive_green),10),'lineprops','-b','transparent',1); hold on;

if run_redcell==1
neg_responsive_cells_red = a_red(isNegResponsive_red,:,:);
neg_mean_responsive_red = mean(neg_responsive_cells_red, 1);
neg_SEM_responsive_red = std(neg_responsive_cells_red)./sqrt(size(neg_responsive_cells_red,1));
 x_red = timestamp(1:length(neg_mean_responsive_red));
 shadedErrorBar(x_red,smooth((neg_mean_responsive_red),10),smooth((neg_SEM_responsive_red),10),'lineprops','-r','transparent',1);

end
title('Mean of All Negatively Responsive Cells');

%% Plot IsResponsive Percentage +/- Across Mice 
isResponsive_green = data.combined.isResponsiveGreen;
isNegResponsive_green = data.combined.isNegResponsive_green;
green_cell_count = 0;

if run_redcell==1
  isResponsive_red = data.combined.isResponsiveRed;
  isNegResponsive_red = data.combined.isNegResponsive_red;
  red_cell_count = 0;
end
figure; 

for a=1:length(mousename) %this will break when add multiple imaging blocks/mice, fix to average across blocks?
        mouseID=mousename{(a)}
        Imaging_Block=Imaging_sets(a,:);
        for i=1:length(Imaging_Block(i))
        Imaging_Num =  sprintf( '%03d', Imaging_Block(i))
        
        length_green_cell = length(data.([mouseID]).(['ImagingBlock' Imaging_Num]).nonredcell);
        green_cells = isResponsive_green(1+green_cell_count:green_cell_count+length_green_cell);
        green_cells_neg = isNegResponsive_green(1+green_cell_count:green_cell_count+length_green_cell);
        percent_isResponsive(1,a) = sum(green_cells)./length(green_cells);
        percent_NegResponsive(1,a) = sum(green_cells_neg)./length(green_cells);
        green_cell_count = green_cell_count + length_green_cell;
        
        %if there are redcells
        if run_redcell==1
            length_red_cell = length(data.([mouseID]).(['ImagingBlock' Imaging_Num]).redcell);
            red_cells = isResponsive_red(1+red_cell_count:red_cell_count+length_red_cell);
            red_cells_neg = isNegResponsive_red(1+red_cell_count:red_cell_count+length_red_cell);
            percent_isResponsive(2,a) = sum(red_cells)./length(red_cells);
            percent_NegResponsive(2,a) = sum(red_cells_neg)./length(red_cells);
            red_cell_count = red_cell_count + length_red_cell;
        end
            
        end
end
 data.([mouseID]).(['ImagingBlock' Imaging_Num]).percentIsResponsive = percent_isResponsive;
  data.([mouseID]).(['ImagingBlock' Imaging_Num]).percentNegResponsive = percent_NegResponsive;
 
 figure;
 x=[1 2];
 y=percent_isResponsive;

SEM = std(y')./sqrt(size(y,2));
errhigh = -SEM;
errlow  = SEM;
ylim([0 1]);

bar(x,mean(y'));     
ylim([0 1]);

hold on;
 
er = errorbar(x,mean(y'),errlow,errhigh);  hold on; 
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
plot(x,y','o'); hold off;
title('Percentage of All Positively Responsive Cells Across Mice, Green, Red');

 figure;
 x=[1 2];
 y=percent_NegResponsive;

SEM = std(y')./sqrt(size(y,2));
errhigh = -SEM;
errlow  = SEM;

bar(x,mean(y'));   
ylim([0 1]);

hold on;
 
er = errorbar(x,mean(y'),errlow,errhigh);  hold on; 
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
plot(x,y','o'); hold off;
ylim([0 1]);
title('Percentage of All Negatively Responsive Cells Across Mice, Green, Red');
 
%% Combine Locomotor Activity - this is getting moved to an earlier stage in analysis

all_cells_isLocoSound_green = []; 
all_cells_noLocoSound_green = []; 

if run_redcell==1
all_cells_isLocoSound_red = []; 
all_cells_noLocoSound_red = [];
end

for a=1:length(mousename)
for i=1:length(Imaging_Block(a))

       mouseID=mousename{(a)};
       Tosca_Session=Session{(a)};
       Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder
       folder = sprintf([path '/2P local data/' user '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]);
       cd(folder)
       Tosca_Run_number = num2str(Tosca_Runs(a));
       Image_Block = Imaging_Block(a);
       date=expt_date{(a)};
       
       %read tosca file
       Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder 
       folder = ([path '2P analysis/2P local data/' user '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]); %direct to specific Tosca folder within a 
       read_loco= [mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt'];
       
       folder = sprintf([path '2P local data/' user '/' mouseID '/' date '/' Image_Block]); %direct to specific Tosca folder within a 
       Imaging_Num =  sprintf( '%03d', Image_Block);
       filename = ['BOT_' mouseID '_noiseburst-' Imaging_Num '_Cycle00001_VoltageRecording_001.csv'];
       M = csvread(filename, 1,0);
       loco_data = dlmread(read_loco);%locomotor data  
       [loco_data,active_time] = locomotor_activity(loco_data,M);
       data.([mouseID]).(['ImagingBlock' Imaging_Num]).locomotion_data = loco_data;
       data.([mouseID]).(['ImagingBlock' Imaging_Num]).active_time = active_time;


     
       redcell = data.([mouseID]).(['ImagingBlock' Imaging_Num]).redcell;
       nonredcell = data.([mouseID]).(['ImagingBlock' Imaging_Num]).nonredcell;
       F7 = data.([mouseID]).(['ImagingBlock' Imaging_Num]).full_trace;
       
       mean_F7_green = squeeze(mean(F7(nonredcell,:),1));
       SEM_F7_green = std(F7(nonredcell,:),1)./sqrt(size(F7(nonredcell,:),1));
       
       mean_F7_red = squeeze(mean(F7(redcell,:),1));
       SEM_F7_red = std(F7(redcell,:),1)./sqrt(size(F7(redcell,:),1));
       
       timestamp = data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
       Sound_Time = data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;
       
       figure; 
       plot(timestamp,smooth(mean_F7_green,10),'-b','LineWidth',3); hold on;
       plot(timestamp,smooth(mean_F7_red,10),'r','LineWidth',3);hold on;
       shadedErrorBar(timestamp,smooth((mean_F7_green),10),smooth((SEM_F7_green),10),'lineprops','-b','transparent',1); hold on;
       shadedErrorBar(timestamp,smooth((mean_F7_red),10),smooth((SEM_F7_red),10),'lineprops','-r','transparent',1); hold on;
    %   plot(loco_data(:,1),loco_data(:,3)*10+mean(mean_F7_green)/3,'-k'); hold on; % loco data is 3 columns (timestamp, activity, active yes (1) or no (0))
       plot(loco_data(:,1),loco_data(:,3)*4,'-k'); hold on;
       %   h=vline(Sound_Time, 'b');
         hold on; 
         
         
        % Find sound times when animal is active
        for time=1:length(Sound_Time)
            sound = Sound_Time(time);
            window = Sound_Time(time)+2; % when is the sound?

            [c closest_frame_sound] = min(abs(loco_data(:,1)-sound));
            [c closest_frame_window] = min(abs(loco_data(:,1)-window));

            length_sound_window(time) = closest_frame_window-closest_frame_sound;
            
             if time>1
                difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
                closest_frame_window = closest_frame_window-difference_length_window(time);
             end
             
             closest_frame_sound_all(time) = closest_frame_sound;
             closest_frame_window_all(time) = closest_frame_window;
            
             
           %  isLocoSound(time) = sum(loco_data(closest_frame_sound:closest_frame_window,3))>0;
             isLocoSound(time) = sum(active_time(closest_frame_sound:closest_frame_window))>0;
        end
        figure; plot(loco_data(:,1),active_time); hold on;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).isLocoSound = isLocoSound;
           
        %now pull out only isLocoSound trials for red and green, average and
        %concatenate across mice
        trace_around_sound_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_green;%trace(dF/F) around sound for each non-VIP (cellxsoundxtime matrix)
        isLocoSound_green = mean(trace_around_sound_green(:,isLocoSound,:),2);
      %  size_isLocogreen=size(trace_around_sound_green(:,isLocoSound,:))
        all_cells_isLocoSound_green = [all_cells_isLocoSound_green; isLocoSound_green]; 
        
        noLocoSound_green = mean(trace_around_sound_green(:,isLocoSound==0,:),2);
      %  size_noLocogreen=size(trace_around_sound_green(:,isLocoSound==0,:))
        all_cells_noLocoSound_green = [all_cells_noLocoSound_green; noLocoSound_green]; 
        
        trace_around_sound_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_red;%trace(dF/F) around sound for each non-VIP (cellxsoundxtime matrix)
        isLocoSound_red = mean(trace_around_sound_red(:,isLocoSound,:),2);
        all_cells_isLocoSound_red = [all_cells_isLocoSound_red; isLocoSound_red]; 
        
        noLocoSound_red = mean(trace_around_sound_red(:,isLocoSound==0,:),2);
        all_cells_noLocoSound_red = [all_cells_noLocoSound_red; noLocoSound_red]; 
    
                                   
end
end
data.combined.trace.isLocoSoundgreen_all=all_cells_isLocoSound_green;
data.combined.trace.isLocoSoundred_all=all_cells_isLocoSound_red;

data.combined.trace.noLocoSoundgreen_all=all_cells_noLocoSound_green;
data.combined.trace.noLocoSoundred_all=all_cells_noLocoSound_red;

%% Plot traces of all cells with locomotion during sounds versus without
%green
mean_loco_green = mean(all_cells_isLocoSound_green,1);
SEM_loco_green = std(all_cells_isLocoSound_green,1)./sqrt(size(all_cells_isLocoSound_green,1));

mean_noloco_green = mean(all_cells_noLocoSound_green,1);
SEM_noloco_green = std(all_cells_noLocoSound_green,1)./sqrt(size(all_cells_noLocoSound_green,1));

%red
mean_loco_red = mean(all_cells_isLocoSound_red,1);
SEM_loco_red = std(all_cells_isLocoSound_red,1)./sqrt(size(all_cells_isLocoSound_red,1));

mean_noloco_red = mean(all_cells_noLocoSound_red,1);
SEM_noloco_red = std(all_cells_noLocoSound_red,1)./sqrt(size(all_cells_noLocoSound_red,1));

timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
x = squeeze(timestamp(1:length(mean_loco_green)));

  figure; 
       shadedErrorBar(x,smooth((mean_loco_green),5),smooth((SEM_loco_green),5),'lineprops','-b','transparent',1); hold on;
       shadedErrorBar(x,smooth((mean_noloco_green),5),smooth((SEM_noloco_green),5),'lineprops','-k','transparent',1); hold on;
       title('Effect of Locomotion - non VIP neurons, green is loco')
     %  legend('locomotion','no locomotion')
       
  figure;     
      
       shadedErrorBar(x,smooth((mean_loco_red),5),smooth((SEM_loco_red),5),'lineprops','-r','transparent',1); hold on;
       shadedErrorBar(x,smooth((mean_noloco_red),5),smooth((SEM_noloco_red),5),'lineprops','-k','transparent',1); 
       title('Effect of Locomotion - VIP neurons, red is loco')
       %legend('locomotion','no locomotion'); hold on;


        
%% individual_traces=data.combined.trace.green_all(:,:,:);
% for iii=1:size(individual_traces,2)
%     figure;
%     x=1:size(individual_traces,3);
%     y=squeeze(individual_traces(2,iii,:));
%     plot(x,y);
% end

     %% save data
folder = 'Z:\Carolyn\Figures for SfN 2019';
addpath(folder);
save('noiseburst_analysis.mat', 'data'); 

  