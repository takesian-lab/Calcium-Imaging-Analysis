function    [data,traces_R,traces_G]=Noiseburst_analysis(a,Frames,Frame_rate,Imaging_Block_String,Imaging_Num,mouseID,date,Sound_Time,...
        timestamp,i,analysis_folder,path_name,length_sound_trial_first,username,data);
    
if strcmp(setup.username, 'Carolyn')
    folder = sprintf([setup.path_name setup.username '/analyzed/' mouseID '/' date '/' analysis_folder '/suite2p/plane0'])
else
    folder = sprintf([setup.path_name setup.username '/' mouseID '/' analysis_folder])
end
cd(folder);
load(['Fall.mat']);
neuropil=Fneu; % all the neuropil traces (should be one for each trace)
cell=F; % all the cell fluorescence data
cell_id=find(iscell(:,1)==1); % finds the rows with actual cells (determined manually in GUI)
redcell_ones=find(redcell(1,:)==1);%which cells are red
idx=find(ismember(redcell_ones,cell_id));%index out red cells that are active
redcell=redcell_ones(idx);%which red cells are active
nonredcell=setdiff(cell_id,redcell);%what are the active non-red cells
F7=cell-0.7*Fneu;%Find neuropil corrected traces

%some of my data sets have double the frame rate, and this is used to get
%it down from 30 to 15 fps
if Frame_rate=='30'
    HELLO=1
    for ii=1:size(F7,1)%run through all the cells
        xx=reshape(F7(ii,:),2,[]);%take out every other frame
        yy=sum(xx,1)./size(xx,1);%average every two frames
        yyy(ii,:)=yy;%put back into a matrix of cells x frames
    end
    F7=yyy; %make it called F7 again so everything is the same
   
end
% 
% redcell=redcell-1;
% nonredcell=nonredcell-1;

%python to matlab correction
cell_number=size(F7,1);
F7=F7(:,(Frames(1,:)));


% %remove '0' values in the event cell#1 = iscell
% for k=1:length(nonredcell);
%     deleted=nonredcell==0;
%     nonredcell(deleted)=[];
% end
% 
% %remove '0' values in the event cell#1 = iscell
% for j=1:length(redcell);
%     deleted=redcell==0;
%     redcell(deleted)=[];
% end



% Pull out the sound traces to each noiseburst
count=0;
% figure;
%define sound window
for time=1:length(Sound_Time)
  %  count=count+1
    sound = Sound_Time(time);
    before = Sound_Time(time)-0.5; % half second before sound
    after = Sound_Time(time)+3; % full trace after sound
    window = Sound_Time(time)+2; % when is the sound?
    start_window = Sound_Time(time)+1;
    
    [c closest_frame_before] = min(abs(timestamp(:)-before));
    [c closest_frame_start_window] = min(abs(timestamp(:)-start_window));
    [c closest_frame_sound] = min(abs(timestamp(:)-sound));
    [c closest_frame_after] = min(abs(timestamp(:)-after));
    [c closest_frame_window] = min(abs(timestamp(:)-window));
    
    length_sound_trial(time) = closest_frame_after-closest_frame_before;
    length_sound_window(time) = closest_frame_window-closest_frame_sound;
    
    if time*a>1 
        difference_length_sound_trial(time) = length_sound_trial(time)-length_sound_trial_first;
        closest_frame_after = closest_frame_after-difference_length_sound_trial(time);
        difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
        closest_frame_window = closest_frame_window-difference_length_window(time);  
    else
        length_sound_trial_first = length_sound_trial(time);
    end
%     
%      if time>1 
%         difference_length_sound_trial(time) = length_sound_trial(time)-length_sound_trial(1);
%         closest_frame_after = closest_frame_after-difference_length_sound_trial(time);
%         difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
%         closest_frame_window = closest_frame_window-difference_length_window(time);  
% % % %     else
% % % %         length_sound_trial_first = length_sound_trial(time);
%     end
%     
   
    for k = 1:length(nonredcell)
        
        %pull out raw trace around sound
        trace_green = (F7(nonredcell(k),:))';%neuropil corrected traces
        raw_trace_around_sound_green(k,time,:) = (trace_green(closest_frame_before:closest_frame_after));
        
        %df/f
        baseline_mean_green = mean(trace_green(closest_frame_before:closest_frame_sound));%take mean of all green traces to get baseline
        trace_around_sound_green(k,time,:) = (bsxfun(@minus, raw_trace_around_sound_green(k,time,:),baseline_mean_green))./baseline_mean_green;
        
        %determine baseline mean and STD
        baseline_trace_green = trace_green(closest_frame_before:closest_frame_sound);
        baseline_df_over_f_green = (bsxfun(@minus, baseline_trace_green,baseline_mean_green))./baseline_mean_green;
        mean_baseline_green(k,time) = mean(baseline_df_over_f_green);
        std_baseline_green(k,time) = std(baseline_df_over_f_green);
        
        %determine average and peak responses to sound
        trace_around_soundwindow_green = trace_green(closest_frame_sound:closest_frame_window);
        trace_around_soundwindow_df_over_f_green = (bsxfun(@minus, trace_around_soundwindow_green,baseline_mean_green))./...
            baseline_mean_green;
        length_t = length(trace_around_soundwindow_df_over_f_green);
        [peak_sound_green(k,time) max_loc] = max(trace_around_soundwindow_df_over_f_green(2:length_t-1));
        max_loc=max_loc+1;
        [negpeak_sound_green(k,time) min_loc] = min(trace_around_soundwindow_df_over_f_green(2:length_t-1));
        min_loc=min_loc+1;
        avg_around_peak_green(k,time) =mean(trace_around_soundwindow_df_over_f_green(max_loc-1:max_loc+1));
        neg_avg_around_peak_green(k,time) =mean(trace_around_soundwindow_df_over_f_green(min_loc-1:min_loc+1));
        avg_sound_green(k,time) = mean(trace_around_soundwindow_df_over_f_green);
        

    end
    

        
    for k = 1:length(redcell)
        
        %pull out raw trace around sound
        trace_red = (F7(redcell(k),:))';%neuropil corrected traces
        raw_trace_around_sound_red(k,time,:) = (trace_red(closest_frame_before:closest_frame_after));
        
        %df/f
        baseline_mean_red = mean(trace_red(closest_frame_before:closest_frame_sound));%take mean of all red traces to get baseline
        trace_around_sound_red(k,time,:) = (bsxfun(@minus, raw_trace_around_sound_red(k,time,:),baseline_mean_red))./baseline_mean_red;
        
        %determine baseline mean and STD
        baseline_trace_red = trace_red(closest_frame_before:closest_frame_sound);
        baseline_df_over_f_red = (bsxfun(@minus, baseline_trace_red,baseline_mean_red))./baseline_mean_red;
        mean_baseline_red(k,time) = mean(baseline_df_over_f_red);
        std_baseline_red(k,time) = std(baseline_df_over_f_red);
        
        %determine average and peak responses to sound
        trace_around_soundwindow_red = trace_red(closest_frame_sound:closest_frame_window);
        trace_around_soundwindow_df_over_f_red = (bsxfun(@minus, trace_around_soundwindow_red,baseline_mean_red))./...
            baseline_mean_red;
        length_t = length(trace_around_soundwindow_df_over_f_red);
        [peak_sound_red(k,time) max_loc] = max(trace_around_soundwindow_df_over_f_red(2:length_t-1));
        max_loc=max_loc+1;
        [negpeak_sound_red(k,time) min_loc] = min(trace_around_soundwindow_df_over_f_red(2:length_t-1));
        min_loc=min_loc+1;
        avg_around_peak_red(k,time) =mean(trace_around_soundwindow_df_over_f_red(max_loc-1:max_loc+1));
        neg_avg_around_peak_red(k,time) =mean(trace_around_soundwindow_df_over_f_red(min_loc-1:min_loc+1));
        avg_sound_red(k,time) = mean(trace_around_soundwindow_df_over_f_red);
      
    end

      
end
traces_G=cat(3,traces_G,raw_trace_around_sound_green);
traces_R=cat(3,traces_R,raw_trace_around_sound_red);
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).raw_trace_around_sound_green = raw_trace_around_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).raw_trace_around_sound_red = raw_trace_around_sound_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).mean_baseline_green = mean_baseline_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).mean_baseline_red = mean_baseline_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).std_baseline_green = std_baseline_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).std_baseline_red = std_baseline_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_green = trace_around_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_red = trace_around_sound_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).peak_sound_green = peak_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).peak_sound_red = peak_sound_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_around_peak_green = avg_around_peak_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_around_peak_red = avg_around_peak_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).neg_avg_around_peak_green = neg_avg_around_peak_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).neg_avg_around_peak_red = neg_avg_around_peak_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_sound_green = avg_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_sound_red = avg_sound_red;
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).redcell = redcell;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).nonredcell = nonredcell;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).full_trace = F7;
end

    
