function    [data]=align_to_stim(save_path,setup);
block_path = save_path;
cd(block_path)
allfiles=dir('*Block*');

% 
% for a=1:length(setup.mousename)
%     mouseID=setup.mousename{a,1};
%     Imaging_Block=setup.Imaging_sets{a,1};
    
    %make matrices to concatenate across imaging blocks
%     traces_G=[];
%     peak_G=[];
%     peak_maxG=[];
%     peak_minG=[];
%     baseline=[];
%     raw_trace_G=[];
%     std_baseline=[];
%     response=[];
    
%    for i=1:length(Imaging_Block)
%         Frame_rate=setup.framerate{a,1}(i);
%         date=char(setup.expt_date{a,1}(i));
%         analysis_folder=char(setup.analysis_name{a,1}(i));
%         if strcmp(setup.username, 'Carolyn')
%             folder = sprintf([setup.path_name setup.username '/analyzed/' mouseID '/' date '/' analysis_folder '/suite2p/plane0'])
%         else
%             folder = sprintf([setup.path_name setup.username '/' mouseID '/' analysis_folder])
%         end
%         cd(folder);
%         load(['Fall.mat']);
    
%         Frames=cell2mat(setup.Frame_set{a,1}(i));
%         Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
%         Imaging_Block_String = num2str(Imaging_Block(i));
%         timestamp =  data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
%         Sound_Time = data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;
% 
%         
%        neuropil=Fneu; % all the neuropil traces (should be one for each trace)
%         cell=F; % all the cell fluorescence data
%         cell_id=find(iscell(:,1)==1); % finds the rows with actual cells (determined manually in GUI)
%         nonredcell=cell_id;%what are the active non-red cells
%         F7=cell-0.7*Fneu;%Find neuropil corrected traces
%         length(Frames)
        
%         %python to matlab correction
%         cell_number=size(F7,1)
%         F7=F7(1:cell_number,(Frames(1,:)));
%         
        % Pull out the sound traces to each noiseburst
        %define sound window
        Sound_Time = block.Sound_Time;
        for time=1:length(Sound_Time)
            sound = Sound_Time(time)
            before = Sound_Time(time)-0.5; % half second before sound
            after = Sound_Time(time)+2; % full trace after sound
            window = Sound_Time(time)+1; % when is the sound?
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
    
            for k = 1:size(block.F,1)
                %pull out raw trace around stim
                trace_block = (F7(nonredcell(k),:))';%neuropil corrected traces
                raw_trace_around_sound_green(k,time,:) = (trace_block(closest_frame_before:closest_frame_after));
                
                %df/f
                baseline_mean_green = mean(trace_block(closest_frame_before:closest_frame_sound));%take mean of all green traces to get baseline
                trace_around_sound_green(k,time,:) = (bsxfun(@minus, raw_trace_around_sound_green(k,time,:),baseline_mean_green))./baseline_mean_green;
                
                %determine baseline mean and STD
                baseline_trace_green = trace_block(closest_frame_before:closest_frame_sound);
                baseline_df_over_f_green = (bsxfun(@minus, baseline_trace_green,baseline_mean_green))./baseline_mean_green;
                mean_baseline_green(k,time) = mean(baseline_df_over_f_green);
                std_baseline_green(k,time) = std(baseline_df_over_f_green);
                
                %determine average and peak responses to sound
                trace_around_soundwindow_green = trace_block(closest_frame_sound:closest_frame_window);
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
        end
        traces_G=cat(2,traces_G,trace_around_sound_green);
        peak_G=[peak_G,peak_sound_green];
        peak_maxG=[peak_maxG,avg_around_peak_green];
        peak_minG=[peak_minG,neg_avg_around_peak_green];
        baseline=cat(2,baseline,mean_baseline_green);
        raw_trace_G=cat(2,raw_trace_G,raw_trace_around_sound_green);
        std_baseline=[std_baseline,std_baseline_green];
        response=[response,avg_sound_green];
       
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).raw_trace_around_sound_green = raw_trace_around_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).mean_baseline_green = mean_baseline_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).std_baseline_green = std_baseline_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_green = trace_around_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).peak_sound_green = peak_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_around_peak_green = avg_around_peak_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).neg_avg_around_peak_green = neg_avg_around_peak_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).avg_sound_green = avg_sound_green;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).nonredcell = nonredcell;
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).full_trace = F7;
    
        clear trace_around_sound_green peak_sound_green neg_avg_around_peak_green mean_baseline_green raw_trace_around_sound_green std_baseline_green avg_sound_green
   end
    
        data.([mouseID]).traces_G=traces_G;
        data.([mouseID]).peak_G=peak_G;
        data.([mouseID]).peak_maxG=peak_maxG;
        data.([mouseID]).peak_minG=peak_minG;
        data.([mouseID]).baseline=baseline;
        data.([mouseID]).raw_trace_G=raw_trace_G;
        data.([mouseID]).std_baseline=std_baseline;
        data.([mouseID]).response=response;
   
   
   
end
end

