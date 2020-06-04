function    [data]=Noiseburst_analysis_greenonly_wisam(data,setup)
% [data]=Noiseburst_analysis_greenonly_wisam(data,setup)
%
% DOCUMENTATION IN PROGRESS
% 
% This function performs noiseburst analysis across blocks,
% mice, and ROIs
% 
% Argument(s): 
%   data (struct)
%   setup (struct)
% 
% Returns:
%   data (struct)
% 
% Notes:
% 
% TODO: Memory allocation needed
% TODO: Magic number
% Search 'TODO'

% Loop over the mice
for a=1:size(setup.mousename,1)
    
    % For a given mouse
    % Loop over the ROIs
	for b=1:size(setup.mousename,2)
    
        % Catch edge case: no ROIs
        if isempty(setup.mousename{a,b})
            continue;
        end
        
        % Grab current mouse ID
        mouseID=setup.mousename{a,b};
        % Grab current Imaging block(s)
        Imaging_Block=setup.Imaging_sets{a,b};
    
        % Make matrices to concatenate across imaging blocks
        % TODO: Memory Allocation
        traces_G=[];
        peak_G=[];
        peak_maxG=[];
        peak_minG=[];
        baseline=[];
        raw_trace_G=[];
        std_baseline=[];
        response=[];
    
        % For a given mouse and ROI
        % Loop over all blocks
        for i=1:length(Imaging_Block)
            
            % Grab unique block string
            unique_block_name = setup.unique_block_names{a,b}(i);
            % Grab block data (struct)
            block = data.([mouseID]).([unique_block_name]);

            % TODO: is this something we could check for once?
            % Check if there is a 'Fall' file
            % return from the function if we don't
            if ismissing(block.setup.suite2p_path)
                disp('Skipping Suite2p data for block...');
                disp(unique_block_name);
                return
            end

            % Grab all time stamps for the entire block
            timestamp =  block.timestamp;
            % Grab all stimulus onset time stamps for the entire block
            % length(Sound_Time) = number of trials in the current block
            Sound_Time = block.Sound_Time;

            % TODO: Get rid of this?
            % Frame_rate = setup.framerate{a,1}(i); % Not used in this script
            
            % Grab all the neuropil traces (there should be one for each ROI)
            % Dim: numROIs : numFrames (in the current block)
            Fneu = block.Fneu; 
            % Grab all the cell fluorescence data (one for each ROI)
            % Dim: numROIs : numFrames (in the current block)
            cell = block.F; 
            
            % Find the rows with actual cells (determined manually in GUI)
            cell_id = find(block.iscell(:,1)==1); 
            % What are the active non-red cells
            nonredcell = cell_id;
            % TODO: Magic number
            % Find neuropil corrected traces
            F7 = cell-0.7*Fneu;

            % TODO: Get rid of this?
            % % python to matlab correction % NO LONGER NEEDED
            % cell_number = size(F7,1);
            % F7 = F7(1:cell_number,:);

            % Pull out the sound traces to each noiseburst define sound window
            for time=1:length(Sound_Time)
                
                % TODO: Magic numbers
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
                    trace_green = (F7(k,:))';%neuropil corrected traces
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
            end
            traces_G=cat(2,traces_G,trace_around_sound_green);
            peak_G=[peak_G,peak_sound_green];
            peak_maxG=[peak_maxG,avg_around_peak_green];
            peak_minG=[peak_minG,neg_avg_around_peak_green];
            baseline=cat(2,baseline,mean_baseline_green);
            raw_trace_G=cat(2,raw_trace_G,raw_trace_around_sound_green);
            std_baseline=[std_baseline,std_baseline_green];
            response=[response,avg_sound_green];

            block.raw_trace_around_sound_green = raw_trace_around_sound_green;
            block.mean_baseline_green = mean_baseline_green;
            block.std_baseline_green = std_baseline_green;
            block.trace_around_sound_green = trace_around_sound_green;
            block.peak_sound_green = peak_sound_green;
            block.avg_around_peak_green = avg_around_peak_green;
            block.neg_avg_around_peak_green = neg_avg_around_peak_green;
            block.avg_sound_green = avg_sound_green;
            block.nonredcell = nonredcell;
            block.full_trace = F7;
            data.([mouseID]).([unique_block_name]) = block;

            clear trace_around_sound_green peak_sound_green neg_avg_around_peak_green mean_baseline_green raw_trace_around_sound_green std_baseline_green avg_sound_green
        end
        
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

