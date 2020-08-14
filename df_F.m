function    [data]=df_F(data);
%
% DOCUMENTATION IN PROGRESS
%
% What does this function do?
% This function takes the concatenated, stim-aligned data and calculates
% dF/Fo based upon the local baseline
%
% Argument(s):
%   data(struct)
%
% Returns:
%   data(struct)
%
% Notes:
%
%
% TODO: Magic numbers, does this work with multiple ROIs?

% Search 'TODO'
%% 

setup = data.setup;
for a=1:size(setup.mousename,1) %Mice
    for b=1:size(setup.mousename,2) %ROIs
        
        if isempty(setup.mousename{a,b})
            continue;
        end
        
       
        mouseID = setup.mousename{a,b};
        
     % concatenated, neuropil-corrected, stim-aligned trace
        F7 = data.([mouseID]).cat.F_cat; 
        unique_block_name = setup.unique_block_names{a,b}(1);
     % settings for trace window
        constants = data.([mouseID]).([unique_block_name]).setup.constant;
        FrameRate = unique(setup.FrameRate{a,b});
        
        if length(FrameRate) > 1
            error('Not all frame rates are the same.')
            %TO DO: Make sure code accounts for different frame rates.
        end

        base_numFrames = round(FrameRate*constants.baseline_length);
        window_frames = (FrameRate*constants.response_window);
        local_base = zeros(size(F7,1),size(F7,2),base_numFrames); 
        
        for i = 1:size(F7,1)
            for j = 1:(size(F7,2))
               local_base(i,j,:) = F7(i,j,1:base_numFrames); %baseline
               mean_locBase = mean(local_base(i,j,:)); %mean local baseline
               stim_trace(i,j,:) = (bsxfun(@minus, F7(i,j,:),mean_locBase))./mean_locBase; % dF/F
               df_base(i,j,:) = stim_trace(i,j,1:base_numFrames);
               
                 %determine average and peak responses to sound
                trace_around_stimwindow(i,j,:) = stim_trace(i,j,base_numFrames+1:window_frames+base_numFrames);
%                 trace_around_soundwindow_df_over_f_green = (bsxfun(@minus, trace_around_soundwindow_green,baseline_mean_green))./...
%                     baseline_mean_green;
                length_t = length(trace_around_stimwindow);
                [peak_stim(i,j) max_loc] = max(trace_around_stimwindow(2:length_t-1));
                max_loc=max_loc+1;
                [negpeak_stim(i,j) min_loc] = min(trace_around_stimwindow(2:length_t-1));
                min_loc=min_loc+1;
                avg_around_peak(i,j) =mean(trace_around_stimwindow(max_loc-1:max_loc+1));
                avg_around_negpeak(i,j) =mean(trace_around_stimwindow(min_loc-1:min_loc+1));
%                 avg_sound_green(i,j) = mean(trace_around_stimwindow);

            end
            
        end
           data.([mouseID]).stim_df_f.baseline = df_base;
           data.([mouseID]).stim_df_f.F7_df_f = stim_trace;
           data.([mouseID]).stim_df_f.max_peak = avg_around_peak;
           data.([mouseID]).stim_df_f.min_peak = avg_around_negpeak;
           data.([mouseID]).stim_df_f.peak_val = peak_stim;
           data.([mouseID]).stim_df_f.window_trace = trace_around_stimwindow;
            
    end
        display(['...stim-aligned DF/F analyzed...'])
    end

