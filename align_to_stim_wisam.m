function    [block]=align_to_stim_wisam(block)
% [block]=align_to_stim_wisam(block)
%
% DOCUMENTATION IN PROGRESS
% 
% This function is applied to a block and loops over trials and ROIs
% And returns the block data windowed around the stimulus for all trials,
% and ROIs
%
% Is the name of this function misleading?
% This function also preforms neuropil correction on the calcium traces.
%
% This function aligns the observation window to the stimulus offset.
% Currently the window is hardcoded to -0.5 to 2.5 seconds after the
% stimulus onset. (TODO: Fix this!)
% 
% Argument(s):
%   block (struct) 
%
% Returns:
%   block -  size(block) =  [numROIs, numTrials/block, window size (frames)]
% 
% Notes:
%
% Needed from block:
%   Sound_Time (vector): from define_sound 
%                        length(Sound_Time) = number of trials in the block
%   timestamp - from define_sound
%   F - from define_suite2p
%   Fneu - from define_suite2p
%
% TODO: What is the difference between 'block.timestamp' and 'block.Sound_Time'?
% TODO: There are magic numbers in this file

% Check that we have the required directories
% We need: 1) Path to Info file, and 2) Path to the blocks
% If we do not have the files, skip this step entirely.
if ismissing(block.setup.suite2p_path) || ismissing(block.setup.block_path)
    disp('Skipping align to stim...');
    return
end

disp('Aligning to stim...');

%% Window Parameters

% TODO: Make these a parameters of the function
time_before_stimulus_onset = 0.5;   % seconds
time_after_stimulus_onset = 2.5;    % seconds

%% Pull out the sound traces to each noiseburst

% Grab the onset times
Sound_Time = block.Sound_Time;

% Loop over stimulus onset times (trials)
% length(Sound_Time) = number of trials in the block
%
% TODO: Memory Allocation
% TODO: What is the purpose of this loop?
for timestamp_stim_onset=1:length(Sound_Time)
    
    % Grab timestamps for the current window
    % These are timestamps from the full time trace (the whole block)
    % [Timestamp] Stimulus onset
    sound           = Sound_Time(timestamp_stim_onset);
    % [Timestamp] Beginning of the window (before the stimulus onset)
    before          = Sound_Time(timestamp_stim_onset) - time_before_stimulus_onset;
    % [Timestamp] End of the window (after the stimulus onset)
    after           = Sound_Time(timestamp_stim_onset) + time_after_stimulus_onset;
    % TODO: why is this different than 'sound'?
    window          = Sound_Time(timestamp_stim_onset) + 1; % when is the sound? 
    % TODO: Why are we adding 1?
    start_window    = Sound_Time(timestamp_stim_onset) + 1;
    
    % TODO: Figure out what's happening here
    % TODO: Memory Allocation
    [c closest_frame_before]        = min(abs(block.timestamp(:)-before));
    [c closest_frame_start_window]  = min(abs(block.timestamp(:)-start_window));
    [c closest_frame_sound]         = min(abs(block.timestamp(:)-sound));
    [c closest_frame_after]         = min(abs(block.timestamp(:)-after));
    [c closest_frame_window]        = min(abs(block.timestamp(:)-window));
    
    % TODO: Memory Allocation
    length_sound_trial(timestamp_stim_onset) = closest_frame_after-closest_frame_before;
    length_sound_window(timestamp_stim_onset) = closest_frame_window-closest_frame_sound;
    
    % TODO: Figure out what's happening here
    % TODO: Memory Allocation
    % if time*a>1
    if timestamp_stim_onset>1
        
        difference_length_sound_trial(timestamp_stim_onset) = length_sound_trial(timestamp_stim_onset)-length_sound_trial_first;
        closest_frame_after = closest_frame_after-difference_length_sound_trial(timestamp_stim_onset);
        difference_length_window(timestamp_stim_onset) = length_sound_window(timestamp_stim_onset)-length_sound_window(1);
        closest_frame_window = closest_frame_window-difference_length_window(timestamp_stim_onset);
        
    else% First time through the outer loop
        % TODO: Is this the time to the first trial?
    	length_sound_trial_first = length_sound_trial(timestamp_stim_onset);
    end
    
    % NOTE: neuropil correction step in a for loop 
    % For a given stimulus onset time (trial)
    % Loop over the ROIs
    for k = 1:size(block.F,1)
        % Generate a neuropil corrected trace
        % TODO: Magic Number: why is this 0.7?
        % length(F7) = number of frames in the entire block
        F7 = block.F(k,:) - 0.7*block.Fneu(k,:);
        % Create windowed 'raw' traces aligned to a stimulus
        % These windows are extracted from the entire block of 'neuropil corrected' 
        % calcium traces
        raw_trace(k,timestamp_stim_onset,:) = (F7(closest_frame_before:closest_frame_after));
        % Mean Baseline activity (dF/F) 
        % Average Onset of the window to the onset of the stimulus
        baseline_mean = mean(F7(closest_frame_before:closest_frame_sound));
        % mean subtracted and normalized by the mean
        % TODO: Is this so we can easily detect positive and negative peaks?
        % TODO: What if there is a big variance in the baseline mean?
        % i.e. An ROI with a high mean baseline will have a smaller dF/F
        trace_around_stim(k,timestamp_stim_onset,:) = (bsxfun(@minus, raw_trace(k,timestamp_stim_onset,:),baseline_mean))./baseline_mean;
    end
    
end

% TODO: Now the window appears to be longer than 30 Hz sampling 
block.aligned_to_stim = trace_around_stim;

end


