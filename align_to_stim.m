function [block] = align_to_stim(block)
% DOCUMENTATION IN PROGRESS
% 
% This function pulls out the trial windows for each sound presentation
% and stores in block
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
%
% Variables needed from block.setup:
% -Sound_Time - from define_sound
% -timestamp - from define_sound
% -F - from define_suite2p
% -Fneu - from define_suite2p
%
%
% TODO: Magic numbers, also pull out neuropil and spike traces
% Search 'TODO'

%%  Skip this function if Suite2p and Bruker data are not available

if ismissing(block.setup.suite2p_path) || ismissing(block.setup.block_path)
    disp('Skipping align to stim...');
    return
end

disp('Aligning to stim...');

%% Pull out the sound traces to each noiseburst
%define sound window
Sound_Time = block.Sound_Time;
for time=1:length(Sound_Time)
    sound = Sound_Time(time);
    before = Sound_Time(time)-0.5; % half second before sound
    after = Sound_Time(time)+2.5; % full trace after sound
    window = Sound_Time(time)+1; % when is the sound?
    start_window = Sound_Time(time)+1;
    
    
    [c closest_frame_before] = min(abs(block.timestamp(:)-before));
    [c closest_frame_start_window] = min(abs(block.timestamp(:)-start_window));
    [c closest_frame_sound] = min(abs(block.timestamp(:)-sound));
    [c closest_frame_after] = min(abs(block.timestamp(:)-after));
    [c closest_frame_window] = min(abs(block.timestamp(:)-window));
    
    
    length_sound_trial(time) = closest_frame_after-closest_frame_before;
    length_sound_window(time) = closest_frame_window-closest_frame_sound;
    
    %     if time*a>1
    if time>1
        difference_length_sound_trial(time) = length_sound_trial(time)-length_sound_trial_first;
        closest_frame_after = closest_frame_after-difference_length_sound_trial(time);
        difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
        closest_frame_window = closest_frame_window-difference_length_window(time);
        
    else
        length_sound_trial_first = length_sound_trial(time);
    end
    
    for k = 1:size(block.F,1)
        %generate a neuropil corrected trace 
        F7 = block.F(k,:) - 0.7*block.Fneu(k,:);
        % pull out the frames aligned to a stim (0.5s before to 3s after)
        raw_trace(k,time,:) = (F7(closest_frame_before:closest_frame_after));
        %df/f
        baseline_mean = mean(F7(closest_frame_before:closest_frame_sound));
        trace_around_stim(k,time,:) = (bsxfun(@minus, raw_trace(k,time,:),baseline_mean))./baseline_mean;
    end
end
 block.aligned_to_stim = trace_around_stim;
end


