function [block] = align_to_stim(block, sound_time_delay)
% DOCUMENTATION IN PROGRESS
% 
% This function pulls out the trial windows for each sound presentation
% and stores in block
% 
% Argument(s): 
%   block (struct)
%   m (struct)
% 
% Returns: 
%   block (struct)
% 
% Notes:
%
%
% 
%
% Variables needed from block.setup:
% -Sound_Time - from define_sound
% -timestamp - from define_sound
% -F - from define_suite2p
% -Fneu - from define_suite2p
% -spks - from define_suite2p
%
%
% TODO:
% Search 'TODO'

%%  Skip this function if Suite2p and Bruker data are not available

if ismissing(block.setup.suite2p_path) || ismissing(block.setup.block_path)
    disp('Skipping align to stim...');
    return
end

disp('Aligning to stim...');


%define sound window
setup = block.setup;
Sound_Time = block.Sound_Time;

%Used for Maryse's behavior stim
if nargin == 2
    Sound_Time = Sound_Time + sound_time_delay;
    block.setup.constant.after_stim = 3;
    setup.constant.after_stim = 3;
end

%Accomodate multiplane data
if isfield(block, 'MultiplaneData')
    multiplaneData = true;
    nPlanes = setup.XML.nPlanes;
else
    multiplaneData = false;
    nPlanes = 1;
end

%% ALIGN DATA

for n = 1:nPlanes
    plane = nPlanes - 1;
    
    if multiplaneData
        planeName = strcat('plane', num2str(n - 1));
        F = block.F.(planeName);
        Fneu = block.Fneu.(planeName);
        spks = block.spks.(planeName);
        timestamp = block.timestamp.(planeName);
    else
        F = block.F;
        Fneu = block.Fneu;
        spks = block.spks;
        timestamp = block.timestamp;
    end
    
    %generate a neuropil corrected trace  = raw F - (neuropil coefficient)*neuropil signal
    F7 = F - setup.constant.neucoeff*Fneu;

    %Preallocate variables and trial length in frames
    %Align everything to closest frame to sound start
    baseline_inFrames = round(setup.constant.baseline_length*block.setup.framerate);
    after_inFrames = round(setup.constant.after_stim*block.setup.framerate);
    duration_inFrames = baseline_inFrames + after_inFrames;

    nanMat = nan(size(F,1),length(Sound_Time),duration_inFrames);
    F7_stim = nanMat;
    F_stim = nanMat;
    Fneu_stim = nanMat;
    spks_stim = nanMat;

    % loop through each stim-presenation
    for time=1:length(Sound_Time)

        sound = Sound_Time(time);
        [~, closest_frame_sound] = min(abs(timestamp(:)-sound));
        A = closest_frame_sound - baseline_inFrames;
        B = closest_frame_sound + after_inFrames - 1;
        a = 1;
        b = duration_inFrames;

        % loop through each "iscell" to find the stim-aligned 1) raw
        % fluoresence 2) neuropil signal 3) neuropil-corrected floresence 4)
        % df/F for the neuropil corrected fluoresence 5) deconvolved spikes

        %If user-defined baseline is before the beginning of the block
        %recording, set A = 1 and the beginning of the block will be nan
        if A < 1
            a = abs(A) + 2;
            A = 1;
        end

        %If user-defined trial is longer than block recording, take portion
        %of trial up to the end of recording, the rest of the frames will be nan
        if B > size(F7,2)
            B = size(F7,2);
            b = length(A:B);
        end
        
        %Catch problems
        if A > B
            error('A should not be greater than B. Check timestamps')
        end
        
        % pull out the frames aligned to a stim (defined in frames)
        F7_stim(:,time,a:b) = F7(:,A:B);
        F_stim(:,time,a:b) =  F(:,A:B);
        Fneu_stim(:,time,a:b) = Fneu(:,A:B);
        spks_stim(:,time,a:b) = spks(:,A:B);
    end
    
    % SAVE ALIGNED DATA TO BLOCK
    if multiplaneData
         block.aligned_stim.F7_stim.(planeName) = F7_stim;
         block.aligned_stim.F_stim.(planeName) = F_stim;
         block.aligned_stim.Fneu_stim.(planeName) = Fneu_stim;
         block.aligned_stim.spks_stim.(planeName) = spks_stim;
    else
         block.aligned_stim.F7_stim = F7_stim;
         block.aligned_stim.F_stim = F_stim;
         block.aligned_stim.Fneu_stim = Fneu_stim;
         block.aligned_stim.spks_stim = spks_stim;
    end
end

 
end

%% PREVIOUS CODE


% % loop through each stim-presenation
% for time=1:length(Sound_Time)
%    
%     sound = Sound_Time(time);
%     
%     % define when baseline starts
%     before = Sound_Time(time)-setup.constant.baseline_length;
%     
%     % define the end time of the trace
%     after = Sound_Time(time)+setup.constant.after_stim;
%     
%     %define the "response window" (stim presentation to end of expected
%     %trace)
%     window = Sound_Time(time)+setup.constant.response_window; 
% 
%     
%     % find the frames that correspond to the above times. This will give
%     % you four frame numbers per stim. 1) start baseline 2) stim 3) end of
%     % the response window 4) end of the trace
%     [c closest_frame_before] = min(abs(block.timestamp(:)-before));
%     [c closest_frame_sound] = min(abs(block.timestamp(:)-sound));
%     [c closest_frame_after] = min(abs(block.timestamp(:)-after));
%     [c closest_frame_window] = min(abs(block.timestamp(:)-window));
%     
%     % how long is the entire trace (in frames)?
%     length_sound_trial(time) = closest_frame_after-closest_frame_before;
%     % how long is the response window (in frames)?
%     length_sound_window(time) = closest_frame_window-closest_frame_sound;
%     
%    % sometimes, the length of the trials are off by one frame due to how we
%    % calculate the frame times. Below, we correct for this by making every
%    % trace the same as trial 1
%     if time>1
%         difference_length_sound_trial(time) = length_sound_trial(time)-length_sound_trial_first;
%         closest_frame_after = closest_frame_after-difference_length_sound_trial(time);
%         difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
%         closest_frame_window = closest_frame_window-difference_length_window(time);
%         
%     else
%         length_sound_trial_first = length_sound_trial(time);
%     end
%     
%     % loop through each "iscell" to find the stim-aligned 1) raw
%     % fluoresence 2) neuropil signal 3) neuropil-corrected floresence 4)
%     % df/F for the neuropil corrected fluoresence 5) deconvolved spikes
%     for k = 1:size(block.F,1)
%         %generate a neuropil corrected trace 
%         % neuropil corrected trace = raw F - (neuropil
%         % coefficient)*neuropil signal
%         F7 = block.F(k,:) - setup.constant.neucoeff*block.Fneu(k,:);
%         
%         % pull out the frames aligned to a stim (defined in frames)
%         if closest_frame_after > size(F7,2)
%             closest_frame_after = size(F7,2);
%         end
%         F7_stim(k,time,:) = (F7(closest_frame_before:closest_frame_after));
%         F_stim(k,time,:) =  block.F(k,closest_frame_before:closest_frame_after);
%         Fneu_stim(k,time,:) = block.Fneu(k,closest_frame_before:closest_frame_after);
%         spks_stim(k,time,:) = block.spks(k,closest_frame_before:closest_frame_after);
% 
%         %df/f
% %         baseline_mean = mean(F7(closest_frame_before:closest_frame_sound));
% %         F7_dfF_stim(k,time,:) = (bsxfun(@minus, raw_trace(k,time,:),baseline_mean))./baseline_mean;
%     end
% end
% %  block.aligned_stim.F7_dfF_stim = F7_dfF_stim;
%  block.aligned_stim.F7_stim = F7_stim;
%  block.aligned_stim.F_stim = F_stim;
%  block.aligned_stim.Fneu_stim = Fneu_stim;
%  block.aligned_stim.spks_stim = spks_stim;
%  


