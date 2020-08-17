function [block] = define_loco_singleblock(block)
% DOCUMENTATION IN PROGRESS
% 
% This function defines trials when the animal is active or not.
% This requires both Tosca and Bruker data.
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
% -block.Sound_Time
% -block.locomotion_data
% -block.active_timef
%
% TODO: Remove magic numbers 
% Search 'TODO'

%% Skip this function if both Tosca and Bruker data are not available

if ismissing(block.setup.Tosca_path)
    disp('Skipping define_loco...');
    return
elseif ismissing(block.setup.block_path) && ismissing(block.setup.VR_path)
    disp('Skipping define_loco...');
    return
end

disp('Finding active trials')

setup = block.setup;
Sound_Time = block.Sound_Time;
loco_data = block.locomotion_data; %TRANSFORMED LOCO DATA
active_time = block.active_time;
trial_data = block.loc_Trial_times;

%% Find sound times when animal is active
for time=1:length(Sound_Time)
    sound = Sound_Time(time);% when is the sound?
    window = Sound_Time(time)+setup.constant.locowindow; 

    [c closest_loc_sound] = min(abs(trial_data{time}(:)-sound));
    [c closest_loc_window] = min(abs(trial_data{time}(:)-window));

    length_sound_window(time) = closest_loc_window-closest_loc_sound;

     if time>1
        difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
        closest_loc_window = closest_loc_window-difference_length_window(time);
     end

     closest_loc_sound_all(time) = closest_loc_sound;
     closest_loc_window_all(time) = closest_loc_window;
     
     % if any times are active, it will count as loco time
    
     %For the case where the closest frame surpasses the length of active_time
     %This should only be the case if we're on the last trial, so let the
     %code throw an error if we're not (no else case)
     if closest_loc_window > length(active_time) && time == length(Sound_Time)
         closest_loc_window = length(active_time);
     end
    
     Loco_1(time) = sum(active_time(closest_loc_sound:closest_loc_window))>1;
end

%% Plot locomotor activity
% figure;
% 
% subplot(2,1,1); hold on
% title('Locomotor activity')
% ylabel('Activity')
% plot(loco_data(:,1), loco_data(:,3));
% 
% subplot(2,1,2); hold on
% ylabel('Considered active')
% xlabel('Seconds')
% plot(loco_data(:,1), active_time > 0); hold on;

%% Save results
block.isLoco = Loco_1;
block.setup = setup;
end