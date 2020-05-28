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

%% Find sound times when animal is active
for time=1:length(Sound_Time)
    sound = Sound_Time(time);% when is the sound?
    window = setup.constant.locowindow; 

    [c closest_frame_sound] = min(abs(loco_data(:,1)-sound));
    [c closest_frame_window] = min(abs(loco_data(:,1)-window));

    length_sound_window(time) = closest_frame_window-closest_frame_sound;

     if time>1
        difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
        closest_frame_window = closest_frame_window-difference_length_window(time);
     end

     closest_frame_sound_all(time) = closest_frame_sound;
     closest_frame_window_all(time) = closest_frame_window;


    % if any times are active, it will count as loco time
     Loco_1(time) = sum(active_time(closest_frame_sound:closest_frame_window))>1;
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