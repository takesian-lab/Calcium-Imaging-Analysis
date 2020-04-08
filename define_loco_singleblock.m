function [block] = define_loco_singleblock(block)

disp('Finding active trials')

%Needed from setup:
%VR_path
%VR_filename
%stim_protocol
%voltage_recording

setup = block.setup;
Sound_Time = block.Sound_Time;
loco_data = block.loco_data; %RAW LOCO DATA

%% Prepare to read VoltageRecording file in locomotor_activity

isLoco = []; 
       
if setup.stim_protocol == 4 || setup.voltage_recording == 0
    cd(setup.VR_path) %Widefield
else
    cd(setup.block_path) %2p
end
filename = setup.VR_filename;

[loco_data,active_time] = locomotor_activity(loco_data,filename);
block.locomotion_data = loco_data; %TRANSFORMED LOCO DATA
block.active_time = active_time;


%% Find sound times when animal is active
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