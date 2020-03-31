
[loco_activity,isLocoSound] = locomotor_frames(abslocdat,Sound_Time 

% when is the animal moving?
% active_time=zeros(i,1);
for i=1:length(abslocdat)
    if abslocdat>0.8
        active_time(i) = locTime2(i);
    end
end

 % Find sound times when animal is active
        for time=1:length(Sound_Time) 
            sound = Sound_Time(time);
            window = Sound_Time(time)+2; % when is the sound?

            [c closest_frame_sound] = min(loco_data(:,1)-sound);
            [c closest_frame_window] = min(loco_data(:,1)-window);

            length_sound_window(time) = closest_frame_window-closest_frame_sound;
            
             if time>1
                difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
                closest_frame_window = closest_frame_window-difference_length_window(time);
             end
             
             closest_frame_sound_all(time) = closest_frame_sound;
             closest_frame_window_all(time) = closest_frame_window;
            
             
           %  isLocoSound(time) = sum(loco_data(closest_frame_sound:closest_frame_window,3))>0;
             isLocoSound(time) = sum(active_time(closest_frame_sound:closest_frame_window))>0;
        end