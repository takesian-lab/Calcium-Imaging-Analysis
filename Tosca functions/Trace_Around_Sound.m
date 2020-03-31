function  [trace_around_sound] = Trace_Around_Sound(All_Images_df_over_f,Sound_Time,timestamp);


size_dim1 = size(All_Images_df_over_f,1);
size_dim2 = size(All_Images_df_over_f,1);
trace_around_sound=zeros(size_dim1,size_dim2,51,120);

for time=1:length(Sound_Time)
     
            sound = Sound_Time(time);
            before = Sound_Time(time)-1;
            after = Sound_Time(time)+4;
            window = Sound_Time(time)+2;
            
            [c closest_frame_before] = min(abs(timestamp(:)-before));
            [c closest_frame_sound] = min(abs(timestamp(:)-sound));
            [c closest_frame_after] = min(abs(timestamp(:)-after));
            [c closest_frame_window] = min(abs(timestamp(:)-window));
           
           length_sound_trial(time) = closest_frame_after-closest_frame_before;
        
            if time>1
                difference_length_sound_trial(time) = length_sound_trial(time)-length_sound_trial(1);
                closest_frame_after = closest_frame_after-difference_length_sound_trial(time);
            end
            time
           % smoothed_trace_around_sound_trial(:,:,:,time) =  smooth(All_Images_df_over_f(:,:,closest_frame_before:closest_frame_after),3);
           for x=1:size_dim1
               for y=1:size_dim2
                trace_around_sound(x,y,:,time) =  All_Images_df_over_f(x,y,closest_frame_before:closest_frame_after);
               end
           end
end