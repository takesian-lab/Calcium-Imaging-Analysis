function [traces]=sound_response_intrinsic(parameters,setup, All_Images_df_over_f);

loops=parameters.loops;
timestamp=parameters.timestamp;
for ll=1:loops
    loop_num=num2str(ll)
    Tile = All_Images_df_over_f.(['Tile' loop_num]);
    x=size(Tile,1);
    y=size(Tile,2);
    z=size(Tile,3);
    for f=1:length(parameters.frequencies);
        fnum=num2str(parameters.frequencies(f));
        for lv=1:length(parameters.levels);
            lvnum=num2str(parameters.levels(lv));
            idx=parameters.stimIDX{f,lv};
            Sound_Time=parameters.adjusted_times(idx); 
            for sound = 1:length(Sound_Time);
               
                sound_time = Sound_Time(sound);
                before_time = Sound_Time(sound)-setup.baseline_window;
                after_time = Sound_Time(sound)+setup.end_window;
                window_time = Sound_Time(sound)+setup.sound_window_end;
                
                [c closest_frame_before] = min(abs(timestamp(:)-before_time));
                [c closest_frame_sound] = min(abs(timestamp(:)-sound_time));
                [c closest_frame_after] = min(abs(timestamp(:)-after_time));
                [c closest_frame_window] = min(abs(timestamp(:)-window_time));
                
                
                length_sound_trial(sound) = closest_frame_after-closest_frame_before;
                %    length_baseline_trial(sound) = closest_frame_sound-closest_frame_before;
                
                %             if matching_value(matching_value)>1 %to make all the sound traces the same
                %                 difference_length_sound_trial(matching_value) = length_sound_trial(matching_value)-length_sound_trial(1);
                %                 closest_frame_after = closest_frame_after-difference_length_sound_trial(matching_value);
                %             end
                %
                %             if sound>1 %to make all the baseline traces the same
                %                 difference_length_baseline(sound) = length_baseline_trial(sound)-length_baseline_trial(1);
                %                 closest_frame_sound = closest_frame_sound-difference_length_baseline(sound);
                %             end
                
                trace_around_sound(:,:,:,sound) =  Tile(:,:,closest_frame_before:closest_frame_after);
%                 stimBaseline(:,:,:,sound)=Tile(:,:,closest_frame_before:closest_frame_sound);
%                 wind(:,:,:,sound)=Tile(:,:,closest_frame_sound:closest_frame_window);
                
            end
            traces.(['Tile' loop_num]){f,lv}=trace_around_sound;
%             stimbase.(['Tile' loop_num]){f,lv}=stimBaseline;
%             window.(['Tile' loop_num]){f,lv}=wind;
%             
            clear trace_around_sound 
%             clear stimBaseline wind 
        end
    end
    
end

end