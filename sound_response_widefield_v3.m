function [traces,parameters]=sound_response_widefield_v3(parameters,data,All_Images_df_over_f);
setup = data.setup;

for i=1:length(setup.Imaging_sets)
    mouseID=setup.mousename{i};
    unique_block_name = setup.unique_block_names{i};
    block = data.([mouseID]).([unique_block_name]);
    
    loops=parameters.loops;
    timestamp=block.timestamp;
    for ll=1:loops
        ll
        loop_num=num2str(ll);
        Tile = All_Images_df_over_f.(['Tile' loop_num]);
        x=size(Tile,1);
        y=size(Tile,2);
        z=size(Tile,3);
        for f=1:length(parameters.frequencies);
            fnum=num2str(parameters.frequencies(f));
            for lv=1:length(parameters.levels);
                lvnum=num2str(parameters.levels(lv));
                
                % plot all trials or only non-loco trials
                if parameters.sort_loco ==0
                    idx=parameters.stimIDX{f,lv};
                    if lv ==1 & f==1 & ll==1
                        display('...creating traces for all trials...')
                    end
                else idx = parameters.loco_0.stimIDX{f,lv};
                    if lv ==1 & f==1 & ll==1
                        display('...creating traces for non-motor trials...')
                    end
                end
            
            
            %adjusted sound times
            if parameters.use_adjusted==1;
                Sound_Time=parameters.adjusted_times(idx);
                if lv ==1 & f==1 & ll==1
                    display('...use adjusted sound times...')
                end
            else Sound_Time=block.Sound_Time(idx);
                if lv==1 & f==1 & ll==1
                    display('...use original sound times...')
                end
            end
            TF = isempty(idx);
            if TF==0;
                for sound = 1:length(Sound_Time);
                    sound_time = Sound_Time(sound);
                    before_time = Sound_Time(sound)-block.setup.constant.baseline_length;
                    after_time = Sound_Time(sound)+block.setup.constant.after_stim;
                    window_time = Sound_Time(sound)+block.setup.constant.response_window;
                    
                    [c closest_frame_before] = min(abs(timestamp(:)-before_time));
                    [c closest_frame_sound] = min(abs(timestamp(:)-sound_time));
                    [c closest_frame_after] = min(abs(timestamp(:)-after_time));
                    [c closest_frame_window] = min(abs(timestamp(:)-window_time));
                    
                    
                    length_sound_trial(sound) = closest_frame_after-closest_frame_before;
                    trace_around_sound(:,:,:,sound) =  Tile(:,:,closest_frame_before:closest_frame_after);
                end
                
            end
            traces.(['Tile' loop_num]){f,lv}=trace_around_sound;
            
        end
       
    end
     clear trace_around_sound
    
end
end
end


