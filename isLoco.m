
function [active_time,isLocoSound,data] = isLoco(setup,data)

for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    Tosca_Session=setup.Session{(a)}
    date=setup.expt_date{(a)};
    Imaging_Block=setup.Imaging_sets(a,:)
    Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder
    Tosca_Run_number = num2str(setup.Tosca_Runs(a));
    folder = sprintf([setup.path_name setup.username '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]);
    cd(folder)
    
    
    isLocoSound=[];
    
    for i=1:length(Imaging_Block(a,:))
        Imaging_Block_String = num2str(Imaging_Block(i))
        Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
        Sound_Time=data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;
        locTime2=data.([mouseID]).(['ImagingBlock' Imaging_Num]).locTime;
        loco_activity=data.([mouseID]).(['ImagingBlock' Imaging_Num]).loco_activity;
        
        % when is the animal moving?
        active_time=zeros(i,1);
%         active_time=zeros(length(loco_activity));
        for i=1:length(loco_activity)
            if loco_activity(i)>0.8
                active_time(i) = locTime2(i);
            end
        end
        
        % Find sound times when animal is active
        for time=1:length(Sound_Time)
            sound = Sound_Time(time);
            window = Sound_Time(time)+2; % when is the sound?
            
            [c closest_frame_sound] = min(abs(locTime2(:,1)-sound));
            [c closest_frame_window] = min(abs(locTime2(:,1)-window));
            
            length_sound_window(time) = closest_frame_window-closest_frame_sound;
            
            if time>1
                difference_length_window(time) = length_sound_window(time)-length_sound_window(1);
                closest_frame_window = closest_frame_window-difference_length_window(time);
            end
            
            closest_frame_sound_all(time) = closest_frame_sound;
            closest_frame_window_all(time) = closest_frame_window;
            
            
            %  isLocoSound(time) = sum(loco_data(closest_frame_sound:closest_frame_window,3))>0;
            loco_1(time) = sum(active_time(closest_frame_sound:closest_frame_window))>0;
        end
        isLocoSound=[isLocoSound,loco_1];
    end
    data.([mouseID]).parameters.loco=isLocoSound;
end
end