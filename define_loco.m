function [data] = define_loco(setup,data)
isLoco = []; 
display('Finding active trials')
for a=1:length(setup.mousename)
    mouseID=setup.mousename{a};
    Tosca_Sessions =setup.Session{a};
    Imaging_Block=setup.Imaging_sets{a,1};
    
for i=1:length(Imaging_Block)
        date= char(setup.expt_date{a,1}(i));
        Imaging_Block_String = num2str(Imaging_Block(i));
        Imaging_Num =  sprintf( '%03d', Imaging_Block(i));
        Tosca_Session = num2str(Tosca_Sessions(i));
        Tosca_Run_number=num2str(setup.Tosca_Runs{a,1}(i));
        Sound_Time=data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;

     
       
       %read tosca file
       Tosca_folder_name = ['Tosca_' mouseID]; %name of the Tosca folder 
       folder = sprintf([setup.path_name setup.username '/' mouseID '/' Tosca_folder_name '/Session ' Tosca_Session]); %direct to specific Tosca folder within a 
       read_loco= [mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt'];
       cd(folder)
       loco_data = dlmread(read_loco);%locomotor data 
       
       %Read BOT file
       folder = data.([mouseID]).(['ImagingBlock' Imaging_Num]).BOT_folder;
       filename = data.([mouseID]).(['ImagingBlock' Imaging_Num]).BOT_filename;
       cd(folder)
       M = csvread(filename, 1,0);
        
       %Prepare to read VoltageRecording file in locomotor_activity
       folder = data.([mouseID]).(['ImagingBlock' Imaging_Num]).VR_folder;
       filename = data.([mouseID]).(['ImagingBlock' Imaging_Num]).VR_filename;
       cd(folder)
       [loco_data,active_time] = locomotor_activity(loco_data,M,read_loco,setup,mouseID,date,Imaging_Block_String,Imaging_Block,filename);
       data.([mouseID]).(['ImagingBlock' Imaging_Num]).locomotion_data = loco_data;
       data.([mouseID]).(['ImagingBlock' Imaging_Num]).active_time = active_time;


     
%        redcell = data.([mouseID]).(['ImagingBlock' Imaging_Num]).redcell;
%        nonredcell = data.([mouseID]).(['ImagingBlock' Imaging_Num]).nonredcell;
%        F7 = data.([mouseID]).(['ImagingBlock' Imaging_Num]).full_trace;
       
%        mean_F7_green = squeeze(mean(F7(nonredcell,:),1));
%        SEM_F7_green = std(F7(nonredcell,:),1)./sqrt(size(F7(nonredcell,:),1));
       
%        mean_F7_red = squeeze(mean(F7(redcell,:),1));
%        SEM_F7_red = std(F7(redcell,:),1)./sqrt(size(F7(redcell,:),1));
%        
%        timestamp = data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
%        Sound_Time = data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;
       
%        figure; 
%        plot(timestamp,smooth(mean_F7_green,10),'-b','LineWidth',3); hold on;
%        plot(timestamp,smooth(mean_F7_red,10),'r','LineWidth',3);hold on;
%        shadedErrorBar(timestamp,smooth((mean_F7_green),10),smooth((SEM_F7_green),10),'lineprops','-b','transparent',1); hold on;
%        shadedErrorBar(timestamp,smooth((mean_F7_red),10),smooth((SEM_F7_red),10),'lineprops','-r','transparent',1); hold on;
%     %   plot(loco_data(:,1),loco_data(:,3)*10+mean(mean_F7_green)/3,'-k'); hold on; % loco data is 3 columns (timestamp, activity, active yes (1) or no (0))
%        plot(loco_data(:,1),loco_data(:,3)*4,'-k'); hold on;
%        %   h=vline(Sound_Time, 'b');
%          hold on; 
         
         
        % Find sound times when animal is active
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
        
        %Plot locomotor activity
        figure;
        
        subplot(2,1,1); hold on
        title('Locomotor activity')
        ylabel('Activity')
        plot(loco_data(:,1), loco_data(:,3));
        
        subplot(2,1,2); hold on
        ylabel('Considered active')
        xlabel('Seconds')
        plot(loco_data(:,1), active_time > 0); hold on;
        
        
        data.([mouseID]).(['ImagingBlock' Imaging_Num]).isLoco = Loco_1;
           
        %now pull out only isLocoSound trials for red and green, average and
        %concatenate across mice
%         trace_around_sound_green = data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_green;%trace(dF/F) around sound for each non-VIP (cellxsoundxtime matrix)
%         isLocoSound_green = mean(trace_around_sound_green(:,isLocoSound,:),2);
      %  size_isLocogreen=size(trace_around_sound_green(:,isLocoSound,:))
        
%         
%         noLocoSound_green = mean(trace_around_sound_green(:,isLocoSound==0,:),2);
      %  size_noLocogreen=size(trace_around_sound_green(:,isLocoSound==0,:))
%         all_cells_noLocoSound_green = [all_cells_noLocoSound_green; noLocoSound_green]; 
%         
%         trace_around_sound_red = data.([mouseID]).(['ImagingBlock' Imaging_Num]).trace_around_sound_red;%trace(dF/F) around sound for each non-VIP (cellxsoundxtime matrix)
%         isLocoSound_red = mean(trace_around_sound_red(:,isLocoSound,:),2);
%         all_cells_isLocoSound_red = [all_cells_isLocoSound_red; isLocoSound_red]; 
%         
%         noLocoSound_red = mean(trace_around_sound_red(:,isLocoSound==0,:),2);
%         all_cells_noLocoSound_red = [all_cells_noLocoSound_red; noLocoSound_red]; 
    
  isLoco = [isLoco, Loco_1];                                  
end
data.([mouseID]).parameters.loco=isLoco;
end
% data.combined.trace.isLocoSoundgreen_all=all_cells_isLocoSound_green;
% data.combined.trace.isLocoSoundred_all=all_cells_isLocoSound_red;
% 
% data.combined.trace.noLocoSoundgreen_all=all_cells_noLocoSound_green;
% data.combined.trace.noLocoSoundred_all=all_cells_noLocoSound_red;
end