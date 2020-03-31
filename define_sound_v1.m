function [Sound_Time,timestamp,sound_v1,sound_v2,locTime2]=define_sound(sound_v1,sound_v2,stim_protocol,i,path_name,username,mouseID,date,Var1,Var2,New_sound_times,voltage_recording,start_time,Imaging_Block_String,Imaging_Num,loco_times,isLocoSound)



%pull out the two variables (not needed for noiseburst)...change variable
%names?
if stim_protocol>1;
    frequency = struct2cell(Var1 (i,:));
    level = struct2cell(Var2 (i,:));
end

%find the file path - I usually name them differently for different types
%of stim... this could be made to be more efficient

%start with the Voltage recording files - this is what we use to
%synchronize BOT data and Tosca data

voltage_number = sprintf('%03d',voltage_recording);
if stim_protocol==4 %widefield
    folder = sprintf([path_name username '/' mouseID '/' date '/VoltageRecording_' mouseID '_widefield_gcamp-' voltage_number]); %direct to specific Tosca folder within a
    cd(folder); %navigate to this folder
    filename = ['VoltageRecording_' mouseID '_widefield_gcamp-' voltage_number '_Cycle00001_VoltageRecording_001.csv'];
elseif stim_protocol==1%noiseburst
    folder = sprintf([path_name username '/' mouseID '/' date '/' Imaging_Block_String]);
    cd (folder); %navigate to this folder
    filename = ['BOT_' mouseID '_noiseburst-' Imaging_Num '_Cycle00001_VoltageRecording_001.csv'];
    %put in here somehting for other stim_protocols
end





%now lets get the relevant data from Voltage recording. This is when Tosca
%sends a 5V burst to the Bruker system. The first time this occurs is t=0
%on Tosca

% folder = sprintf([path_name username '/' mouseID '/' date '/' Imaging_Block_String]);
% cd(folder)
    M = csvread(filename, 1,0);
    start = M(find( M(:,2) > 3, 1 ) )./1000;% find the first time that tosca sends a signal to VoltageRecording (in seconds)
    Sound_Time(1,:)=(New_sound_times-start_time)+start;
    if stim_protocol>1
        sound_1 = squeeze(struct2cell(Var1 (:,:)));
        sound_2 = squeeze(struct2cell(Var2 (:,:)));
        %                         Sound_Frequency(sound) = frequency(:,:,sound);
        %                         Sound_Level(sound) = level(:,:,sound);
    else sound_1=0;
         sound_2=0;
    end
   locTime2 = start+loco_times;
    
    
    
    % Let's find the time stamp for each frame. This requires to pull out
    % the BOT data and correct for the short (>5ms delay) between Voltage
    % recording and BOT
    folder = sprintf([path_name username '/' mouseID '/' date '/' Imaging_Block_String]); %direct to BOT data folder
    cd(folder);
    if stim_protocol==1
        filename = ['BOT_' mouseID '_noiseburst-' Imaging_Num '_Cycle00001-botData.csv'];
    elseif stim_protocol==4
        filename = ['BOT_' mouseID '_widefield_gcamp-00' Imaging_Num '_Cycle00001-botData.csv'];
    end
    frame_data = csvread(filename, 1,0);
    timestamp = frame_data(:,1)-frame_data(1,1);% this is where we make that small correction   
    
    sound_v1=[sound_v1;sound_1];
    sound_v2=[sound_v2;sound_2];
    
    
    
end



