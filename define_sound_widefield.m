function [parameters,setup,test]=define_sound_widefield(parameters,setup)

for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    date=setup.expt_date{(a)};
    Imaging_Block=setup.BOT_maps(a,:)
    
    
    for i=1:length(setup.BOT_maps)
        BOT_number = num2str(setup.BOT_maps(i));
        voltage_number = sprintf('%03d',setup.voltage_recording);
        folder = sprintf([setup.path_name '/' mouseID '/' date '/VoltageRecording_' mouseID '_widefield_RF_630-' voltage_number]); %direct to specific Tosca folder within a
        cd(folder);
        filename = ['VoltageRecording_' mouseID '_widefield_RF_630-' voltage_number '_Cycle00001_VoltageRecording_001.csv'];
        M = csvread(filename, 1,0);
        start = M(find( M(:,2) > 3, 1 ) )./1000;% find the first time that tosca sends a signal to VoltageRecording (in seconds)
        parameters.adjusted_times(1,:)=(parameters.New_sound_times-parameters.start_time)+start;
        test.M_col2=find(M(:,2) > 3);
        test.M_col4=find(M(:,4) > 3);
%         test.M_col6=find(M(:,6) > 1);
    end
    
    
    for i=1:length(setup.BOT_maps)
        % Let's find the time stamp for each frame
         BOT_number = num2str(setup.BOT_maps(i));
        folder = sprintf([setup.path_name '/' mouseID '/' date '/' BOT_number]); %direct to BOT data folder
        cd(folder);
       
        filename = ['BOT_' mouseID '_widefield_RF_630-00' BOT_number '_Cycle00001-botData.csv'];
        frame_data = csvread(filename, 1,0);
        timestamp = frame_data(:,1);
        test.timestamp=timestamp;
        time_adjust = timestamp(1,1);
        parameters.timestamp = timestamp;%-time_adjust; %there is a small delay between voltage recording and BOT timestamps. This corrects for it.
        
    end
    
end
end




