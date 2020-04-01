function [block] = define_sound_singleblock(block,setup)

display('Pulling out Bruker data...');

%Needed from setup:
%block_path
%mousename
%Imaging_block
%stim_protocol

%%

Imaging_Block=setup.imaging_set;             
Imaging_Block_String = num2str(Imaging_Block);
Imaging_Num =  sprintf( '%03d', Imaging_Block);

New_sound_times=block.New_sound_times;
start_time=block.start_time;

%find the file path - this should be agnostic to the naming system
%of the user as long as the files are in format:
% path/username/mouseID/date/block/file

folder = sprintf([setup.path_name setup.username '/' mouseID '/' date]);
cd(setup.block_path)
filedir = dir('*/'); %LOOKS DOWN ONE FOLDER
filenames = {filedir(:).name};

VR_files = filenames(contains(filenames,'VoltageRecording'));
BOT_files = filenames(contains(filenames,'botData'));
BOT_num = Imaging_Block_String;

if setup.stim_protocol == 4 %widefield - look for files containing voltage_number
    voltage_number = sprintf('%03d',setup.voltage_recording{a,1}(i)); %For widefield
    VR_num = voltage_number;
else %2p - look for files containing Imaging_Block_String
    VR_num = Imaging_Block_String;
end

%Find the BOT or VoltageRecording files by number
for f = 1:2 %Bot, VR
    if f == 1
        temp = BOT_files;
        NUM = BOT_num;
        nLettersToRemove = numel('_Cycle00001-botData.csv');
    elseif f == 2
        temp = VR_files;
        NUM = VR_num;
        nLettersToRemove = numel('_Cycle00001_VoltageRecording_001.csv');
    end

    matchingFiles = zeros(size(temp));
    for ff = 1:length(temp)
        %First delete the end of the file name
        shortFileName = temp{ff}(1:end-nLettersToRemove);
        %Then look for NUM. This method is agnostic to the number of
        %digits in NUM (e.g. '9' vs. '10')
        nDigits = numel(NUM) - 1;
        if numel(shortFileName) > nDigits
            matchingFiles(ff) = strcmp(shortFileName(end-nDigits:end),NUM);
        end
    end

    possibleFiles = [temp(matchingFiles == 1)];
    if f == 1
        BOT_filename = possibleFiles{contains(possibleFiles,'.csv')};
        BOT_folder = filedir(find(strcmp(filenames, BOT_filename))).folder;
    elseif f == 2
        VR_filename = possibleFiles{contains(possibleFiles,'.csv')};
        VR_folder = filedir(find(strcmp(filenames, BOT_filename))).folder;
    end
end

%start with the Voltage recording files - this is what we use to
%synchronize BOT data and Tosca data

%now lets get the relevant data from Voltage recording. This is when Tosca
%sends a 5V burst to the Bruker system. The first time this occurs is t=0
%on Tosca

cd(VR_folder)
display(['Loading ' VR_filename])
M = csvread(VR_filename, 1,0);
start = M(find( M(:,2) > 3, 1 ) )./1000;% find the first time that tosca sends a signal to VoltageRecording (in seconds)
Sound_Time(1,:)=(New_sound_times-start_time)+start;
loco_times=data.([mouseID]).(['ImagingBlock' Imaging_Num]).loco_times;
locTime2 = start+loco_times;

% Let's find the time stamp for each frame. This requires to pull out
% the BOT data and correct for the short (>5ms delay) between Voltage
% recording and BOT

cd(BOT_folder) %The only time this might be different for VR_folder is for widefield
display(['Loading ' BOT_filename])
frame_data = csvread(BOT_filename, 1,0);
timestamp = frame_data(:,1)-frame_data(1,1);% this is where we make that small correction
data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp=timestamp;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time=Sound_Time;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).locTime=locTime2;

%Record filepaths for later
data.([mouseID]).(['ImagingBlock' Imaging_Num]).VR_folder = VR_folder;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).VR_filename = VR_filename;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).BOT_folder = BOT_folder;
data.([mouseID]).(['ImagingBlock' Imaging_Num]).BOT_filename = BOT_filename;


end



