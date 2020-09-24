
function [loco_data,active_time] = locomotor_activity(loco_data,filename,setup)
% read_loco= [mouseID '-Session' Tosca_Session '-Run' Tosca_Run_number '.loco.txt'];
%loco_data = dlmread(read_loco);%locomotor data
r=loco_data.t(:);% I am only looking at column 1
r=r(:,1)-r(1,1);
B=r;


% loco_timedif=diff(r);%difference between time stamps
% B = padarray(loco_timedif, [1,0], 'pre');%add a 0 to first time stamp (it is lost during the diff function...?)
%%
% load the voltage recording data
%Imaging_Num =  sprintf( '%03d', Imaging_Block);
%       folder = sprintf([setup.path_name setup.username '/' mouseID '/' date '/' Imaging_Block_String]); %direct to specific Tosca folder within a 
%         addpath (folder); %navigate to this folder

% filename = ['BOT_' mouseID '_ReceptiveField' Imaging_Num '_Cycle00001_VoltageRecording_001.csv'];
M = csvread(filename, 1,0);
time_zero=find(M(:,2)>3,1,'first')-1;%time on VoltageRecording that corresponds with t=0 on locomotor data
time_zero=time_zero./1000;
B=r;
% B=B*1000; %change sec to msec
% B(B==0)=time_zero;%Put loco timescale on Voltage recording timescale


new_times = B;
for i=1:length(new_times);
    new_times(i) = new_times(i) + time_zero;
end
% scaled_times = (new_times./1000);
loco_data.t(:) = new_times(:,1); %put the new timestamps into loco_data


%I used this step to take the absolute value of the data, but whether the
%mouse is running forward or backward could be interesting in the future
loco_data.speed(:)=(abs(loco_data.speed(:))); %take absolute value of data


% figure; hold on
% title('Locomotor activity')
% ylabel('Activity')
% xlabel('Seconds')
% plot(loco_data(:,1), loco_data(:,3));
%% when is the animal moving?
active_time=zeros(i,1);
for i=1:length(loco_data.speed)
    if loco_data.speed(i)>setup.constant.locoThresh
        active_time(i) = loco_data.t(i);
    end
end
%% save

% Imaging_Num2 = num2str(Imaging_Block);

% filename=[mouseID 'locomotor' Imaging_Num '.mat']
% folder = ['C:\2P analysis\Suite2P analyzed data\' mouseID '\' date '\' Imaging_Num2];
% addpath(folder);
% save([folder '\' filename]); 
% c=smooth(active_time, 10);
% hold on
% ax1 = gca; % current axes
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% x2 = loco_data(:,1);
% y2 = loco_data(:,3);
% plot(x2, y2, 'Parent',ax2,'Color','k'),

end
