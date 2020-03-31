
%% Sort 2-P Data by Locomotor Activity
    %Anne Takesian - 3/4/2019


%Define analysis-specific info here:        
        user='Carolyn';
        mouseID='VxAA060118M2';
        date='2018-11-07';
        Imaging_Block = [1]; 

 for block = 1:length(Imaging_Block)       
    Imaging_Num =  sprintf( '%03d', Imaging_Block(block));
    Imaging_Num2 = num2str(Imaging_Block(block));
    cd('C:\');        
    addpath(genpath('2P analysis'));
    folder = ['C:\2P analysis\Suite2P analyzed data\' mouseID '\' date '\' Imaging_Num2];
    load([mouseID 'noiseburst_analysis' Imaging_Num '.mat'])
    load([mouseID 'locomotor' Imaging_Num '.mat'])
    
    
    timestamp = data.([mouseID]).(['ImagingBlock' Imaging_Num]).timestamp;
    Sound_Time = data.([mouseID]).(['ImagingBlock' Imaging_Num]).Sound_Time;
    
   neuropil=dat.FcellNeu{1,1}; % all the neuropil traces (should be one for each trace)
    cell=dat.Fcell{1,1}; % all the cell fluorescence data
    deconvolution = dat.sp{1,1}; % pulls out all the deconvultion traces
%     statistics = dat.stat{1,1}; % pulls out all the stats
    cell_id=find([dat.stat.iscell]==1); % finds the rows with actual cells (determined manually in GUI)
    
  figure; 
  nonactive_VIP_cells = [];
  active_VIP_cells = [];
  nonactive_nonVIP_cells = []; 
  active_nonVIP_cells = []; 
     for h=1:length(cell_id);
         row_num = cell_id(h);
         cell_trace = cell(h,:);
         mean_gCAMP = mean(cell_trace);
         df_f = (cell_trace-mean_gCAMP)./mean_gCAMP;
         A = smooth(df_f,10);   
         A =   squeeze(A);
         length(timestamp)
         length(A)
         length(cell)
         plot(timestamp(1:length(timestamp)),A +1*h,'LineWidth',2);
         line=vline(Sound_Time);
         hold on; plot(Sound_Time,1,'o');
         %average active or non-active times
         hold on;
         plot(loco_data(:,1), loco_data(:,3)./15,'k'); 
    
         VIP_cell_id = VIP_cell_numbers{block,:};
         r = cell_id(h);
             active_time_index = find(active_time~=0);
            active = active_time(active_time_index);w
            non_active = timestamp;
            non_active(active_time_index) = [];
            non_active_time_index = find(~ismember(timestamp,active_time));
         
         if ismember (r,VIP_cell_id)  
            nonactive_VIP_cells = [nonactive_VIP_cells mean(cell_trace(:,non_active_time_index))];
            active_VIP_cells = [active_VIP_cells mean(cell_trace(:,active_time_index))];
         else
            nonactive_nonVIP_cells = [nonactive_nonVIP_cells mean(cell_trace(:,non_active_time_index))];
            active_nonVIP_cells = [active_nonVIP_cells mean(cell_trace(:,active_time_index))];
         end
     end
     
Amean=mean(cell(cell_id,:),1);
Amean=smooth(Amean);
figure; plot(timestamp(1:length(timestamp)),Amean);
   h=vline(Sound_Time);
         hold on; plot(Sound_Time,1,'o');
      
%average either VIP or non-VIP and plot
VIP_cell_id = VIP_cell_numbers{block,:};
VIP_mean=mean(cell(VIP_cell_id,:),1);

non_VIP_cell_id = non_VIP_cell_numbers{block,:};
non_VIP_mean=mean(cell(non_VIP_cell_id,:),1);

VIP_mean_smooth=smooth(VIP_mean);
non_VIP_mean_smooth=smooth(non_VIP_mean);

figure; plot(timestamp(1:length(timestamp)),VIP_mean_smooth,'r'); hold on;
plot(timestamp(1:length(timestamp)),non_VIP_mean_smooth,'b');
   h=vline(Sound_Time);
   hold on;
plot(loco_data(:,1), loco_data(:,3)*10+mean(VIP_mean_smooth),'g'); 



 end
    




