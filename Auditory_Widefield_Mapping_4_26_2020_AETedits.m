
%% Generating Auditory Cortical Maps
%Anne Takesian (AET) - 2/22/2019
%small updates CGS 4/12/19
%AET updated 6/20/19 to be consistent with Ed's analysis
%AET and CGS updated 4/20/20

%added variable image_f to define imaging freq (SETUP)
%changed filtered images to be 256x256 frames (EXTRACT IMAGES)
%changed timestamp to auto subtract delay to first frame(TIMESTAMP)

% Carolyn updated in January 2020
%   two new functions 'behavior_RF' and 'define_sound'
%   combined steps duing detrend and df/f to help save memory

%This program outputs large sound-evoked activity maps across the cranial window to
%identify A1 and secondary auditory cortices using widefield imaging.

%The program should proceed through several functions (see associated .m
%files) that are described below. In order to run this program, the
%following m files should be in the same folder.
%tosca_trace.m
%tosca_read_trial.m
%tosca_read_run.m
%test_analysisEDITAD.m
%runline.m
%combineVoltageBehaviorTime.m
%A1_mapping.m
%A1_map_tile.m
%locdetrend.m
%Sound responses
%behavior_RF
%defome_spound
%pixel_detrend_widefield

%To run, the data needs to be set up in the same path everytime.
%Define experiment-specific data below.

%% SETUP
% set up your files in folders in D drive using the following
setup.username = 'Carolyn'
setup.mousename={'CnJ'};
setup.expt_date={'2020-03-06'};

setup.Session ={'2'};
setup.Tosca_Runs = [1];
setup.BOT_maps = [1];
setup.voltage_recording = [1];
setup.path_name = 'C:/'
temp_spat_analysis = [0]; %do you want temporal and spatial filtering? adds about 30 hours 
detrend_filter = [300 10];

%define which number the BOT file starts (n) and the number BOT
setup.BOT_start = [1];
parameters.image_f = 10;

setup.ReceptiveField = 1;

setup.run_water=0;
p = gcp('nocreate'); 
if isempty(p)
    %%% If no pool, do not create new one.
    parpool('local',10);
end

%% Behavioral data to generate sound times
%this works right now as long as it is one animal and one BOT file - it may
%break with more.
[parameters] = behavior_RF_widefield(parameters,setup);
[parameters,setup,test]=define_sound_widefield(parameters,setup);

%% load images and crop window
[imageData,parameters] = crop_window(setup,parameters);
 
 %% anne's "old code" to test for window quality and check sound stimulation
 Full_Tile_Mean = mean(mean(imageData.Cropped_Imaging_Data,1),2);
 Full_Tile_Mean_Detrend = locdetrend(Full_Tile_Mean(1,1,:),1,detrend_filter); %
 timestamp = parameters.timestamp; %change this number
      figure;   
      hold on;
       
      fo =squeeze(Full_Tile_Mean)-squeeze(Full_Tile_Mean_Detrend);
      f = squeeze(Full_Tile_Mean);  
      plot(timestamp(1:length(timestamp)), squeeze(Full_Tile_Mean_Detrend),'c');
      plot(timestamp(1:length(timestamp)), smooth(f-1750,3),'b');
      plot(timestamp(1:length(timestamp)), fo,'m');
      
      h=vline(parameters.adjusted_times(:),'k'); hold on
  %    plot(timestamp(1:length(timestamp)), df_over_fo,'r');
      
      hold on; plot(parameters.adjusted_times(:),6500,'o');
      title(sprintf('Mean Response across all Trials from Tile %d', i));
      %xlim([0 1000])
      pause;
  
 %  average trace around sound 

      total_average = 1;
      Sound_Time = parameters.adjusted_times(:);
      
      for y = 1:size(Sound_Time,1)
        before = Sound_Time(y)-1;
        after = Sound_Time(y)+2.5;
        sound = Sound_Time(y);
        [c closest_frame_before] = min(abs(timestamp(:,1)-before));
        [c closest_frame_sound] = min(abs(timestamp(:,1)-sound));
        [c closest_frame_after] = min(abs(timestamp(:,1)-after));
 
        length_sound_trial(y) = closest_frame_after-closest_frame_before;
        
        if y>1
                difference_length_sound_trial(y) = length_sound_trial(y)-length_sound_trial(1);
                closest_frame_after = closest_frame_after-difference_length_sound_trial(y);
        end
        
        baseline_mean = mean(Full_Tile_Mean_Detrend(closest_frame_before:closest_frame_sound),1);
        average = Full_Tile_Mean_Detrend (closest_frame_before:closest_frame_after);
        smooth_All_Images= smooth(average);  
        all_trials(:,y) = smooth_All_Images;
        

      end
      
 
      mean_across_all_trials = mean(all_trials,2);
%       data.([mouseID]).(['Tile' BOT_number]).Average_Sound_Response = mean_across_all_trials;
      pause; 
      

      figure; 
      plot(1:length(average), squeeze(mean_across_all_trials),'k')
      title(sprintf('Mean Response across all Trials from Tile %d', i));
      
      %AT added 4/15/20 to center window around peak response across window
      
      [amp_average_response,time_average_response] = max(mean_across_all_trials);
      estimated_peak = parameters.image_f*1.8; %the peak should be about 0.8sec after stim (1s)
      estimated_time = (time_average_response-estimated_peak)*(1/parameters.image_f); 
      parameters.adjusted_times_est = parameters.adjusted_times+estimated_time;
   
   %  repeat this average trace around sound

      total_average = 1;
      Sound_Time = parameters.adjusted_times_est(:);
      
      for y = 1:size(Sound_Time,1)
        before = Sound_Time(y)-1;
        after = Sound_Time(y)+3.5;
        sound = Sound_Time(y);
        [c closest_frame_before] = min(abs(timestamp(:,1)-before));
        [c closest_frame_sound] = min(abs(timestamp(:,1)-sound));
        [c closest_frame_after] = min(abs(timestamp(:,1)-after));
 
        length_sound_trial(y) = closest_frame_after-closest_frame_before;
        
        if y>1
                difference_length_sound_trial(y) = length_sound_trial(y)-length_sound_trial(1);
                closest_frame_after = closest_frame_after-difference_length_sound_trial(y);
        end
        
        baseline_mean = mean(Full_Tile_Mean_Detrend(closest_frame_before:closest_frame_sound),1);
        average = Full_Tile_Mean_Detrend (closest_frame_before:closest_frame_after);
        smooth_All_Images= smooth(average);  
        all_trials(:,y) = smooth_All_Images;
      end
      
      mean_across_all_trials = mean(all_trials,2);
      figure; 
      plot(1:length(average), squeeze(mean_across_all_trials),'k')
      title(sprintf('Mean Response after adjustment across all Trials from Tile %d', i));
      
      mean_across_all_trials = mean(all_trials,2);
%       data.([mouseID]).(['Tile' BOT_number]).Average_Sound_Response = mean_across_all_trials;
      pause;    
      
      figure;
      index_high_levels = find(parameters.level_list>60);
      mean_across_all_high_trials = mean(all_trials(:,index_high_levels),2);
      plot(1:length(average), mean_across_all_high_trials,'r')
      
      hold on;
      index_mid_levels = find(parameters.level_list>20 & parameters.level_list<70);
      mean_across_all_mid_trials = mean(all_trials(:,index_mid_levels),2);
      plot(1:length(average), mean_across_all_mid_trials,'y')
      
      hold on;
      index_low_levels = find(parameters.level_list<20);
      mean_across_all_low_trials = mean(all_trials(:,index_low_levels),2);
      plot(1:length(average), mean_across_all_low_trials,'g')
       title(sprintf('Mean Response by Level across all Trials from Tile %d', i));
      pause;
      
      figure;
      index_high_freq = find(parameters.frequency_list>23);
      mean_across_all_highf_trials = mean(all_trials(:,index_high_freq),2);
      plot(1:length(average), mean_across_all_highf_trials,'r')
      
      hold on;
      index_mid_freq = find(parameters.frequency_list>10 & parameters.frequency_list<32);
      mean_across_all_midf_trials = mean(all_trials(:,index_mid_freq),2);
      plot(1:length(average), mean_across_all_midf_trials,'y')
      
      hold on;
      index_low_freq = find(parameters.frequency_list<10);
      mean_across_all_lowf_trials = mean(all_trials(:,index_low_freq),2);
      plot(1:length(average), mean_across_all_lowf_trials,'g')
      title(sprintf('Mean Response by Freq across all Trials from Tile %d', i));
      pause;
      
      figure;
      
      mean_across_all_first_half_trials = mean(all_trials(:,1:300),2);
      plot(1:length(average), mean_across_all_first_half_trials,'r')
      
      hold on;
       mean_across_all_first_half_trials = mean(all_trials(:,300:500),2);
      plot(1:length(average), mean_across_all_first_half_trials,'g')
      
      hold on;
      mean_across_all_second_half_trials = mean(all_trials(:,500:868),2);
      plot(1:length(average), mean_across_all_second_half_trials,'y ')
      title(sprintf('Mean Response by Time across all Trials from Tile %d', i));
  
      
parameters.adjusted_times = parameters.adjusted_times_est;   %AT 
      
%% check for memory
[loops] = memorycheck(imageData);

%% DETREND: Grab a coffee - this will take approx. 2 hours
t = cputime;
parameters.loops=loops;
for i=1:length(setup.BOT_maps)
    BOT_number = num2str(setup.BOT_maps(i));
   Cropped = imageData.Cropped_Imaging_Data;
    [All_Images_df_over_f] = Pixel_Detrend_Widefield_v2(Cropped,loops);   
end
clear Cropped
%
e = cputime-t
clear t e i BOT_number

%% view df_over_f
for i=1:length(setup.BOT_maps)
    %mean of all pixels across time
    BOT_number = num2str(setup.BOT_maps(i));
    FullTile_df= squeeze(mean(mean(All_Images_df_over_f.Tile4,1),2));
    timestamp=parameters.timestamp;
    Sound_Time=parameters.adjusted_times; %changed AT on 4/11/20 to adjust for timing issue
    figure;
    plot(timestamp,FullTile_df,'b');
    hold on; vline(Sound_Time);
    title(sprintf('Mean Response across all Trials from Tile %d', i));
end
clear BOT_number FullTile_df i Sound_Time timestamp


%% create frequency and level indicies and find responses to sound across stim
[parameters] = indexStimuli(parameters,setup);
[traces]=sound_response_widefield_v2(parameters,All_Images_df_over_f);

%% view the traces with sounds

%raw trace
raw = imageData.Cropped_Imaging_Data;
raw = squeeze(mean(mean(raw,1),2));
x=parameters.timestamp;
Sound_Time = parameters.adjusted_times;
figure; plot(x,raw); hold on; vline(Sound_Time);

%% pull out baseline and window
lenghth_trial=size(traces.Tile1{1,1},3);
baseline=1:(0.5*parameters.image_f);
window=(length(baseline)+1):(parameters.image_f*2);

for ll=1:loops
    loop_num=num2str(ll)
     for f=1:length(parameters.frequencies);
        fnum=num2str(parameters.frequencies(f));
        for lv=1:length(parameters.levels);
            lvnum=num2str(parameters.levels(lv));
            idx=parameters.stimIDX{f,lv};
            base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
%             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
        end
     end
end

% get mean and STD baseline



%% mean across stimuli
loops=parameters.loops;
count=1;
for ll=1:loops
    loop_num=num2str(ll);
    for f=1:length(parameters.frequencies);
        for lv=1:length(parameters.levels);
            idx=parameters.stimIDX{f,lv};
            m=traces.(['Tile' loop_num]){f,lv};
            %mean response per stim
            mean_stim(:,:,:)=mean(m,4);
            avgTrace.(['Tile' loop_num]){f,lv}= mean_stim;
            clear mean_stim
            %mean baseline 
%             b=base.(['Tile' loop_num]){f,lv};
%             mean_base=squeeze(mean(mean(mean(b,4),2),1));
%             tempBase(count,:)=mean_base;
            
            %mean window
%             d=stimbase.(['Tile' loop_num]){f,lv};
%              mean_win=squeeze(mean(mean(mean(d,4),2),1));
%              tempWin(count,:)=mean_win;
%              count=count+1;
        end
    end
    
end
% 
% parameters.avgBaseline=mean(tempBase,1);
% parameters.stdBaseline=std(tempBase,1);
clear tempBase idx f b count ll loop_num lv m mean_base

%% put tiles back together
 rejoin_tiles=[];
  %loop through
    
 
      for f=1:length(parameters.frequencies)
          numF=num2str(round(parameters.frequencies(f)))
          for lv=1:length(parameters.levels);
              numLV=num2str(parameters.levels(lv))
             
              rejoin_tiles=[];
               for ll=1:loops
      loop_num=num2str(ll)
              
              mm=avgTrace.(['Tile' loop_num]){f,lv};
              rejoin_tiles=double(cat(1, rejoin_tiles, mm));
          end
          stimAverages.(['kHz' numF ]){lv}=rejoin_tiles;
          clear rejoin_tiles
      end
  end
  clear mm avgTrace

%% what do the stim averages look like?
for f=1:length(parameters.frequencies)
          numF=num2str(round(parameters.frequencies(f)))
          subplot(2,4,f)
         
          for lv=1:length(parameters.levels);
              numLV=num2str(parameters.levels(lv))
              a1 = stimAverages.(['kHz' numF ]){lv};
               a2 = squeeze(mean(mean(a1,1),2));
               
               plot(smooth(a2,5)); hold on
               
               
              
          end
end



%% convert to tif? and then store as individual file. 
    %folder = 'D:\2P analysis\2P local data\Wisam\YD111219F3\2020-03-06\AnalyzedTiffs\';
    
    folder = 'C:\Anne\';
    cd(folder)
    tic
    for f=1:length(parameters.frequencies);
        toc
        numF=num2str(round(parameters.frequencies(f)));
        for lv=1:length(parameters.levels);
            numLV=num2str(parameters.levels(lv));
            idx=parameters.stimIDX{f,lv};
            outputFileName= (['avgStim' numF 'khz' numLV 'db.tif']);
            for k = 1:size(stimAverages.(['kHz' numF ]){lv},3)
                imwrite(stimAverages.(['kHz' numF ]){lv}(:, :, k), outputFileName, 'WriteMode', 'append')
            end
        end
    end

%% temporal and spatial denoise- takes approx. 30 hours

if temp_spat_analysis == 1
t=cputime
    %folder = 'D:\2P analysis\2P local data\Wisam\YD111219F3\2020-03-06\AnalyzedTiffs\';
    folder = 'C:\Anne\';
    cd(folder)
    d = dir([folder '/*.tif']);%extract tiffs
    
   parfor f=1:length(parameters.frequencies)
        numF=num2str(round(parameters.frequencies(f)))
        for lv=1:length(parameters.levels)
            numLV=num2str(parameters.levels(lv))
            idx=parameters.stimIDX{f,lv};
            fname = (['avgStim' numF 'khz' numLV 'db.tif']);
            info = imfinfo(fname);
            num_images = numel(info);
            for k = 1:num_images
                AA = imread(fname, k);
                stack(:,:,k)=AA(:,:);
                stack=double(stack);
            end
            outputFileName= (['TempDenoise' numF 'khz' numLV 'db.tif']);
            disp(['Temporally deconvolving for, ' num2str(parameters.frequencies(f))  ' KHz, ' num2str(parameters.levels(lv)) ' dB'])
           
            
            %             tempD=zeros(size(mm,1),length(y),size(mm,3));
%             clear tempD
            tempD=zeros(size(stack,1),(size(stack,2)),size(stack,3));
            for x=1:size(stack,1);
             
                
                for y=1:size(stack,2);
                    %temporal filter using Paninski code
                    A =(stack(x,y,:));
                    [c,b,c1,g,sn,sp]=constrained_foopsi_new(A); %,b,c1,g,sn,options);
                    tempD(x,y,:) = c;
                end
            end
            for kk = 1:size(tempD,3);
                imwrite(tempD(:, :, kk), outputFileName, 'WriteMode', 'append');
            end
            %
        end
        
    end
e=t-cputime
clear tempD x y stack t sp sn outsputFileName lumLB numF...
    num_images lv ll loops loop_num kk k infor info idx g...
    e d

% spatial deconvolve

t=cputime;
PSF=fspecial('gaussian',[3 3],0.5);
%folder = 'D:\2P analysis\2P local data\Wisam\YD111219F3\2020-03-06\AnalyzedTiffs\';
folder = 'C:\Anne\';
cd(folder)
d = dir([folder '/*.tif']);%extract tiffs
parfor f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)));
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
        idx=parameters.stimIDX{f,lv};
        fname= (['TempDenoise' numF 'khz' numLV 'db.tif'])
        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:num_images
            AA = imread(fname, k);
            stack(:,:,k)=AA(:,:);
            stack=double(stack);
        end
        outputFileName= (['SpatDenoise' numF 'khz' numLV 'db.tif']);
        disp(['Spatially deconvolving for, ' num2str(parameters.frequencies(f))  ' KHz, ' num2str(parameters.levels(lv)) ' dB'])
       
        spatD=zeros(size(stack,1),(size(stack,2)),size(stack,3));
       
                for time=1:size(stack,3)
                    spatD(:,:,time) = deconvlucy(stack(:,:,time), PSF); %spatial filter for each pixel
                end
        
        for kk = 1:size(spatD,3);
            imwrite(spatD(:, :, kk), outputFileName, 'WriteMode', 'append');
        end
    end
end
e=t-cputime
end
%% plot temporal noise responses
%folder = 'D:\2P analysis\2P local data\Wisam\YD111219F3\2020-03-06\AnalyzedTiffs\';
folder = 'C:\Anne';
    cd(folder)
    d = dir([folder '/*.tif']);%extract tiffs

for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)));
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));

        fname1= (['TempDenoise' numF 'khz' numLV 'db.tif'])
        
            info1 = imfinfo(fname1);
            num_images1 = numel(info1);
            for k1 = 1:num_images1
                AA1 = imread(fname1, k1);
                stack1(:,:,k1)=AA1(:,:);
                tempD_stack=double(stack1);
            end
       
        %plot avg stim    
        mean_dff0=squeeze(mean(mean(avgStim_stack,1),2)); 
        x=1:length(mean_dff0);
        figure; plot (x, smooth(mean_dff0,10), '-r'); hold on;
            
            
        if temp_spat_analysis == 1 
            fname2 = (['avgStim' numF 'khz' numLV 'db.tif'])
            fname3 = (['SpatDenoise' numF 'khz' numLV 'db.tif'])
        
            info2 = imfinfo(fname2);
            num_images2 = numel(info2);
            for k2 = 1:num_images2
                AA2 = imread(fname2, k2);
                stack2(:,:,k2)=AA2(:,:);
                avgStim_stack=double(stack2);
            end
            
            info3 = imfinfo(fname3);
            num_images3 = numel(info3);
            for k3 = 1:num_images3
                AA3 = imread(fname3, k3);
                stack3(:,:,k3)=AA3(:,:);
                SpatDenoise_stack=double(stack3);
            end
            mean_temp=squeeze(mean(mean(tempD_stack,1),2));
            mean_spat=squeeze(mean(mean(SpatDenoise_stack,1),2));
            plot(x, smooth(mean_temp,10), 'b');hold on
            plot(x, smooth(mean_spat,10), 'g');
            
            %plot normalized
            figure; plot (x, smooth(mean_dff0,10)./max(smooth(mean_dff0,10)), '-r'); hold on;
            plot(x, smooth(mean_temp,10)./max(smooth(mean_temp,10)), 'b');hold on
            plot(x, smooth(mean_spat,10)./max(smooth(mean_spat,10)), 'g');
        end
    end
end

%% find cumulative baseline and the response window
baseline=1:(0.5*parameters.image_f);
window=(length(baseline)):(parameters.image_f*3);
%folder = 'D:\2P analysis\2P local data\Wisam\YD111219F3\2020-03-06\AnalyzedTiffs\';
folder = 'C:\Anne';
cd(folder)
d = dir([folder '/*.tif']);%extract tiffs
count = 1;
image = imageData.Cropped_Imaging_Data;
total_stim = length(parameters.levels)*length(parameters.frequencies);
accumBase = zeros(size(image,1),size(image,2),total_stim*length(baseline));
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
       % fname = (['SpatDenoise' numF 'khz' numLV 'db.tif']);
       % fname = (['TempDenoise' numF 'khz' numLV 'db.tif']);
       fname = (['avgStim' numF 'khz' numLV 'db.tif']);
        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:num_images
            AA = imread(fname, k);
            stack(:,:,k)=AA(:,:);
            stack=double(stack);
        end
        ResponseWindow{f,lv}=stack(:,:,window);
        DFF0_mean{f,lv} = stack (:,:,:);
        baseLoc = stack(:,:,baseline);
        accumBase = cat(3, accumBase, baseLoc); % why is this such a large number?
        
    
        %         base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
        %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
    clear stack AA
    end
end

stdBaseline=std(accumBase,0,3);
meanAccumBaseline = mean(accumBase,3);

%plot response window
y=DFF0_mean{8,8};
%y=y(150:180,170:200,:);
y= squeeze(mean(mean(y,2),1));
x=1:length(y);
figure;
plot(x,y)

%what does a whole response window look like
 g = mean(ResponseWindow{8,8},3);
%         CLIM = [0 350];
   imagesc(g);
        

%% plot all frequencies and amplitudes
     figure;
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        y = mean(ResponseWindow{f,lv},3); 
        n = ((f-1)*length(parameters.frequencies))+lv;
        subplot(length(parameters.levels),length(parameters.frequencies),n);
        imagesc(y);
    %    caxis([250 300]);
        title(sprintf('Freq %d', round(parameters.frequencies(f),2)));
        axis image;
        set(gca,'XTick',[], 'YTick', [])
    end
end


%% plot frequencies averaged across all levels
     figure;
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
        for k = 1:length(parameters.levels);
            AA = mean(ResponseWindow{f,lv},3);
            stack(:,:,k)=AA(:,:);
            stack=double(stack);
        end
        
        y = mean(stack,3); 
        subplot(1,length(parameters.frequencies),f);
        imagesc(y);
        title(sprintf('Freq %d', f));
        axis image;
        set(gca,'XTick',[], 'YTick', [])
end





                %% 
%% zscore data
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
    figure;
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
        idx=parameters.stimIDX{f,lv};
        mainResponse = ResponseWindow{f,lv};%temporally/spatially filtered response (2s post  sound)
        Df_f0 = DFF0_mean{f,lv};%temporally/spatially filtered mean response (entire trace i.e. baseline and 3s post sound)
        %z score of every frame of response window
        zscR = (bsxfun(@minus, mainResponse, (repmat(meanAccumBaseline,...
            [1 1 length(window)]))))./(repmat(stdBaseline,...
            [1 1 length(window)]));
        %zscore of every frame for entire trace
        zscResponse{f,lv}=zscR;
        zscS = (bsxfun(@minus,Df_f0, (repmat(meanAccumBaseline,...
            [1 1 size(Df_f0,3)]))))./(repmat(stdBaseline,...
            [1 1 size(Df_f0,3)]));
        zscSignal{f,lv} =  zscS;
        
        
        [ d1 d2 frames ] =size(zscR);%size of response window data
        
        
        maxResponse{f,lv} = max(mainResponse,[],3);%max of response window (not zscored)
        meanR{f,lv} = mean(mainResponse);%average of sound response window
        meanZresponse{f,lv} = mean(zscR(:,:,1:15),3);
        plot(squeeze(mean(mean(zscR(:,:,:),1),2)));
        title(sprintf('Freq %d', round(parameters.frequencies(f),2)));
        
        hold on;

    end
end

%% average z-scores across all levels
     figure;
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
        for k = 1:length(parameters.levels);
            AA = meanZresponse{f,lv};
            stack(:,:,k)=AA(:,:);
            stack=double(stack);
        end
        
        y = mean(stack,3); 
        subplot(1,length(parameters.frequencies),f);
        imagesc(y);
        title(sprintf('Freq %d', f));
        axis image;
        set(gca,'XTick',[], 'YTick', [])
     %   caxis([2 15]);
        colormap jet
end

% %% average the DFF0 for each f/lv and see what it looks like
% figure;
% for f=1:length(parameters.frequencies);
%     numF=num2str(round(parameters.frequencies(f)))
%     for lv=1:length(parameters.levels);
%         numLV=num2str(parameters.levels(lv));
%         mainResponse = ResponseWindow{f,lv};%temporally/spatially filtered response (2s post  sound)
%         Df_f0 = DFF0_mean{f,lv};%temporally/spatially filtered
%         g = mean(Df_f0,3);
%       %  CLIM = [0 350];
%      %   imagesc(g,CLIM);
%         imagesc(g);
%         
%         subplotSpot=f+(length(parameters.frequencies))*(length(parameters.levels)-lv);
%         subplot(length(parameters.levels),length(parameters.frequencies),subplotSpot);
%         
%         title({sprintf('%s dB',numLV);sprintf('%s kHz',numF)});
%         
%     end
% end
% 

% figure;
% 
% for V1=1:length(data.([mouseID]).parameters.Var1List)%loop through stim#1
%     for V2=1:length(data.([mouseID]).parameters.Var2List)%loop through stim#2
%         stim_list=data.([mouseID]).parameters.stimIDX{V1,V2};
%         stimIDXpos=(data.([mouseID]).parameters.isRespPosStim(V1,V2,:));%index of responsive cells by stim
%         
%         meanPosResp=squeeze(mean(mean(data.([mouseID]).traces_G(stimIDXpos,stim_list,:),2),1));%avg response of positive responsive cells by stim
%         std_resp=squeeze(std(mean(data.([mouseID]).traces_G(stimIDXpos,stim_list,:),2),1));
%         subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2))
%         subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot),
%         shadedErrorBar(x_green,smooth(meanPosResp,10),smooth(std_resp,10),'lineprops','b')
%         stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
%         stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
%         if setup.stim_protocol==2
%             title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
%         end
%          axis([0 length(x_green) 0 70])
%         
%     end
% end
%% peak - new - still in progress
%find peak of response window
% [dim1 dim2 dim3] = size(imageData.Cropped_Imaging_Data);
% peak_factor=0.2;
% peakWidth=peak_factor*parameters.image_f./2;
% 
% 
% for f=1:length(parameters.frequencies);
%     numF=num2str(round(parameters.frequencies(f)))
%     for lv=1:length(parameters.levels);
%         numLV=num2str(parameters.levels(lv));
%         response=ResponseWindow{f,lv};
%                 for x=1:dim1
%                     for y = 1:dim2
%                         [maxVal max_loc]=max(response(x,y,:),[],3);
%                         if max_loc<3
%                             peak = response(x,y,max_loc);
%                         elseif max_loc>28
%                             peak = response(x,y,max_loc);
%                         else
%                         peak_frames=(max_loc-peakWidth):(max_loc+peakWidth);
%                         peak=mean(response(x,y,peak_frames),3);
%                         end
%                 peakPx(x,y,:) = peak;
%                     end
%                 end
%   
%              stim_peak{f,lv}=peakPx;
%     end
% end
% 
% 





%% determine CFs for each pixel
[dim1 dim2 dim3] = size(imageData.Cropped_Imaging_Data);
CF=NaN(dim1,dim2,1);           
response_threshold = 3;


for x= 1:dim1 % go through all x pixels
    disp(['CF mapping for x pixel ', num2str(x)])
    for y = 1:dim2 %go through all y pixels
        threshold_reached = 0; %reset threshold flag
        threshold_level=length(parameters.levels); %reset threshold level
        for lv = 1:length(parameters.levels)-1  % if level is less than max level and threshold has not been found
            for f = 1:length(parameters.frequencies)
                peak_response = meanZresponse{f,lv}(x,y,:);% 
                if peak_response>response_threshold % if threshold has been found - set by user in parameter file above
                    threshold_level = lv;
                    threshold_reached=1;
                    break
                end
            end
            if threshold_reached == 1
                break
            end
        end
       
         if threshold_level < length(parameters.levels)
                level_level_plusone = threshold_level+1;
            for f = 1:length(parameters.frequencies)
                frequency_num = num2str(round(parameters.frequencies(f)));
                signal_threshold(f) = meanZresponse{f,threshold_level}(x,y,:);
                signal_threshold_plus_one(f) = meanZresponse{f,level_level_plusone}(x,y,:);
                avg_threshold_responses(f) = (signal_threshold(f)+signal_threshold_plus_one(f))/2;
            end
            options = fitoptions('gauss1');
            options.Lower = [0 1 0];
            gauss_fit = fit(parameters.frequencies',avg_threshold_responses', 'gauss1',options);
            
            gauss_curve = gauss_fit(parameters.frequencies(1):0.1:parameters.frequencies(length(parameters.frequencies)));
            x_freq = [parameters.frequencies(1):0.1:parameters.frequencies(length(parameters.frequencies))];
            [peak_amplitude,characteristic_frequency] = max(gauss_fit(x_freq));
       %      figure; %plot gaussian fits
        %     plot(parameters.frequencies,avg_threshold_responses); % use
%              hold on; plot(parameters.frequencies,signal_threshold, 'b'); % use
%              hold on; plot(parameters.frequencies,signal_threshold_plus_one, 'c'); % use
%              hold on; plot(x_freq, gauss_curve);
  %           pause;
%       
            
            CF(x,y,:)=characteristic_frequency*0.1+parameters.frequencies(1);
           imageData.CF=CF;
         end
    end
end



 %% RFs - determine CFs for each pixel
% trace_around_sound = data.([mouseID]).(['Tile' BOT_number]).All_Imaging_Around_Sound;
% sz1=size(trace_around_sound,1);
% sz2=size(trace_around_sound,2);
% CF=zeros(sz1,sz2,1);
% frequencies_new = [4 5.7 8 11.3 16 22.6 32];
% for x=1:sz1
%     x
%     for y=1:sz2  %loop through the pixels
%         threshold_reached = 0;
%         for l = 1:length(levels)-1 & threshold_reached<1
%             threshold=120;
%             for f = 1:length(frequencies_new)
%                 level_num = num2str(round(levels(l)));
%                 frequency_num = num2str(round(frequencies_new(f)));
%                 
%                 peak_response=data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num]).(['Frequency' frequency_num]).peakresponse(x,y,:);
%                 if peak_response>0.00005
%                     threshold=levels(l);
%                     threshold_num=l;
%                     threshold_reached=1;
%                 end
%             end
%         end
%         
%         %data.([mouseID]).(['Tile' BOT_number]).(['Frequency' frequency_num]).threshold=threshold;
%         
%         if threshold < levels(length(levels))
%             level_num_threshold = num2str(round(levels(threshold_num)));
%             level_num_plusone = num2str(round(levels(threshold_num+1)));
%             for f = 1:length(frequencies_new)
%                 frequency_num = num2str(round(frequencies_new(f)));
%                 response_threshold = data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num_threshold]).(['Frequency' frequency_num]).peakresponse(x,y,:);
%                 response_threshold_plus_one = data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num_plusone]).(['Frequency' frequency_num]).peakresponse(x,y,:);
%                 avg_threshold_responses(f) = (response_threshold+response_threshold_plus_one)/2;
%             end
%             
%             gauss_fit = fit(frequencies_new',avg_threshold_responses', 'gauss2');
%             %                             gauss_fit = max(frequencies_new',avg_threshold_responses');
%             [peak_amplitude,characteristic_frequency] = max(gauss_fit(frequencies_new(1):0.1:frequencies_new(length(frequencies_new))));
%             %  figure; plot(frequencies(f),levels(l),peak_response(x,y,:);
%             
%             CF(x,y,:)=characteristic_frequency*0.1+frequencies_new(1);
%             data.([mouseID]).(['Tile' BOT_number]).CF=CF;
%         end
%     end
%end
%%
%         for f = 3:4%1:length(frequencies)
%             for l = 4:5%1:length(levels)
%                 level_num = num2str(round(levels(l)));
%                 frequency_num = num2str(round(frequencies(f)));
%                  for x=50:50%1:sz1
%                         for y=50:50%1:sz2
%                              mean_for_zscore = mean(baselines(x,y,:),3); % should be x pixel, y pixel, mean
%                              sd_for_zscore = std(baselines(x,y,:),3);
%                              peak_response=data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num]).(['Frequency' frequency_num]).peakresponse;
%                              zscore(x,y) = (peak_response-mean_for_zscore(x,y,:))./sd_for_zscore(x,y,:);
%                         end
%                  end
%             end
%         end
%
%         data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num]).(['Frequency' frequency_num]).zscore=zscore;
%
%

%   %% BASELINE HISTOGRAM: Pull out the responses to sound
%     for i=1:length(BOT_maps)
%       BOT_number = num2str(BOT_maps(i));
%       [trace_around_sound, baseline_histogram] = Sound_Responses_Widefield(data.([mouseID]).(['Tile' BOT_number]).All_Images_df_over_f,data.([mouseID]).(['Tile' BOT_number]).Sound_Time,data.([mouseID]).(['Tile' BOT_number]).Imaging_Timestamp);
%       data.([mouseID]).(['Tile' BOT_number]).All_Imaging_Around_Sound = trace_around_sound;
%       sz1 = size(baseline_histogram,1);
%       sz2 = size(baseline_histogram,2);
%       sz3 = size(baseline_histogram,3);
%       sz4 = size(baseline_histogram,4);
%       baseline_histogram_all_sounds = reshape(baseline_histogram, [sz1,sz2,sz3*sz4]);
%       figure;
%       histogram(baseline_histogram_all_sounds(1,1,:));% x pixel,y pixel,histogram of baselines per pixel
%       data.([mouseID]).(['Tile' BOT_number]).histogram = baseline_histogram_all_sounds;
%       data.([mouseID]).(['Tile' BOT_number]).histogram_by_sound = baseline_histogram;
%     end
%
%   %% MAX RESPONSES: Define max responses for each sound
%   for i=1:length(BOT_maps)
%       BOT_number = num2str(BOT_maps(i));
%          for l = 1:length(levels)
%            level_num = num2str(round(levels(l)));
%            figure;
%            histogram = [];
%            suptitle(sprintf('Average Responses at level %s',level_num));
%                for f = 1:length(frequencies)
%                     frequency_num = num2str(round(frequencies(f)));
%                     mean_response = data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num]).(['Frequency' frequency_num]).response;
%                     subplot(length(levels),length(frequencies),f);
%                     plot(squeeze(mean_response(50,50,:))); hold on;
%                     title(sprintf('%s',frequency_num));
%                     %spatial deconvolution for each pixel % add gCAMP6s
%                     %kernel here
%                     %temporal deconvolution filter for each pixel
%                     response_spatial_filter = deconvlucy(response_temporal_filter, psf); %spatial filter for each pixel
%                     plot(squeeze(response_spatial_filter(50,50,:)));
%
%                     for x=1:length(response_spatial_filter,1)
%                            for y = 1:length(response_spatial_filter,2)
%                             %    [peak, peak_loc] = max(response_spatial_filter(x,y,0.5*image_f:1.25*image_f),3);
%                             %    peak_response(x,y,:) = mean(response_spatial_filter(x,y,peak_loc-1:peak_loc+1);
%                                 data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num]).(['Frequency' frequency_num]).response=peak_response;
%                                 histogram = [response_spatial_filter(x,y,1:0.5*image_f) histogram];
%                            end
%                     end
%                end
%                end
%                 data.([mouseID]).(['Tile' BOT_number]).histogram=histogram;
%                end



%% PLOT MEAN RESPONSES: Plot mean responses to each frequency/level
% for i=1:length(BOT_maps)
%     BOT_number = num2str(BOT_maps(i));
%     BW = data.([mouseID]).(['Tile' BOT_number]).Window_Mask;
%     temp = num2cell(data.([mouseID]).(['Tile' BOT_number]).Window_Postion);
%     [x_min,x_max,y_min,y_max] = deal(temp{:});
%     BW=BW(y_min:y_max,x_min:x_max);
%     for l = 1:length(levels)
%         figure;
%         level_num = num2str(round(levels(l)));
%         suptitle(sprintf('Average Responses at level %s',level_num));
%         for f = 1:length(frequencies)
%             subplot(3,3,f);
%             frequency_num = num2str(round(frequencies(f)));
%             mean_response_temp = data.([mouseID]).(['Tile' BOT_number]).(['Level' level_num]).(['Frequency' frequency_num]).peakresponse;
%             A = medfilt2(mean_response_temp,[10 10]);
%             A(BW==0)=-0.1;
%             %   A=A(y_min:y_max,x_min:x_max);
%             myColorMap = jet(256);
%             
%             imagesc(A);  colormap jet;
%             
%             caxis([0 0.02]);
%             pbaspect([1 1 1]);
%             title(sprintf('%s',frequency_num));
%         end
%     end
%     
% end


%% MAP FREQUENCIES: Plot CFs

BW = mean(imageData.Cropped_Imaging_Data,3);

%make elliptical mask to eliminate values outside window
size_x = size(BW,1);
size_y = size(BW,2);
[col row] = meshgrid(1:size_x,1:size_y);
center_x = round(size_x/2);
center_y = round(size_y/2);
rad_x = round(size_x/2);
rad_y = round(size_y/2);
ellipse = (row-center_y).^2./rad_y^2+(col-center_x).^2./rad_x^2 <=1;
mask = ellipse';

%display window with masked edges
figure;
ax1=axes;
BW_mask = BW.*mask;
imagesc(BW_mask);
colormap(ax1,'gray');
hold on;

%display window with masked edges with CF overlay
ax2=axes;
CF = imageData.CF;
CF = medfilt2(CF,[5 5]);
CF_mask = CF.*mask;
CF_color_map = load('CFColormap.mat');
alpha(size_x,size_y)=0;
alpha(CF_mask>1)=0.3;
alpha(isnan(CF_mask))=0;
%alpha(isnan)=0;
im = imagesc(ax2,CF_mask,'alphadata',alpha);
colormap(ax2, CF_color_map.mymap);
caxis(ax2,[min(nonzeros(CF_mask)) max(nonzeros(CF_mask))]);
%alpha(0.3);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
colorbar;
%CF_color = CF_color_map.mymap;
%colormap(CF_color_map.mymap);
%RGB = ind2rgb(CF_mask, CF_color_map.mymap);
%CF_color_map.mymap);
%im2 = imagesc(CF_mask);
%image(RGB);
%pbaspect([1 1 1]);
%set(im2,'AlphaData',~isnan(CF_mask))
%colormap color_map.mymap;
%colormap(CF_color_map.mymap);

%colormap jet
%caxis([4 45]);
%image_CF = imagesc(CF_mask);
%set(h,'AlphaData',CF_mask);
%alpha(im,0.5);

%set(im,'AlphaData',~isnan(CF_mask))

figure;
im2 = imagesc(CF_mask);
%image(RGB);
pbaspect([1 1 1]);
set(im2,'AlphaData',~isnan(CF_mask))




%title(sprintf('CF Map'));
%hold off;


%     figure;
%     frequency_num = num2str(round(frequencies(f)));
%     CF = data.([mouseID]).(['Tile' BOT_number]).CF;s
%     A = medfilt2(CF,[2 2]);
%     A(BW==0)=-0.1;
%     %   A=A(y_min:y_max,x_min:x_max);
%     imagesc(A); colormap jet;
%     caxis([4 32]);


%% SAVE: save whole data file
save(['YD111219F3_test'], 'test', '-v7.3')
saved=2