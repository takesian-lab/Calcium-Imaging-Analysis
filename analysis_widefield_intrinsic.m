clear all;

%%
%Anne Takesian (AET) - 2/22/2019
%small updates CGS 4/12/19
%AET updated 6/20/19 to be consistent with Ed's analysis
%AET and CGS updated 4/20/20
%made compatible with "compile_blocks_from_info" July 2020 CGS

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

%CFcolormap
%define_sound_widefield.m
%behavior_RF_widefield.m
%crop_window.m
%locdetrend.m
%memorycheck.m
%Pixel_Detrend_Widefield_v2.m
%indexStimuli.m
%sound_response_widefield_v2.m
%constrained_foopsi_new.m
%vline.m


%To run, the data needs to be set up in the same path everytime.
%Define experiment-specific data below.

%% Load Data if it already exists, otherwise create new Data struct

loadPreviousData = 0;

if loadPreviousData
    %Load data
    [FileName,PathName] = uigetfile('.mat');
    load([PathName '/' FileName])
else
    
    %% Magic numbers and define what type of analysis you are doing
    %stim protocol code is:
    %noiseburst = 1
    %ReceptiveField = 2
    %FM sweep = 3
    %SAM = 6
    %widefield = 4
    %SAM freq = 6
    %Behavior = 7 and 8
    %Random H20 = 9
    %Noiseburst_ITI = 10
    %Random air puff = 11
    
    
    stim_protocol = 4; % this is the code for "widefield"
    
    % which channel was the data collected with?
    imaging_chan = 'Ch1';
    BOT_start = [1];
    detrend_filter = [300 10];
    
    
    
    %% Load Info.mat
    % Make setup and data structure out of all blocks that correspond to stim_protocol
    % Later we can also add other things like groups
    
    PC_name = getenv('computername');
    
    switch PC_name
        case 'RD0366' %Maryse
            info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            compiled_blocks_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocks';
            save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            info_filename = 'Info_NxDB092719M2';
        case 'RD0332' %Carolyn
            info_path = 'D:\2P analysis\2P local data\Carolyn';
            %             compiled_blocks_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Daily Imaging';
            compiled_blocks_path = 'Z:\Nick\Microglia Project\COVID_timeline\Compiled blocks\widefield_intrinsic';
            save_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Widefield';
            info_filename = 'Info_widefield_intrinsic_forCarolyn';
        case 'RD0386' %Wisam
            % INSERT PATHS HERE
            info_filename = 'Info';
        otherwise
            disp('Computer does not match known users')
    end
    
    cd(info_path)
    Info = importfile(info_filename);
    
    %Create data structure for files corresponding to stim_protocol
    [data] = fillSetupFromInfoTable_v2(Info, compiled_blocks_path, stim_protocol);
    data.setup.imaging_chan = imaging_chan;
    data.setup.BOT_start = BOT_start;
    
    %     data.setup.run_redcell = run_redcell;
    p = gcp('nocreate');
    if isempty(p)
        %%% If no pool, do not create new one.
        parpool('local',4);
    end
end
%% crop the window
[imageData,data,parameters] = crop_window(data);


%% anne's "old code" to test for window quality and check sound stimulation
% we are still trying to determine whether the sounds are off by
% approximately 500ms. Set this value here, to visually test the shift.
% (the adjust_factor is in seconds)
adjust_factor = 0.5; % in seconds


for i=1:length(data.setup.Imaging_sets)
    mouseID=data.setup.mousename{i};
    unique_block_name = data.setup.unique_block_names{i};
    block = data.([mouseID]).([unique_block_name]);
    
    Full_Tile_Mean = mean(mean(imageData.Cropped_Imaging_Data,1),2);
    Full_Tile_Mean_Detrend = locdetrend(Full_Tile_Mean(1,1,:),1,detrend_filter); %
    timestamp = block.timestamp; %change this number
   block.adjusted_times(:)= block.Sound_Time(:)-adjust_factor; % this '5' accounts for an apparent 500 ms shift in the data
    figure;
    hold on;
    
    % plot trend, detrend, trace, and adjusted times
    fo =squeeze(Full_Tile_Mean)-squeeze(Full_Tile_Mean_Detrend);
    f = squeeze(Full_Tile_Mean);
    plot(timestamp(1:length(timestamp)), squeeze(Full_Tile_Mean_Detrend),'c');
    plot(timestamp(1:length(timestamp)), smooth(f-1750,3),'b');
    plot(timestamp(1:length(timestamp)), fo,'m');
    
    h=vline(block.adjusted_times(:),'k'); hold on % adjusted sound times
    title(sprintf('Mean Response, 500ms adjusted,all Trials', i));
    %xlim([0 1000])
    
    % plot sound times determined from 'compile_blocks_from_info'
    figure;
    hold on;
    fo =squeeze(Full_Tile_Mean)-squeeze(Full_Tile_Mean_Detrend);
    f = squeeze(Full_Tile_Mean);
    plot(timestamp(1:length(timestamp)), squeeze(Full_Tile_Mean_Detrend),'c');
    plot(timestamp(1:length(timestamp)), smooth(f-1750,3),'b');
    plot(timestamp(1:length(timestamp)), fo,'m');
    
    h=vline(block.Sound_Time(:),'k'); hold on %sound times
    title(sprintf('Mean Response, Sound times (not adjusted),all Trials', i));
    %xlim([0 1000])
    pause;
    
    %  average trace around sound
    
    total_average = 1;
    Sound_Time = block.Sound_Time(:);
    
     for y = 1:size(Sound_Time,1)
        before = Sound_Time(y)-block.setup.constant.baseline_length;
        after = Sound_Time(y)+block.setup.constant.after_stim;
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
    pause;
    
    
    figure;
   plot(1:length(average), squeeze(mean_across_all_trials),'k')
    title(sprintf('Mean Response, all Trials', i));
    
    %AT added 4/15/20 to center window around peak response across window
    clear all_trials
    [amp_average_response,time_average_response] = min(mean_across_all_trials);% its min for intrinsic
    estimated_peak = data.setup.FrameRate{i}*5; %the peak min response should be about 5 sec after stim (1s)
    estimated_time = (time_average_response-estimated_peak)*(1/data.setup.FrameRate{i});
    parameters.adjusted_times_est = block.Sound_Time+estimated_time;
    
    %  repeat this average trace around sound
     clear all_trials baseline_mean average smooth_All_Images
    total_average = 1;
    Sound_Time = parameters.adjusted_times_est(:);
    
    for y = 1:size(Sound_Time,1)
         before = Sound_Time(y)-block.setup.constant.baseline_length;
        after = Sound_Time(y)+block.setup.constant.after_stim;
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
    title(sprintf('Mean Response after adjustment,all Trials', i));
    
    mean_across_all_trials = mean(all_trials,2);
    %       data.([mouseID]).(['Tile' BOT_number]).Average_Sound_Response = mean_across_all_trials;
    pause;
    
    figure;
    index_high_levels = find(data.([mouseID]).parameters.variable2>60);
    mean_across_all_high_trials = mean(all_trials(:,index_high_levels),2);
    plot(1:length(average), mean_across_all_high_trials,'r')
    
    hold on;
    %       index_mid_levels = find(parameters.level_list>20 & parameters.level_list<70);
    %       mean_across_all_mid_trials = mean(all_trials(:,index_mid_levels),2);
    %       plot(1:length(average), mean_across_all_mid_trials,'y')
    
    hold on;
    %       index_low_levels = find(parameters.level_list<20);
    %       mean_across_all_low_trials = mean(all_trials(:,index_low_levels),2);
    %       plot(1:length(average), mean_across_all_low_trials,'g')
    %        title(sprintf('Mean Response by Level across all Trials from Tile %d', i));
    pause;
    
    figure;
    index_high_freq =find(data.([mouseID]).parameters.variable1>23);
    mean_across_all_highf_trials = mean(all_trials(:,index_high_freq),2);
    plot(1:length(average), mean_across_all_highf_trials,'r')
    
    hold on;
    index_mid_freq = find(data.([mouseID]).parameters.variable1>10 & data.([mouseID]).parameters.variable1<32);
    mean_across_all_midf_trials = mean(all_trials(:,index_mid_freq),2);
    plot(1:length(average), mean_across_all_midf_trials,'y')
    
    hold on;
    index_low_freq =find(data.([mouseID]).parameters.variable1<10);
    mean_across_all_lowf_trials = mean(all_trials(:,index_low_freq),2);
    plot(1:length(average), mean_across_all_lowf_trials,'g')
    title(sprintf('Mean Response by Freq across all Trials from Tile %d', i));
    pause;
    
    %       figure;
    %
    %       mean_across_all_first_half_trials = mean(all_trials(:,1:300),2);
    %       plot(1:length(average), mean_across_all_first_half_trials,'r')
    %
    %       hold on;
    %        mean_across_all_first_half_trials = mean(all_trials(:,300:500),2);
    %       plot(1:length(average), mean_across_all_first_half_trials,'g')
    %
    %       hold on;
    %       mean_across_all_second_half_trials = mean(all_trials(:,500:868),2);
    %       plot(1:length(average), mean_across_all_second_half_trials,'y ')
    %       title(sprintf('Mean Response by Time across all Trials from Tile %d', i));
    
    
    parameters.adjusted_times = parameters.adjusted_times_est;   %AT
end
%% check for memory
[loops] = memorycheck(imageData);

%% DETREND: Grab a coffee - this will take approx. 2 hours
t = cputime;
parameters.loops=loops;
for i=1:length(data.setup.Imaging_sets)
    %     BOT_number = num2str(data.setup.Imaging_sets(i));
    Cropped = imageData.Cropped_Imaging_Data;
    [All_Images_df_over_f] = Pixel_Detrend_Widefield_v2(Cropped,loops);
end
clear Cropped
%
e = cputime-t
clear t e i BOT_number

%% view df_over_f
tile_to_view = [2];
tilenum = num2str(tile_to_view);

for i=1:length(data.setup.Imaging_sets)
    %mean of all pixels across time
    %     BOT_number = num2str(setup.BOT_maps(i));
    FullTile_df= squeeze(mean(mean(All_Images_df_over_f.(['Tile' tilenum]),1),2));
    mouseID=data.setup.mousename{i};
    unique_block_name = data.setup.unique_block_names{i};
    block = data.([mouseID]).([unique_block_name]);
    timestamp=block.timestamp;
    Sound_Time=parameters.adjusted_times(:); %estimated sound times
    figure;
    plot(timestamp,FullTile_df,'b');
    hold on; vline(Sound_Time);
    title(sprintf('Mean Response across all Trials from Tile %d', tile_to_view));
end
clear BOT_number FullTile_df i Sound_Time timestamp

%% create frequency and level indicies and find responses to sound across stim
% [parameters] = indexStimuli(parameters,setup);


%we are still testing out why the stimulus-evoked peaks seem off from where
%they are expected. Above, Anne wrote a code to determine the response
%window based upon peak location. This "use_adjusted" variable just lets
%you choose whether you want to use the original or adjusted times in the
%rest of the analysis
parameters.use_adjusted=0;


for i=1:length(data.setup.Imaging_sets)
    mouseID=data.setup.mousename{i};
    unique_block_name = data.setup.unique_block_names{i};
    block = data.([mouseID]).([unique_block_name]);
    FreqList=unique(data.([mouseID]).parameters.variable1);%what is Variable#1
    LevList=unique(data.([mouseID]).parameters.variable2);%what is Variable#
    parameters.frequencies=FreqList;%store in a list for later
    parameters.levels=LevList;%store in a list for later
    n1=data.([mouseID]).parameters.variable1;%pull out variable#1 in order presented in experiment
    n2=data.([mouseID]).parameters.variable2;%pull out variable#2 in order presented in experiment
    for m=1:length(FreqList)%loop through variable1 (frequency for TRF
        p=find(n1==FreqList(m));%pull out a particular stimulus (Var #1) (i.e. 4kHz)
        for q=1:length(LevList)
            r=find(n2==LevList(q));%pull out a particular stimulus (Var #2) (i.e. 60dB)
            [s]=intersect(p,r); %find specific stim types (i.e. 4khz, 60dB)
            parameters.stimIDX(m,q)={s};%stim index (Var#1xVar#2, i.e. freq x level)
        end
    end
end



[traces]=sound_response_widefield_v3(parameters,data,All_Images_df_over_f);
%% pull out baseline and window
length_trial=size(traces.Tile1{1,1},3);
baseline=1:(0.5*data.setup.FrameRate{1});% TODO: magic number
window=(length(baseline)+1):(data.setup.FrameRate{1}*2);

for ll=1:loops
    loop_num=num2str(ll)
    for f=1:length(parameters.frequencies);
        fnum=num2str(parameters.frequencies(f));
        for lv=1:length(parameters.levels);
            lvnum=num2str(parameters.levels(lv));
            idx=parameters.stimIDX{f,lv};
            TF = isempty(idx); %the index of freq x level does not always have a value for every combination
            if TF==0
                base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
                %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
            end
        end
    end
end

%% mean across stimuli
loops=parameters.loops;
count=1;
for ll=1:loops
    loop_num=num2str(ll);
    for f=1:length(parameters.frequencies);
        for lv=1:length(parameters.levels);
            idx=parameters.stimIDX{f,lv};
            TF = isempty(idx);
            if TF == 0
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
          fieldName = matlab.lang.makeValidName(['kHz' numF ]); %Replace invalid characters from fieldname, like -
          stimAverages.(fieldName){lv}=rejoin_tiles;
          clear rejoin_tiles
      end
  end
%    clear mm avgTrace

%% what do the stim averages look like?
figure;
for f=1:length(parameters.frequencies)
    numF=num2str(round(parameters.frequencies(f)))
    subplot(2,4,f)
    
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        fieldName = matlab.lang.makeValidName(['kHz' numF ])
        a1 = stimAverages.(fieldName){lv};
        lv
        TF = isempty(a1);
        if TF ==0
            a2 = squeeze(mean(mean(a1,1),2));
            plot(smooth(a2)); hold on
        end     
    end
end
%% convert to tif? and then store as individual file.
folder = 'D:\2P analysis\2P local data\Carolyn\Widefield\VxDD033120F2_gcamp';

%     folder = 'C:\Anne\';
 cd(folder)
    tic
    for f=1:length(parameters.frequencies);
        toc
        numF=num2str(round(parameters.frequencies(f)));
        fieldName = matlab.lang.makeValidName(['kHz' numF]);
        for lv=1:length(parameters.levels);
            numLV=num2str(parameters.levels(lv));
            idx=parameters.stimIDX{f,lv};
            fileName = matlab.lang.makeValidName(['avgStim_' numF 'kHz_' numLV]);
            outputFileName= ([fileName 'db.tif']);
            for k = 1:size(stimAverages.(fieldName){lv},3)
                imwrite(stimAverages.(fieldName){lv}(:, :, k), outputFileName, 'WriteMode', 'append')
            end
        end
    end

%% %% find cumulative baseline and the response window
baseline=1:(0.5*data.setup.FrameRate{1});
window=(length(baseline)):(data.setup.FrameRate{1}*3);
folder = 'D:\2P analysis\2P local data\Carolyn\Widefield\VxDD033120F2_gcamp';
% folder = 'C:\Anne';
cd(folder)
d = dir([folder '/*.tif']);%extract tiffs
count = 1;
image = imageData.Cropped_Imaging_Data;
total_stim = length(parameters.levels)*length(parameters.frequencies);
accumBase = zeros(size(image,1),size(image,2),total_stim*length(baseline));
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
    fieldName = matlab.lang.makeValidName(['kHz' numF]);
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
       % fname = (['SpatDenoise' numF 'khz' numLV 'db.tif']);
       % fname = (['TempDenoise' numF 'khz' numLV 'db.tif']);
       fileName = matlab.lang.makeValidName(['avgStim_' numF 'kHz_' numLV]);
       fname= ([fileName 'db.tif']);
%        fname = (['avgStim' numF 'khz' numLV 'db.tif']);
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
y=DFF0_mean{4,1};
%y=y(150:180,170:200,:);
y= squeeze(mean(mean(y,2),1));
x=1:length(y);
figure;
plot(x,y)

%what does a whole response window look like
 g = mean(ResponseWindow{5,2},3);
%         CLIM = [0 350];
   imagesc(g);


%% plot all frequencies and amplitudes
figure;
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
     fieldName = matlab.lang.makeValidName(['kHz' numF]);
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


 %% Save data
% 
% if loadPreviousData
%     cd(PathName) %Save in the same place you loaded data from
%     save([FileName(1:end-4) '_reload'])
% else
%     cd(save_path)
%     d = datestr(now,'yyyymmdd-HHMMSS');
%     save(['Data_' d '.mat'], 'data');
% end