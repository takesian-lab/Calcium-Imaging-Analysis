clear all;

%% Noiseburst all cells from suite2p
%Anne Takesian - 2/22/2019
%updated Carolyn, compatible with Python version of Suite2p. Also does Red vs Green cell 7/23/19
%Updated Feb 2020, CGS - put most of the analysis into functions.
%Updated April 2020, MET - V3 created to load compiled blocks

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
    
    
    parameters.stim_protocol = 4; % widefield RF = 4, noiseburst ITI = 10
    imaging_chan = 'Ch2'; %was the data collected Ch1 or Ch2?
    BOT_start = [1];
    detrend_filter = [300 10];
    parameters.sort_loco =[1]; % 0 = all trials, 1 = non motor trials only
    
    
    
    %% Load Info.mat
    % Make setup and data structure out of all blocks that correspond to stim_protocol
    % Later we can also add other things like groups
    
    PC_name = getenv('computername');
    
    switch PC_name
        case 'RD0366' %Maryse
            info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            compiled_blocks_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledWidefieldBlocks';
            save_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
            info_filename = 'Info_widefield';
        case 'RD0332' %Carolyn
            info_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Info Sheets';
            %             compiled_blocks_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Daily Imaging';
            compiled_blocks_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\Compiled Blocks';
            save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\analyzed widefield\NxDE102320M2\Gcamp RF\adjusted trials';
            info_filename = 'Info_NxDE102320M2';
            
        case 'RD-6-TAK2' %Esther's computer
            info_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis';
            save_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\Widefield\YE083020F2';
            compiled_blocks_path = '\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\CompiledWidefieldBlocks';
            info_filename = 'Info_YE083020F2';
        case 'RD0386' %Wisam
            % INSERT PATHS HERE
            info_filename = 'Info';
        otherwise
            disp('Computer does not match known users')
    end
    
    cd(info_path)
    Info = importfile(info_filename);
    
    %Create data structure for files corresponding to stim_protocol
    [data] = fillSetupFromInfoTable_v2(Info, compiled_blocks_path, parameters.stim_protocol);
    data.setup.imaging_chan = imaging_chan;
    data.setup.BOT_start = BOT_start;
    
    %     data.setup.run_redcell = run_redcell;
    % p = gcp('nocreate');
    % if isempty(p)
    %     %%% If no pool, do not create new one.
    %     parpool('local',10);
    % end
end




%% register the image - 
% i registered the Tiffs in suite2p and then used those to load into the
% crop window function - can I somehow make this part of "compile blocks"
% eventually this code will run inside the 'crop_window' function. In the
% mean time - I am starting here from Full_Tile_Matrix (which is saved in
% the save path for this mouse)

% this (for now) has everything you need to crop the window. Dont run the
% next section
cd(save_path)
load('Full_Tile_Matrix.mat')
first_frame = Full_Tile_Matrix(:,:,1);
%make a 3d stack of the first frame?
fixed=zeros(size(Full_Tile_Matrix,1),size(Full_Tile_Matrix,2),size(Full_Tile_Matrix,3));%pre-allocate the matrix
for d = 1:size(Full_Tile_Matrix,3)
    fixed(:,:,d) = first_frame(:,:);
end

moving = Full_Tile_Matrix;
[optimizer, metric] = imregconfig('monomodal'); %can I optimize these settings? (currently default settings)
moving_reg = imregister(moving,fixed,'translation',optimizer,metric);
donereg =1




figure;
       reg_image = mat2gray(mean(moving_reg,3));
       imshow(reg_image);

        title(sprintf('Mean Window Tile; registered'));
         figure;
        h_im = imshow(reg_image); hold on;
        message = sprintf('Draw a circle and then hit the space bar,');
        uiwait(msgbox(message));
        h = imellipse;
        pause;
        BW = createMask(h);
        pos_window = getPosition(h); %returns [x min, y min, width, height]
        %csvwrite('[mouseID] 'window_position'']) to mouse folder
        parameters.x_min = pos_window(1)
        parameters.y_min = pos_window(2)
        parameters.x_max = pos_window(3)+parameters.x_min
        parameters.y_max = pos_window(4)+parameters.y_min
        Tile_Mean_ROI = reg_image;
        Tile_Mean_ROI(BW==0)=0;
        figure; imshow(Tile_Mean_ROI);title(sprintf('Masked Window'));
        pause;
        
        Tile_Mean_ROI=Tile_Mean_ROI(parameters.y_min:parameters.y_max,parameters.x_min:parameters.x_max);
        figure; imshow(Tile_Mean_ROI); title(sprintf('Cropped Window'));
        pause;
        parameters.Window_Postion = [parameters.x_min parameters.x_max parameters.y_min parameters.y_max];
        imageData.Cropped_Imaging_Data = moving_reg(parameters.y_min:parameters.y_max,parameters.x_min:parameters.x_max,:);
%         clear Full_Tile_Matrix
        clear Cropped
        clear Tile_Mean_ROI



 %% crop the window - I will take registerd images and add it to this function once I am done editing
[imageData,data,parameters,Full_Tile_Matrix] = crop_window(data,parameters);
% cd(save_path)
% save('Full_Tile_Matrix.mat','Full_Tile_Matrix','-v7.3');
% clear Full_Tile_Matrix;



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
    [amp_average_response,time_average_response] = max(mean_across_all_trials);
    estimated_peak = data.setup.FrameRate{i}*1.2; %the peak should be about 0.8sec after stim (1s)
    estimated_time = (time_average_response-estimated_peak)*(1/data.setup.FrameRate{i});
    parameters.adjusted_times_est = block.Sound_Time+estimated_time;
    
    %  repeat this average trace around sound with the estimated times
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
        all_trials(:,y) = smooth_All_Images(:,:);
    end
    
    mean_across_all_trials = mean(all_trials,2);
    figure;
    plot(1:length(average), squeeze(mean_across_all_trials),'k')
    title(sprintf('Mean Response after adjustment,all Trials', i));
    
    mean_across_all_trials = mean(all_trials,2);
    %       data.([mouseID]).(['Tile' BOT_number]).Average_Sound_Response = mean_across_all_trials;
    pause;
    % look at sound responses in high/low/med frequency trials
    if parameters.stim_protocol==4;
        figure;
        index_high_levels = find(data.([mouseID]).parameters.variable2>60);
        mean_across_all_high_trials = mean(all_trials(:,index_high_levels),2);
        plot(1:length(average), mean_across_all_high_trials,'r')
        
        hold on;
        index_mid_levels = find(data.([mouseID]).parameters.variable2>20 & data.([mouseID]).parameters.variable2<70);
        mean_across_all_mid_trials = mean(all_trials(:,index_mid_levels),2);
        plot(1:length(average), mean_across_all_mid_trials,'y')
        
        hold on;
        index_low_levels = find(data.([mouseID]).parameters.variable2<20);
        mean_across_all_low_trials = mean(all_trials(:,index_low_levels),2);
        plot(1:length(average), mean_across_all_low_trials,'g')
        title(sprintf('Mean Response by Level across all Trials from Tile %d', i));
        pause;
        
        
        figure;
        index_high_freq = find(data.([mouseID]).parameters.variable1>23);
        mean_across_all_highf_trials = mean(all_trials(:,index_high_freq),2);
        plot(1:length(average), mean_across_all_highf_trials,'r')
        
        hold on;
        index_mid_freq = find(data.([mouseID]).parameters.variable1>10 & data.([mouseID]).parameters.variable1<32);
        mean_across_all_midf_trials = mean(all_trials(:,index_mid_freq),2);
        plot(1:length(average), mean_across_all_midf_trials,'y')
        
        hold on;
        index_low_freq = find(data.([mouseID]).parameters.variable1<10);
        mean_across_all_lowf_trials = mean(all_trials(:,index_low_freq),2);
        plot(1:length(average), mean_across_all_lowf_trials,'g')
        title(sprintf('Mean Response by Freq across all Trials from Tile %d', i));
        pause;
        
        figure;
        % THis only runs if there are a lot of trials...
        
        %         mean_across_all_first_half_trials = mean(all_trials(:,1:300),2);
        %         plot(1:length(average), mean_across_all_first_half_trials,'r')
        %
        %         hold on;
        %         mean_across_all_first_half_trials = mean(all_trials(:,300:500),2);
        %         plot(1:length(average), mean_across_all_first_half_trials,'g')
        %
        %         hold on;
        %         mean_across_all_second_half_trials = mean(all_trials(:,500:720),2);
        %         plot(1:length(average), mean_across_all_second_half_trials,'y ')
        %         title(sprintf('Mean Response by Time across all Trials from Tile %d', i));
    end
    
    
    parameters.adjusted_times = parameters.adjusted_times_est;   %AT
end

%% check for memory
[loops] = memorycheck(imageData); 

%% DETREND: Grab a coffee - this will take approx. 2 hours (for full ReceptiveField, 45mins for new/reduced field)
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

% the data are divided into multiple 'tiles' to reduce memory issues. Tiles
% go from rostral to caudal when going from 1:length(Tiles) i.e.
% Tile1==most rostral
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
parameters.use_adjusted=1;

for i=1:length(data.setup.Imaging_sets)
    mouseID=data.setup.mousename{i};
    unique_block_name = data.setup.unique_block_names{i};
    block = data.([mouseID]).([unique_block_name]);
    
    if parameters.stim_protocol == 10;
        LevList=unique(data.([mouseID]).parameters.variable1); %noiseburst level
          parameters.levels=LevList;%store in a list for later
%         FreqList=1; % they are all the same sound
        n1=data.([mouseID]).parameters.variable1;%pull out variable#1 in order presented in experiment
        %         n2=data.([mouseID]).parameters.variable2;%pull out variable#2 in order presented in experiment
        for m=1:length(LevList)%loop through variable1 (frequency for TRF
            p=find(n1==LevList(m));%pull out a particular stimulus (Var #1) (i.e. 4kHz)
            [s]=p; %find specific stim types (i.e. 4khz, 60dB)
            parameters.stimIDX(m)={s};%stim index (Var#1xVar#2, i.e. freq x level)
        end
        
    else
        FreqList=unique(data.([mouseID]).parameters.variable1);%what is Variable#1 (freq)
        LevList=unique(data.([mouseID]).parameters.variable2);%what is Variable# (lev)
        n1=data.([mouseID]).parameters.variable1;%pull out variable#1 in order presented in experiment
        n2=data.([mouseID]).parameters.variable2;%pull out variable#2 in order presented in experiment
    
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
end

for i=1:length(data.setup.Imaging_sets)
    mouseID=data.setup.mousename{i};
    unique_block_name = data.setup.unique_block_names{i};
    block = data.([mouseID]).([unique_block_name]);
    
    
    if parameters.stim_protocol==10;
            for lv=1:length(parameters.levels);
                lvnum=num2str(parameters.levels(lv));
                idx=parameters.stimIDX{lv};
                parameters.loco_1.stimIDX{lv} = idx(find(block.active_trials(idx)==1));
                parameters.loco_0.stimIDX{lv} = idx(find(block.active_trials(idx)==0));
            end
            else
                for f=1:length(parameters.frequencies);
                    fnum=num2str(parameters.frequencies(f));
                    for lv=1:length(parameters.levels);
                        lvnum=num2str(parameters.levels(lv));
                        idx=parameters.stimIDX{f,lv};
                        parameters.loco_1.stimIDX{f,lv} = idx(find(block.active_trials(idx)==1));
                        parameters.loco_0.stimIDX{f,lv} = idx(find(block.active_trials(idx)==0));
                    end
                end
        end
    end


[traces,parameters]=sound_response_widefield_v3(parameters,data,All_Images_df_over_f);

%% pull out baseline and window
length_trial=size(traces.Tile1{1,1},3);
baseline=1:(block.setup.constant.baseline_length*data.setup.FrameRate{1});
% window=(length(baseline)+1):(data.setup.FrameRate{1}*block.setup.constant.response_window);
window=(length(baseline)+1):(data.setup.FrameRate{1}*1.5);
for ll=1:loops
    loop_num=num2str(ll);
    if parameters.stim_protocol ==10;
        for lv=1:length(parameters.levels);
            lvnum=num2str(parameters.levels(lv));
            if parameters.sort_loco ==0
                idx=parameters.stimIDX{lv};
                if lv ==1 & f==1 & ll==1
                    display('...creating baseline for all trials, noiseburst...')
                end
            else idx = parameters.loco_0.stimIDX{lv};
                if lv ==1 & f==1 & ll==1
                    display('...creating baseline for non-motor trials, noiseburst...')
                end
            end
            
            TF = isempty(idx); %the index of freq x level does not always have a value for every combination
            if TF==0
                base.(['Tile' loop_num]){lv}(:,:,:,:) = traces.(['Tile' loop_num]){lv}(:,:,baseline,:);
                %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
            end
        end
%       
        
    else
    for f=1:length(parameters.frequencies);
        fnum=num2str(parameters.frequencies(f));
        for lv=1:length(parameters.levels);
            lvnum=num2str(parameters.levels(lv));
            if parameters.sort_loco ==0
                idx=parameters.stimIDX{f,lv};
                if lv ==1 & f==1 & ll==1
                    display('...creating baseline for all trials...')
                end
            else idx = parameters.loco_0.stimIDX{f,lv};
                if lv ==1 & f==1 & ll==1
                    display('...creating baseline for non-motor trials...')
                end
            end
            
            TF = isempty(idx); %the index of freq x level does not always have a value for every combination
            if TF==0
                base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
                %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
            end
        end
        end
    end
end

%% mean across stimuli
loops=parameters.loops;
count=1;
for ll=1:loops
    loop_num=num2str(ll);
    if parameters.stim_protocol ==10;
                for lv=1:length(parameters.levels);
            if parameters.sort_loco ==0
                idx=parameters.stimIDX{lv};
                if lv ==1 & ll==1
                    display('...averaging across all trials, noiseburst...')
                end
            else idx = parameters.loco_0.stimIDX{lv};
                if lv ==1 & ll==1
                    display('...averaging across non-motor trials, noiseburst...')
                end
            end
            % occasionally, a stim will be empty, and this will correct
            % for when this occurs
            TF = isempty(idx);
            if TF == 0
                % do we want to pull out locomotor trials? If no (old
                % version of the code) run this step, if yes, see below.
                
                m=traces.(['Tile' loop_num]){lv};
                %mean response per stim
                mean_stim(:,:,:)=mean(m,4);
                avgTrace.(['Tile' loop_num]){lv}= mean_stim;
                clear mean_stim
                
            end
        end
   
    else
    for f=1:length(parameters.frequencies);
        for lv=1:length(parameters.levels);
            if parameters.sort_loco ==0
                idx=parameters.stimIDX{f,lv};
                if lv ==1 & f==1 & ll==1
                    display('...averaging across all trials...')
                end
            else idx = parameters.loco_0.stimIDX{f,lv};
                if lv ==1 & f==1 & ll==1
                    display('...averaging across non-motor trials...')
                end
            end
            % occasionally, a stim will be empty, and this will correct
            % for when this occurs
            TF = isempty(idx);
            if TF == 0
                % do we want to pull out locomotor trials? If no (old
                % version of the code) run this step, if yes, see below.
                
                m=traces.(['Tile' loop_num]){f,lv};
                %mean response per stim
                mean_stim(:,:,:)=mean(m,4);
                avgTrace.(['Tile' loop_num]){f,lv}= mean_stim;
                clear mean_stim
                
            end
        end
    end
    end
end


%
% parameters.avgBaseline=mean(tempBase,1);
% parameters.stdBaseline=std(tempBase,1);
clear tempBase idx b count ll loop_num lv m mean_base

%% what do the individual traces look like for a given tile?
%pick a tile to look at here, and set level to look at that intensity
%1=10...8=80dB
loop_num=num2str(1);


if parameters.stim_protocol==10;
  
    for lv=1:length(parameters.levels);
        if parameters.sort_loco ==0
            idx=parameters.stimIDX{lv};
        else idx = parameters.loco_0.stimIDX{lv};
       
        end
        % occasionally, a stim will be empty, and this will correct
        % for when this occurs
        TF = isempty(idx);
        if TF == 0
            m=traces.(['Tile' loop_num]){lv};
            size(traces.(['Tile' loop_num]){lv},4)
            figure;
            for mm = 1:size(traces.(['Tile' loop_num]){lv},4);
                sTrace = traces.(['Tile' loop_num]){lv}(:,:,:,mm);
                msTrace = squeeze(mean(mean(sTrace,1),2));
                plot(smooth(msTrace,5)); hold on;
                
            end
            
            
        end
    end

else
for f=1:length(parameters.frequencies);
    for lv=1;%1:length(parameters.levels);
        if parameters.sort_loco ==0
            idx=parameters.stimIDX{f,lv};
        else idx = parameters.loco_0.stimIDX{f,lv};
        end
        % occasionally, a stim will be empty, and this will correct
        % for when this occurs
        TF = isempty(idx);
        if TF == 0
            m=traces.(['Tile' loop_num]){f,lv};
            figure;
            for mm = 1:size(traces.(['Tile' loop_num]){f,lv},4);
                sTrace = traces.(['Tile' loop_num]){f,lv}(:,:,:,mm);
                msTrace = squeeze(mean(mean(sTrace,1),2));
                plot(smooth(msTrace,5)); hold on;
                
            end
            
            
        end
    end
end
end

%% put tiles back together
rejoin_tiles=[];
%loop through
if parameters.stim_protocol==10;
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        
        rejoin_tiles=[];
        for ll=1:loops
            loop_num=num2str(ll);
            
            mm=avgTrace.(['Tile' loop_num]){lv};
            rejoin_tiles=double(cat(1, rejoin_tiles, mm));
        end
        fieldName = matlab.lang.makeValidName(['db' numLV ]); %Replace invalid characters from fieldname, like -
        stimAverages.(fieldName)=rejoin_tiles;
        clear rejoin_tiles
    end

else

for f=1:length(parameters.frequencies)
    numF=num2str(round(parameters.frequencies(f)))
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        
        rejoin_tiles=[];
        for ll=1:loops
            loop_num=num2str(ll);
            
            mm=avgTrace.(['Tile' loop_num]){f,lv};
            rejoin_tiles=double(cat(1, rejoin_tiles, mm));
        end
        fieldName = matlab.lang.makeValidName(['kHz' numF ]); %Replace invalid characters from fieldname, like -
        stimAverages.(fieldName){lv}=rejoin_tiles;
        clear rejoin_tiles
    end
end
end
%    clear mm avgTrace

%% what do the stim averages look like? 
figure;
if parameters.stim_protocol ==10
  
    for lv=1:length(parameters.levels);
          
        numLV=num2str(parameters.levels(lv))
         figure;
        fieldName = matlab.lang.makeValidName(['db' numLV ])
        a1 = stimAverages.(fieldName);
        TF = isempty(a1);
        if TF ==0
            a2 = squeeze(mean(mean(a1,1),2));
            plot(smooth(a2)); %hold on
        end
    end

else
for f=1:length(parameters.frequencies)
    numF=num2str(round(parameters.frequencies(f)))
    subplot(2,ceil(length(parameters.frequencies)/2),f)
    
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        fieldName = matlab.lang.makeValidName(['kHz' numF ])
        a1 = stimAverages.(fieldName){lv};
        TF = isempty(a1);
        if TF ==0
            a2 = squeeze(mean(mean(a1,1),2));
            plot(smooth(a2)); hold on
        end
    end
end
end
%% convert to tif? and then store as individual file.
folder = save_path;
cd(folder)
tic
%noiseburst widefield
if parameters.stim_protocol==10
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
         fieldName = matlab.lang.makeValidName(['db' numLV]);
        idx=parameters.stimIDX{lv};
        fileName = matlab.lang.makeValidName(['avgStim_' numLV]);
        a1 = stimAverages.(fieldName);
        TF = isempty(a1);
        if TF ==0
            outputFileName= ([fileName 'db.tif']);
            for k = 1:size(stimAverages.(fieldName),3)
                imwrite(stimAverages.(fieldName)(:, :, k), outputFileName, 'WriteMode', 'append')
            end
        end
    end
else
    % RF and FMsweep widefield
for f=1:length(parameters.frequencies);
    toc
    numF=num2str(round(parameters.frequencies(f)));
    fieldName = matlab.lang.makeValidName(['kHz' numF]);
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
        idx=parameters.stimIDX{f,lv};
        fileName = matlab.lang.makeValidName(['avgStim_' numF 'kHz_' numLV]);
        a1 = stimAverages.(fieldName){lv};
        TF = isempty(a1);
        if TF ==0
            outputFileName= ([fileName 'db.tif']);
            for k = 1:size(stimAverages.(fieldName){lv},3)
                imwrite(stimAverages.(fieldName){lv}(:, :, k), outputFileName, 'WriteMode', 'append')
            end
        end
    end
end
end

%% normalize traces to mean 0 and std 1 - in progress....
% folder = save_path;
% cd(folder)
% d = dir([folder '/*.tif']);%extract tiffs
% total_stim = length(parameters.levels)*length(parameters.frequencies);
% accumBase = zeros(size(imageData.Cropped_Imaging_Data,1),size(imageData.Cropped_Imaging_Data,2),total_stim*length(baseline));
% for f=1:length(parameters.frequencies);
%     numF=num2str(round(parameters.frequencies(f)));
%     fieldName = matlab.lang.makeValidName(['kHz' numF]);
%     for lv=1:length(parameters.levels);
%         numLV=num2str(parameters.levels(lv));
%         % fname = (['SpatDenoise' numF 'khz' numLV 'db.tif']);
%         % fname = (['TempDenoise' numF 'khz' numLV 'db.tif']);
%         fileName = matlab.lang.makeValidName(['avgStim_' numF 'kHz_' numLV]);
%         fname= ([fileName 'db.tif']);
%          a1 = stimAverages.(fieldName){lv};
%         TF = isempty(a1);
%         if TF ==0
%         %        fname = (['avgStim' numF 'khz' numLV 'db.tif']);
%         info = imfinfo(fname);
%         num_images = numel(info)
%         for k = 1:num_images
%             AA = imread(fname, k);
%             stack(:,:,k)=AA(:,:);
%             stack=double(stack);
%         end
%         dimx = size(stack,1);
%         dimy = size(stack,2)
%         for x = 1:dimx
%             for y =  1:dimy
%                 normz = zscore(stack,0,3); %cgs left off here
%             end
%         end
%
%         %TODO: CGS - put the zscore here?
%
%
%         ResponseWindow{f,lv}=stack(:,:,window);
%         DFF0_mean{f,lv} = stack (:,:,:);
%         baseLoc = stack(:,:,baseline);
%         accumBase = cat(3, accumBase, baseLoc); % why is this such a large number?
%         end
%
%         %         base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
%         %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
%         clear stack AA
%     end
% end


%% %% find cumulative baseline and the response window
baseline=1:(block.setup.constant.baseline_length*data.setup.FrameRate{1});
window=(length(baseline)):(data.setup.FrameRate{1}*block.setup.constant.response_window);
folder = save_path;
% folder = 'C:\Anne';
cd(folder)
d = dir([folder '/*.tif']);%extract tiffs
count = 1;
% image = imageData.Cropped_Imaging_Data;
if parameters.stim_protocol==10
total_stim = length(parameters.levels);
else total_stim = length(parameters.levels)*length(parameters.frequencies);
end
accumBase = zeros(size(imageData.Cropped_Imaging_Data,1),size(imageData.Cropped_Imaging_Data,2),total_stim*length(baseline));


if parameters.stim_protocol==10
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
         fieldName = matlab.lang.makeValidName(['db' numLV]);
        % fname = (['SpatDenoise' numF 'khz' numLV 'db.tif']);
        % fname = (['TempDenoise' numF 'khz' numLV 'db.tif']);
        fileName = matlab.lang.makeValidName(['avgStim_' numLV]);
        fname= ([fileName 'db.tif']);
        TF = exist(fname);
        %          a1 = stimAverages.(fieldName){lv};
        %         TF = isempty(a1);
        if TF ==2 %sometimes there is a stim that is missing in the matrix of variables, this will only run if it exists
            %        fname = (['avgStim' numF 'khz' numLV 'db.tif']);
            info = imfinfo(fname);
            num_images = numel(info);
            for k = 1:num_images
                AA = imread(fname, k);
                stack(:,:,k)=AA(:,:);
                stack=double(stack);
            end
            ResponseWindow{lv}=stack(:,:,window);
            DFF0_mean{lv} = stack (:,:,:);
            baseLoc = stack(:,:,baseline);
            accumBase = cat(3, accumBase, baseLoc); % why is this such a large number?
        end
        
        %         base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
        %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
        clear stack AA
    end
else

for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)));
    fieldName = matlab.lang.makeValidName(['kHz' numF]);
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
        % fname = (['SpatDenoise' numF 'khz' numLV 'db.tif']);
        % fname = (['TempDenoise' numF 'khz' numLV 'db.tif']);
        fileName = matlab.lang.makeValidName(['avgStim_' numF 'kHz_' numLV]);
        fname= ([fileName 'db.tif']);
        TF = exist(fname);
        %          a1 = stimAverages.(fieldName){lv};
        %         TF = isempty(a1);
        if TF ==2 %sometimes there is a stim that is missing in the matrix of variables, this will only run if it exists
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
        end
        
        %         base.(['Tile' loop_num]){f,lv}(:,:,:,:) = traces.(['Tile' loop_num]){f,lv}(:,:,baseline,:);
        %             wind.(['Tile' loop_num]){f,lv}(:,:,:,:)=traces.(['Tile' loop_num]){f,lv}(:,:,window,:);
        clear stack AA
    end
end
end

stdBaseline=std(accumBase,0,3);
meanAccumBaseline = mean(accumBase,3);

%plot response window
% y=DFF0_mean{5,2};%RF
y=DFF0_mean{1,9};
%y=y(150:180,170:200,:);
y= squeeze(mean(mean(y,2),1));
x=1:length(y);
figure;
plot(x,y)

%what does a whole response window look like
% g = mean(ResponseWindow{5,2},3);
% %         CLIM = [0 350];
% imagesc(g);


%% plot all frequencies and amplitudes
figure;
n = 0;
if parameters.stim_protocol==10

    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        lv
        y = mean(ResponseWindow{lv},3);
        
        %         n = ((f-1)*length(parameters.frequencies))+lv
        n=n+1
        subplot(1,length(parameters.levels),n);
        imagesc(y);
        %    caxis([250 300]);
        title(sprintf('dB %d', round(parameters.levels(lv),2)));
        axis image;
        set(gca,'XTick',[], 'YTick', [])
    end

else
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)));
    fieldName = matlab.lang.makeValidName(['kHz' numF]);
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv))
        y = mean(ResponseWindow{f,lv},3);
        
        %         n = ((f-1)*length(parameters.frequencies))+lv
        n=n+1
        subplot(length(parameters.frequencies),length(parameters.levels),n);
        imagesc(y);
        %    caxis([250 300]);
        title(sprintf('Freq %d', round(parameters.frequencies(f),2)));
        axis image;
%         set(gca,'XTick',[], 'YTick', [])
lim = caxis
caxis([0 3])

    end
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
%     set(gca,'XTick',[], 'YTick', [])
lim = caxis
caxis([0 3])
end

%% zscore data
for f=1:length(parameters.frequencies);
    numF=num2str(round(parameters.frequencies(f)))
    figure;
    for lv=1:length(parameters.levels);
        numLV=num2str(parameters.levels(lv));
        mainResponse = ResponseWindow{f,lv};%temporally/spatially filtered response (2s post  sound)
        Df_f0 = DFF0_mean{f,lv};%temporally/spatially filtered mean response (entire trace i.e. baseline and 3s post sound)
        %z score of every frame of response window
        zscR = (bsxfun(@minus, mainResponse, (repmat(meanAccumBaseline,...
            [1 1 length(window)]))))./(repmat(stdBaseline,...
            [1 1 length(window)]));
        %         zscR = zscore(Df_f0);
        
        %zscore of every frame for entire trace
        zscResponse{f,lv}=zscR;
        zscS = (bsxfun(@minus,Df_f0, (repmat(meanAccumBaseline,...
            [1 1 size(Df_f0,3)]))))./(repmat(stdBaseline,...
            [1 1 size(Df_f0,3)]));
        zscSignal{f,lv} =  zscS;
        
        
        [ d1 d2 frames ] =size(zscR);%size of response window data
        
        
        maxResponse{f,lv} = max(mainResponse,[],3);%max of response window (not zscored)
        meanR{f,lv} = mean(mainResponse);%average of sound response window
        meanZresponse{f,lv} = mean(zscR(:,:,:),3);
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
%         if setup.parameters.stim_protocol==2
%             title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
%         end
%          axis([0 length(x_green) 0 70])
%
%     end
% end

%% determine CFs for each pixel
[dim1 dim2 dim3] = size(imageData.Cropped_Imaging_Data);
CF=NaN(dim1,dim2,1);
response_threshold = 2;


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


figure;
im2 = imagesc(CF_mask);
%image(RGB);
pbaspect([1 1 1]);
set(im2,'AlphaData',~isnan(CF_mask))

%% test "hot spots" for activity
message = sprintf('Draw a circle on completed map and then hit the space bar,');
uiwait(msgbox(message));
h = imellipse;
pause;
BW_2 = createMask(h);
pos_window = getPosition(h); %returns [x min, y min, width, height]
%csvwrite('[mouseID] 'window_position'']) to mouse folder
x_minSMwin = pos_window(1)
y_minSMwin = pos_window(2)
x_maxSMwin = pos_window(3)+x_minSMwin
y_maxSMwin = pos_window(4)+y_minSMwin
Tile_Mean_ROIsmWIN = imageData.Cropped_Imaging_Data;
Tile_Mean_ROIsmWINI(BW_2==0)=0;
figure; imshow(Tile_Mean_ROIsmWIN);title(sprintf('Masked small Window'));
pause;

Tile_Mean_ROIsmWIN=Tile_Mean_ROIsmWIN(y_minSMwin:y_maxSMwin,x_minSMwin:x_maxSMwin);
figure; imshow(Tile_Mean_ROIsmWIN); title(sprintf('Cropped small Window'));
pause;
Window_Postion_small = [x_minSMwin x_maxSMwin y_minSMwin y_maxSMwin];

%apply the new cropping to the mean tiles, and look to see if the
%tuning is as expected - note to Carolyn- this is still the code for
%accumBase, but it will serve as an easy template for what you are doing
%here
smallWindow_data = Full_Tile_Matrix(y_minSMwin:y_maxSmwin,x_minSMwin:x_maxSMwin,:);
folder = save_path;
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


%% Save data

% if loadPreviousData
%     cd(PathName) %Save in the same place you loaded data from
%     save([FileName(1:end-4) '_reload'])
% else
%     cd(save_path)
%     d = datestr(now,'yyyymmdd-HHMMSS');
%     save(['Data_' d '.mat'], 'data');
% end
cd(save_path)
save('imageData.mat','imageData','-v7.3');