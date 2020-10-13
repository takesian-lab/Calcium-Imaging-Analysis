function [loco_cor_cell] = visualize_locomotor(block)
% Carolyn 5/29/2020


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
    run_redcell = 0;
    std_level = 1.5;
    std_level_byStim = 1.5;



%% loco_trace
    
%DF/F for fluorescent trace
            F = block.F; %fluoresence
            Fneu = block.Fneu; % neuropil
            F7 = F - block.setup.constant.neucoeff * Fneu; % neuropil corrected trace
            
            for k = 1:size(F7,1)
                baseline_Fo = mean(F7(k,:),2);
                DF_F0(k,:) = F7(k,:)-baseline_Fo./baseline_Fo;
             locdetrend_temp = locdetrend(F7(k,:),1,[300 10]);
             DF_F0_dt(k,:) = locdetrend_temp'./(F7(k,:)-locdetrend_temp');
%              DF_F0_dt(k,:) = F7(k,:);
            end
            
            locoTrace = block.loco_activity;
            locoTime = block.loco_times;
            redcell = block.redcell;
           


%% plot detrended data...(ten traces)....

% note:the loco-trace is on a slightly different scale as the Ca-trace.
% This will be corrected below
count = 0;
subplot (2,1,1)
for i = 1:size(DF_F0_dt,1) % plot first 10 traces
    timestamp = block.timestamp;
    y = DF_F0_dt(i,:)+count;
    count = count+3;
    plot(timestamp,smooth(y,10)), hold on
end
subplot(2,1,2)
plot(locoTime,locoTrace)


%% crop the Ca++ trace to only look at when we have locomotor data
locoStart = locoTime(1);
locoEnd = locoTime(end);
timestamp = block.timestamp;
[c closest_frame_locoStart] = min(abs(timestamp(:,1)-locoStart));
[c closest_frame_locoEnd] = min(abs(timestamp(:,1)-locoEnd));

% get traces that align with locoStart to locoEnd
for i = 1:size(DF_F0_dt,1)
    DF_F0_loco(i,:) = DF_F0_dt(i,closest_frame_locoStart:closest_frame_locoEnd);
end

%% upsample loco trace
% remove NaNs by interpolation
% paint = inpaint_nans(locoTrace); 

% upsample by interplation to make same size as Ca++ trace
v = 1:numel(locoTrace);
ve = 1:size(DF_F0_loco,2);
vr = linspace(min(v), max(v), size(DF_F0_loco,2));
locoTrace_up = interp1(v, locoTrace, vr);

% threshold the locomotor trace
locoTrace_up(locoTrace_up<block.setup.constant.locoThresh)=0;

%% look at detrended gcamp and locomotor trace on same graph, same scale

% get traces that align with locoStart to locoEnd

% count = 0;
% subplot (2,1,1)
% for i =1;%:20;%size(DF_F0,1)
%     y = DF_F0_loco(i,:)+count;
%     plot(smooth(y,10)); hold on
%     count = count+3;    
% end
% subplot(2,1,2)
%      plot(smooth(locoTrace_up(1,:),10))
    
%% Plot locomotor trace (now, thresholded), against gcamp signal (detrended)
% Make the first axes in its figure
% figure;
% x = locoTrace_up;
% x2 = 1:numel(locoTrace_up);
% for i = 1:size(DF_F0_loco,1)
%     y= DF_F0_loco(i,:);
%   
%     if redcell(i) ==1
%          subplot(8, 5, i);
%    
% %     scatter(x,y,'ro');
%     mygca(i) = gca;
% %     yb = scatstat1(x,y,15,@median);
% %     plot(x,yb,'ro')
% % x = 1:50; 
% % y = -0.3*x + 2*randn(1,50); 
% p = polyfit(x,y,1); 
% f = polyval(p,x); 
% plot(x,y,'o',x,f,'-') 
% legend('data','linear fit') 
%     else
%          hAx(i)=subplot(8, 8, i);
% %          scatter(x,y,'go');
%            mygca(i) = gca;
%            p = polyfit(x,y,1); 
% f = polyval(p,x); 
% plot(x,y,'o',x,f,'-') 
% legend('data','linear fit') 
%     end
% end
% yl = cell2mat(get(mygca, 'Ylim'))
% ylnew = [min(yl(:,1)) max(yl(:,2))];
% set(mygca, 'Ylim', ylnew)

% scatterhist(locoTrace_up,DF_F0_loco(1,:));

%% 
% % find data associated with running versus data associated with not running
% 
% noloco_boutIDX =locoTrace_up==0;
% loco_boutIDX = ~noloco_boutIDX;
% 
% for i = 1:size(DF_F0_loco,1)
%     mean_loco(i) = mean(DF_F0_loco(i,loco_boutIDX));
%     mean_noloco(i) =mean(DF_F0_loco(i,noloco_boutIDX));
% end
% 
% mLoco = mean(mean_loco);
% semLoco = std(mean_loco)./sqrt(numel(mean_loco));
% 
% 
% % scatter(x,mean_loco); hold on
% % scatter(x2, mean_noloco)
% 
% mnoLoco =  mean(mean_noloco);
% semnoLoco = std(mean_noloco)./sqrt(numel(mean_noloco));
%% do the correlation  - not actually z-scored, I just didnt change the name of the variables
figure;

z_locoTrace = locoTrace_up;
% z_DF_F0 = (DF_F0_loco);
z_locoDiff = zscore(diff(locoTrace_up));


count = 1;
for i = 1:size(DF_F0_loco,1)
%     if redcell(i)==1
z_DF_F0(i,:) = (DF_F0_loco(i,:));
[coef,lags] = xcorr(z_DF_F0(i,:),z_locoTrace(:,:),'coeff',60); %magic number
% [coef_diff,lags_diff] = xcorr(z_locoDiff,z_DF_F0(i,:));
[max_cor,idx] = max(coef);
cor_idx(i) = idx;
% [max_cor_diff,idx_diff] = max(coef);
max_c(i) = max_cor;
% max_c_diff(i) = max_cor_diff;
max_lag(i) = lags(idx);
[r,p] = corrcoef(z_DF_F0(i,:), z_locoTrace(:,:));
pval{i}=p;
if pval{i}(1,2)<0.05
    loco_cor_cell(count) = block.cell_number(i)
    count = count+1;
end

rr(i,:)=r(1,2);
% pp(i,:)=p(1,2)

figure;


plot(lags,coef,'-b')
title(['cell number ',num2str(block.cell_number(i))])
% max_lag_diff = lags_diff(idx_diff);
%     end
end
%% trace in response to loco onset/offset
% loco_trace_filt = smooth(locoTrace_up,10);
% for i = 2:length(loco_trace_filt)
%     if loco_trace_filt(i-1) ==0 && loco_trace_filt(i)>0 
%         onsetIDX(i)=1;
%     end
%     if loco_trace_filt(i-1)>0 && loco_trace_filt(i)==0
%         offsetIDX(i)=1;
%     end
% end
%  onFrames = find(onsetIDX);
%  offFrames = find(offsetIDX);
% % define the onset window
% pre_onset = 2*30;
% onset_end = 2*30;
% on_trace_wind = zeros(size(onFrames,1),pre_onset+onset_end+1);
% figure;

% for k = 1:size(DF_F0_loco,1)
%     for i = 1:length(onFrames)
%         if onFrames(i)>pre_onset
%             on_trace_wind(i,:) = DF_F0_loco(k,onFrames(i)-pre_onset:onset_end+onFrames(i));
% %             plot(smooth(on_trace_wind(i,:),10)), hold on
%         end
%        x=1:size(onFrames(i)-pre_onset:onset_end+onFrames(i),2);
%     end
%      
%         a_trace = mean(on_trace_wind,1);
%         onset_traces(k,:) = a_trace;
%         plot(x,smooth(onset_traces(k,:),10));hold on
% end
% % figure;
% % x=1:size(onFrames(i)-pre_onset:onset_end+onFrames(i),2);
% % mean_onset = mean(on_trace_wind,1);
% % sem_onset = std(on_trace_wind,1)./sqrt(size(on_trace_wind,1));
% % shadedErrorBar(x,smooth(mean_onset,10),smooth(sem_onset,10),'lineProps','-b');
% 
% % plot(mean_onset)
% 
% %% find mean loco onset traces
% 
% 
% %%
% figure;
% scatter(lags,coef)
%% plot loco vs non loco trials

actIDX = block.active_trials;
stimTrace = block.aligned_stim.F7_stim;
baseline = block.setup.constant.baseline_length * block.setup.framerate;
base_frames = stimTrace(:,:,1:baseline);
mean_base_frames = mean(base_frames,3);
for i = 1:size(stimTrace,1)
    for j = 1:size(stimTrace,2)
       ftemp = stimTrace(i,j,:);
       btemp = mean_base_frames(i,j);
       dftemp = ftemp-btemp./btemp;
       df_f(i,j,:) = dftemp;
    end
end

% now pull out the loco/non loco trials
acttrial = find(actIDX==1);
inacttrial =find(actIDX~=1);
active_traces = df_f(:,acttrial,:);
inactive_traces = df_f(:,inacttrial,:);

%avg across cells
for i = 1:size(active_traces,1) 
    a = squeeze(mean(active_traces(i,:,:),2));
    a_err = squeeze(std(active_traces(i,:,:),0,2)/sqrt(size(active_traces,2)));
    b = squeeze(mean(inactive_traces(i,:,:),2));
    b_err = squeeze(std(inactive_traces(i,:,:),0,2)/sqrt(size(inactive_traces,2)));
    x = 1:size(active_traces,3);
    figure; 
    shadedErrorBar(x,smooth(a,10),a_err,'lineprops','-g'); hold on
    shadedErrorBar(x,smooth(b,10),b_err,'lineprops','-m'); hold on
    vline(baseline)
    
    active_cell(i,:,:) =a; 
    inactive_cell(i,:,:) = b;
    title(['cell number ',num2str(block.cell_number(i))])
    legend('active','inactive')
end



end
