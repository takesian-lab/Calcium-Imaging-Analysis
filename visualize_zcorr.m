function visualize_zcorr(block, plane)
% visualize_zcorr(block, plane)
%
% DOCUMENTATION IN PROGRESS
%
% This function allows you to preview the output of ops.zcorr from a single block
% 
% Argument(s): 
%   block (struct)
%   plane (int) for multiplane data, number (0,1,2...) of the plane you want to anaylze
%
% Returns:
%   
% 
% Notes:
%
%
% TODO: 
% Search 'TODO'

%% Setup

if nargin > 1
    multiplaneData = true;
    planeName = ['plane' num2str(plane)];
    
    if block.setup.XML.bidirectionalZ
        error('Code is not set up for bidirectional data yet')
    end
        
elseif nargin == 1 && isfield(block,'MultiplaneData')
    error('Please choose plane number: visualize_zcorr(block,plane)')
else
    multiplaneData = false;
end

max_drift = 10; %Number of planes that the imaging plane can drift in +/- Z

%Extract zcorr from block and establish best imaging plane
if multiplaneData
    zcorr_raw = block.zcorr.(planeName);
    nPlanes = block.ops.nplanes;
else
    zcorr_raw = block.zcorr;
end
[~, best_z_raw] = max(zcorr_raw,[],1);
imaging_plane = mode(best_z_raw);
A = imaging_plane - max_drift;
B = imaging_plane + max_drift;
if A < 1
    A = 1;
end
zcorr_cropped = zcorr_raw(A:B,:);
[~, best_z] = max(zcorr_cropped,[],1);
imaging_plane_cropped = mode(best_z);

%locomotor activity
if multiplaneData
    timestamp = block.timestamp.(planeName);
else
    timestamp = block.timestamp;
end

loco_time = block.locomotion_trace; %Matched to Bruker timestamps
loco_speed = block.loco_activity;
loco_speed(loco_speed < block.setup.constant.locoThresh) = 0; %Correct for floor of loco readout
        
if length(loco_speed) ~= length(loco_time)
    shorter_loco = min(length(loco_speed),length(loco_time));
    loco_time = loco_time(1:shorter_loco);
    loco_speed = loco_speed(1:shorter_loco);
    warning('locomotion_trace and loco_activity are not the same length')
end    
%% PLOT

if timestamp(end) < 200
    unit = 10;
else
    unit = 100;
end

%Make timestamp for plotting
plotLabels = 0:unit:(timestamp(end) - mod(timestamp(end),100));
plotInds = nan(size(plotLabels));
for i = 1:length(plotLabels)
    [~, plotInds(i)] = min(abs(timestamp-plotLabels(i)));
end

figure;
subplot(4,1,1)
imagesc(zcorr_raw)
hline(A, 'w')
hline(B, 'w')
set(gca, 'XTick', plotInds)
set(gca, 'XTickLabel', plotLabels)
ylabel('Z position')

subplot(4,1,2)
imagesc(zcorr_raw(A:B,:))
hline(imaging_plane_cropped, 'r')
set(gca, 'YTick', 1:5:(B-A)+1)
set(gca, 'YTickLabel', A:5:B+1)
set(gca, 'XTick', plotInds)
set(gca, 'XTickLabel', plotLabels)
ylabel('Cropped Z position')

subplot(4,1,3)
plot(best_z)
hline(imaging_plane_cropped, 'r')
set(gca, 'YTick', 1:5:(B-A)+1)
set(gca, 'YTickLabel', A:5:B+1)
set(gca, 'XTick', plotInds)
set(gca, 'XTickLabel', plotLabels)
ylabel('Best Z stack position')
xlim([1, length(best_z)])

subplot(4,1,4)
plot(loco_time, loco_speed, 'r')
xlim([timestamp(1) timestamp(end)])
xlabel('Time (s)')
ylabel('Loco activity (cm/s)')

suptitle(block.setup.block_supname)

%% XCORR OF LOCO ACTIVITY vs. Z STACK BEST POSITION

%Trim and resample Z or loco vector to match the other
%Z is typically longer during single plane recordings (because Tosca ends before BOT)
%loco is typically longer during multiplane recordings

%Trim timepoints that only exist in one vector
A = max([timestamp(1), loco_time(1)]);
Z = min([timestamp(end), loco_time(end)]);
keep_Z = intersect(find(timestamp >= A), find(timestamp <= Z));
keep_loco = intersect(find(loco_time >= A), find(loco_time <= Z));
timestamp_trimmed = timestamp(keep_Z);
best_z_trimmed = best_z(keep_Z);
loco_time_trimmed = loco_time(keep_loco);
loco_speed_trimmed = loco_speed(keep_loco);

%Find shorter vector
resample_necessary = true;
if length(timestamp_trimmed) < length(loco_time_trimmed)
    loco_longer = true;
    shorter_x = timestamp_trimmed;
    shorter_y = best_z_trimmed;
    longer_x = loco_time_trimmed;
    longer_y = loco_speed_trimmed;
elseif length(loco_time_trimmed) < length(timestamp_trimmed)
    loco_longer = false;
    shorter_x = loco_time_trimmed;
    shorter_y = loco_speed_trimmed;
    longer_x = timestamp_trimmed;
    longer_y = best_z_trimmed;   
else
    resample_necessary = false;
    xcorr_x = timestamp_trimmed;
    xcorr_loco = loco_speed_trimmed;
    xcorr_Z = -(best_z_trimmed - imaging_plane_cropped); %centered
end
    
%Resample longer vector
if resample_necessary
    new_longer_ind = nan(size(shorter_x));
    for i = 1:length(shorter_x)
        [~, new_longer_ind(i)] = min(abs(longer_x-shorter_x(i)));
    end
    new_longer_y = longer_y(new_longer_ind);
    xcorr_x = shorter_x;
    
    %Sanity check
%     figure
%     subplot(2,1,1)
%     plot(longer_x, longer_y)
%     subplot(2,1,2)
%     plot(new_longer_x, new_longer_y)

    if loco_longer
        xcorr_loco = new_longer_y;
        xcorr_Z = -(shorter_y - imaging_plane_cropped); %centered
    else
        xcorr_loco = shorter_y;
        xcorr_Z = -(new_longer_y - imaging_plane_cropped); %centered
    end
end

%Define lag time in samples
Fs = round(1/((xcorr_x(end) - xcorr_x(1))/length(xcorr_x)));
nSeconds = 10;
maxlag = Fs*nSeconds;

nNans = sum(isnan(xcorr_loco));
if nNans > 0
    xcorr_loco(isnan(xcorr_loco)) = 0;
    warning([num2str(nNans) ' NaNs replaced with zeros'])
end

%xcorr
[c, lags] = xcorr(xcorr_loco, xcorr_Z, maxlag, 'coeff');

%% FIGURE
figure;
subplot(4,1,1)
plot(xcorr_Z)
hline(0, 'r')
ylabel('Z difference')
xlim([1 length(xcorr_x)])

subplot(4,1,2)
plot(xcorr_loco, 'r')
xlabel('Time (s)')
ylabel('Loco activity')
xlim([1 length(xcorr_x)])

subplot(4,1,3:4)
plot(lags,c)
ylabel('Corr coeff')
set(gca, 'XTick', -maxlag:Fs:maxlag)
set(gca, 'XTickLabel', -maxlag/Fs:1:maxlag/Fs)
xlabel('Lag (s)')

suptitle(block.setup.block_supname)

%% XCORR OF CELL ACTIVITY vs. Z STACK BEST POSITION

if multiplaneData
    F7 = block.F.(planeName) - block.setup.constant.neucoeff*block.Fneu.(planeName);
    spks = block.spks.(planeName);
else
    F7 = block.F - block.setup.constant.neucoeff*block.Fneu;
    spks = block.spks;
end

%Compute DF/F
mean_F = mean(F7,2);% average green for each cell
F7_df_f = (F7 - mean_F)./mean_F;%(total-mean)/mean
%A = smooth(F7_df_f,10);

%Trim and resample if necessary
xcorr_F = F7_df_f(:,keep_Z);
xcorr_spks = spks(:,keep_Z);
if resample_necessary && ~loco_longer
    xcorr_F = xcorr_F(:,new_longer_ind);
    xcorr_spks = xcorr_spks(:,new_longer_ind);
end

%Compute lags and xcorr
cmat_loco_F = [];
cmat_Z_F = [];
cmat_loco_spks = [];
cmat_Z_spks = [];

for i = 1:size(xcorr_F,1)
    [c1, lags] = xcorr(xcorr_loco', xcorr_F(i,:), maxlag,  'coeff');
    [c2, lags] = xcorr(xcorr_Z, xcorr_F(i,:), maxlag, 'coeff');
    [c3, lags] = xcorr(xcorr_loco', xcorr_spks(i,:), maxlag,  'coeff');
    [c4, lags] = xcorr(xcorr_Z, xcorr_spks(i,:), maxlag, 'coeff');
    cmat_loco_F = [cmat_loco_F; c1];
    cmat_Z_F = [cmat_Z_F; c2];
    cmat_loco_spks = [cmat_loco_spks; c3];
    cmat_Z_spks = [cmat_Z_spks; c4];
end

%%

for f = 1:2
    
    if f == 1 %F
        cmat_loco = cmat_loco_F;
        cmat_Z = cmat_Z_F;
        activityType = 'F';
    elseif f == 2 %Spikes
        cmat_loco = cmat_loco_spks;
        cmat_Z = cmat_Z_spks;
        activityType = 'Spks';
    end

    figure

    subplot(3,2,1)
    plot(lags,cmat_loco)
    ylabel('Corr coeff')
    title(['xcorr(Loco, ' activityType ')'])

    subplot(3,2,2)
    plot(lags,cmat_Z)
    ylabel('Corr coeff')
    title(['xcorr(Z difference, ' activityType ')'])

    subplot(3,2,3)
    plot(lags,nanmean(cmat_loco))
    ylabel('Avg corr coeff')

    subplot(3,2,4)
    plot(lags,nanmean(cmat_Z))
    ylabel('Avg corr coeff')

    subplot(3,2,5)
    imagesc(cmat_loco)
    caxis([0 0.6])
    ylabel('Cells')
    set(gca, 'XTick', 0:Fs:2*maxlag)
    set(gca, 'XTickLabel', -maxlag/Fs:1:maxlag/Fs)
    xlabel('Lag (s)')

    subplot(3,2,6)
    imagesc(cmat_Z)
    caxis([0 0.6])
    ylabel('Cells')
    set(gca, 'XTick', 0:Fs:2*maxlag)
    set(gca, 'XTickLabel', -maxlag/Fs:1:maxlag/Fs)
    xlabel('Lag (s)')

    for i = 1:4
        subplot(3,2,i)
        set(gca, 'XTick', -maxlag:Fs:maxlag)
        set(gca, 'XTickLabel', -maxlag/Fs:1:maxlag/Fs)
        xlabel('Lag (s)')
    end
    suptitle(block.setup.block_supname)


    %Compute +/- peak of cross-correlogram
    [loco_peak, loco_peak_ind] = max(abs(cmat_loco),[],2);
    [Z_peak, Z_peak_ind] = max(abs(cmat_Z),[],2);
    loco_sign = double(cmat_loco(loco_peak_ind) > 0);
    Z_sign = double(cmat_Z(Z_peak_ind) > 0);
    loco_sign(loco_sign == 0) = -1;
    Z_sign(Z_sign == 0) = -1;
    loco_peak = loco_peak.*loco_sign;
    Z_peak = Z_peak.*Z_sign;

    figure
    scatter(loco_peak, Z_peak)
    xlabel(['xcorr(Loco, ' activityType ')'])
    ylabel(['xcorr(Z difference, ' activityType ')'])
    xlim([-0.4 0.4])
    ylim([-0.4 0.4])
    hline(0)
    vline(0)

    if ~isempty(find(abs(loco_peak) > 0.4)) || ~isempty(find(abs(Z_peak) > 0.4))
        warning('Some points exceed xy limits of scatterplot')
    end

    suptitle(block.setup.block_supname)
end

%% Plot figure of neural activity compared to Z and loco
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% SF = 0.25; %Shrinking factor for traces to appear more spread out (for visualization purposes)
% if multiplaneData
%     cellnum = block.cell_number.(planeName);
% else
%     cellnum = block.cell_number;
% end
% count = 1; %for staggering plot lines
% 
% for c = 1:20%length(cellnum)
%     current_cellnum = cellnum(c);
% 
%     subplot(4,4,1:8); hold on
% 
%     row_num = find(cellnum == current_cellnum);
%     cell_trace = xcorr_F7(row_num,:);%pull out the full trace for each cell
%     
%     plot(xcorr_x, cell_trace*SF + count,'LineWidth',1);
%     count = count + 1;
% end
% 
% suptitle(block.setup.block_supname)
% title('DF/F')
% xlim([0 xcorr_x(end)])
% xlabel('Time (s)')
% set(gca, 'YTick', [1:1:count-1])
% set(gca, 'YTickLabel', [cellnum(1:count-1)])
% ylabel('Cell')
% 
% %Plot best Z plane
% subplot(4,4,9:12); hold on %loco
% plot(xcorr_x, xcorr_Z);
% title('Best Z plane')
% ylabel('Best Z plane')
% xlim([0 xcorr_x(end)])
% xlabel('Time (s)')
% 
% %Plot locomotor activity
% subplot(4,4,13:16); hold on %loco
% plot(xcorr_x, xcorr_loco);
% title('Locomotor activity')
% ylabel('Activity (cm/s)')
% xlim([0 xcorr_x(end)])
% xlabel('Time (s)')

%% Plot histogram of F and Spikes

figure('units','normalized','outerposition',[0 0 1 1])

suptitle(block.setup.block_supname)

%F PSTH
subplot(4,4,1:4); hold on
psth = sum(xcorr_F);
area(xcorr_x, psth); 
title('DF/F PSTH')
xlim([0 xcorr_x(end)])
ylabel('DF/F')

%Spikes PSTH
subplot(4,4,5:8); hold on
psth = sum(xcorr_spks);
area(xcorr_x, psth);
title('Spikes PSTH')
xlim([0 xcorr_x(end)])
ylabel('Spikes')

%Plot best Z plane
subplot(4,4,9:12); hold on %loco
plot(xcorr_x, xcorr_Z);
title('Best Z plane')
ylabel('Best Z plane')
xlim([0 xcorr_x(end)])

%Plot locomotor activity
subplot(4,4,13:16); hold on %loco
plot(xcorr_x, xcorr_loco);
title('Locomotor activity')
ylabel('Activity (cm/s)')
xlim([0 xcorr_x(end)])
xlabel('Time (s)')
end