function visualize_zcorr(block, plane)
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
else
    multiplaneData = false;
end

max_drift = 10; %Number of planes that the imaging plane can drift in +/- Z

%Extract zcorr from block and establish best imaging plane
if multiplaneData
    zcorr_raw = block.ops.zcorr.(planeName);
    nPlanes = numel(fieldnames(block.ops.zcorr));
else
    zcorr_raw = block.ops.zcorr;
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
timestamp = block.timestamp;
if multiplaneData
    planeInd = [1:nPlanes:length(timestamp)] + plane;
    planeInd(planeInd > length(timestamp)) = [];
    timestamp = timestamp(planeInd);
end

loco_time = block.locomotion_trace;
loco_speed = block.loco_activity;
        
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
xlim([timestamp(1) timestamp(end)]) %loco will be shorter because tosca ends before PV
xlabel('Time (s)')
ylabel('Loco activity (cm/s)')

suptitle(block.setup.block_supname)

%% XCORR OF LOCO ACTIVITY vs. Z STACK BEST POSITION

%Only take Z information where loco activity is available
%And resample points in Z vector to match loco vector

best_z_index = nan(size(loco_time));
for i = 1:length(loco_time)
    [~, best_z_index(i)] = min(abs(timestamp-loco_time(i)));
end
best_z_loco = best_z(best_z_index);

%Center and normalize
best_z_loco_centered = -(best_z_loco - imaging_plane_cropped);

%Define lag time in samples
loco_sampling_rate = round(1/((loco_time(end) - loco_time(1))/length(loco_time)));
nSeconds = 10;
maxlag = loco_sampling_rate*nSeconds;

%xcorr
[c, lags] = xcorr(loco_speed, best_z_loco_centered, maxlag, 'coeff');

%% FIGURE
figure;
subplot(4,1,1)
plot(best_z_loco_centered)
hline(0, 'r')
ylabel('Z difference')
xlim([0 length(loco_time)])

subplot(4,1,2)
plot(loco_speed, 'r')
xlabel('Time (s)')
ylabel('Loco activity')
xlim([0 length(loco_time)])

subplot(4,1,3:4)
plot(lags,c)
ylabel('Corr coeff')
set(gca, 'XTick', -maxlag:loco_sampling_rate:maxlag)
set(gca, 'XTickLabel', -maxlag/loco_sampling_rate:1:maxlag/loco_sampling_rate)
xlabel('Lag (s)')

suptitle(block.setup.block_supname)

%% XCORR OF CELL ACTIVITY vs. Z STACK BEST POSITION

if multiplaneData
    F7 = block.F.(planeName) - block.setup.constant.neucoeff*block.Fneu.(planeName);
else
    F7 = block.F - block.setup.constant.neucoeff*block.Fneu;
end
        
F7_loco = F7(:,best_z_index); %Resampled to match loco

%Compute DF/F
mean_F = mean(F7_loco,2);% average green for each cell
F7_df_f = (F7_loco - mean_F)./mean_F;%(total-mean)/mean
%A = smooth(F7_df_f,10);

cmat_loco = [];
cmat_bestz = [];

for i = 1:size(F7_loco,1)
    [c1, lags] = xcorr(loco_speed', F7_df_f(i,:), maxlag,  'coeff');
    [c2, lags] = xcorr(best_z_loco_centered, F7_df_f(i,:), maxlag, 'coeff');
    cmat_loco = [cmat_loco; c1];
    cmat_bestz = [cmat_bestz; c2];
end

%%
figure

subplot(3,2,1)
plot(lags,cmat_loco)
ylabel('Corr coeff')
title('xcorr(Loco, F)')

subplot(3,2,2)
plot(lags,cmat_bestz)
ylabel('Corr coeff')
title('xcorr(Z difference, F)')

subplot(3,2,3)
plot(lags,mean(cmat_loco))
ylabel('Avg corr coeff')

subplot(3,2,4)
plot(lags,mean(cmat_bestz))
ylabel('Avg corr coeff')

subplot(3,2,5)
imagesc(cmat_loco)
caxis([0 0.6])
ylabel('Cells')
set(gca, 'XTick', 0:loco_sampling_rate:2*maxlag)
set(gca, 'XTickLabel', -maxlag/loco_sampling_rate:1:maxlag/loco_sampling_rate)
xlabel('Lag (s)')
    
subplot(3,2,6)
imagesc(cmat_bestz)
caxis([0 0.6])
ylabel('Cells')
set(gca, 'XTick', 0:loco_sampling_rate:2*maxlag)
set(gca, 'XTickLabel', -maxlag/loco_sampling_rate:1:maxlag/loco_sampling_rate)
xlabel('Lag (s)')

for i = 1:4
    subplot(3,2,i)
    set(gca, 'XTick', -maxlag:loco_sampling_rate:maxlag)
    set(gca, 'XTickLabel', -maxlag/loco_sampling_rate:1:maxlag/loco_sampling_rate)
    xlabel('Lag (s)')
end
suptitle(block.setup.block_supname)

figure
scatter(max(cmat_loco,[],2), max(cmat_bestz,[],2))
xlabel('xcorr(Loco, F)')
ylabel('xcorr(Z difference, F)')

suptitle(block.setup.block_supname)
end