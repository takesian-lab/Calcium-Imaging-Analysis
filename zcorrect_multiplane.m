function zcorrect_multiplane(block, plotFigureOnly)
% zcorrect_multiplane(block, plotFigureOnly)
%
% DOCUMENTATION IN PROGRESS
%
% Compute zcorr for multiplane data and make Z-corrected blocks
% 
% Argument(s): 
%   block (struct)
%   plotFigureOnly (optional 0 or 1)
% 
% Returns:
%   
% 
% Notes:
%   Run split_BOTs_by_plane to separate TIFFs into different folders based
%   on plane number
%
% TODO: 
% Search 'TODO'

%% Setup

if nargin < 2
    plotFigureOnly = 1;
end

max_drift = 10; %Number of planes that the imaging plane can drift in +/- Z

%locomotor activity
timestamp = block.timestamp;
loco_time = block.locomotion_trace;
loco_speed = block.loco_activity;
        
if length(loco_speed) ~= length(loco_time)
    shorter_loco = min(length(loco_speed),length(loco_time));
    loco_time = loco_time(1:shorter_loco);
    loco_speed = loco_speed(1:shorter_loco);
    warning('locomotion_trace and loco_activity are not the same length')
end  

%% Extract zcorr from block and establish best imaging plane
zcorr = block.ops.zcorr;
planes = fieldnames(zcorr);
nPlanes = numel(planes);

%Standardize Z matrix size
min_frames = size(zcorr.(planes{1}), 2);
for p = 2:nPlanes
    if size(zcorr.(planes{p}), 2) < min_frames
        min_frames = size(zcorr.(planes{p}), 2);
    end
end

%Concatenate
zcorr_concat = [];
concat_borders = ones(1,nPlanes + 1);
for p = 1:nPlanes
    if p > 1
        concat_borders(p) = size(zcorr_concat,2) + 1;
    end
    zcorr_concat = [zcorr_concat,zcorr.(planes{p})(:,1:min_frames)];
end
concat_borders(end) = size(zcorr_concat,2) + 1;

figure; hold on
    
for p = 1:nPlanes
    B1 = concat_borders(p);
    B2 = concat_borders(p+1) - 1;
    
    [~, best_z_raw] = max(zcorr_concat,[],1);
    imaging_plane = mode(best_z_raw);
    A = imaging_plane - max_drift;
    B = imaging_plane + max_drift;
    if A < 1
        A = 1;
    end
    zcorr_cropped = zcorr_concat(A:B,:);
    [best_z_val, best_z] = max(zcorr_cropped,[],1);
    imaging_plane_cropped = mode(best_z);
 
    % PLOT

    subplot(4,nPlanes,p)
    imagesc(zcorr_concat(:,B1:B2))
    hline(A, 'w')
    hline(B, 'w')
    ylabel('Z position')
    title(['Plane' num2str(p-1)])
    caxis([min(min(zcorr_concat)), max(max(zcorr_concat))])

    subplot(4,nPlanes,p+nPlanes)
    imagesc(zcorr_concat(A:B,B1:B2))
    hline(imaging_plane_cropped, 'r')
    set(gca, 'YTick', 1:5:(B-A)+1)
    set(gca, 'YTickLabel', A:5:B+1)
    ylabel('Cropped Z position')
    caxis([min(min(zcorr_concat)), max(max(zcorr_concat))])

    subplot(4,nPlanes,p+nPlanes*2)
    plot(best_z(B1:B2))
    hline(imaging_plane_cropped, 'r')
    hline(mode(best_z(B1:B2)), 'k')
    set(gca, 'YTick', 1:5:(B-A)+1)
    set(gca, 'YTickLabel', A:5:B+1)
    ylabel('Best Z stack position')
    xlim([1, length(best_z(B1:B2))])
    ylim([min(best_z) max(best_z)])
    xlabel('Frames')

    subplot(4,nPlanes,p+nPlanes*3)
    plot(loco_time, loco_speed, 'r')
    xlim([timestamp(1) timestamp(end)]) %loco will be shorter because tosca ends before PV
    xlabel('Time (s)')
    ylabel('Loco activity (cm/s)')
end

suptitle(block.setup.block_supname)

%% Determine best plane for Z-corrected block

%Restack
best_z_stack = [];
best_z_val_stack = [];
for p = 1:nPlanes
    B1 = concat_borders(p);
    B2 = concat_borders(p+1) - 1;
    best_z_stack = [best_z_stack; best_z(B1:B2)];
    best_z_val_stack = [best_z_val_stack; best_z_val(B1:B2)];
end
best_z_stack_relative = abs(best_z_stack - imaging_plane_cropped);

%Choose best plane based on:
% 1. Frame is most correlated with imaging plane
% 2. (If multiple options after #1) Frame closest to the most accurate plane
% 3. (If mutliple options after #1 & #2) Frame with highest z-corr value

equal_to_imaging_plane = zeros(size(best_z_stack_relative));
equal_to_imaging_plane(best_z_stack_relative == 0) = 1;
[~, most_accurate_plane] = max(sum(equal_to_imaging_plane,2));

best_plane = nan(1,min_frames);

for f = 1:min_frames
    [I, ~] = find(best_z_stack_relative(:,f) == min(best_z_stack_relative(:,f)));
    
    if length(I) == 1
        best_plane(f) = I;  %Criteria #1
    else
        [II, ~] = find(abs(I-most_accurate_plane) == min(abs(I-most_accurate_plane)));            
        if length(II) == 1
            best_plane(f) = I(II);  %Criteria #2
        else
            [~, I_max] = max(best_z_val_stack(I,f));
            best_plane(f) = I(I_max);  %Criteria #3
        end
    end
    
    if isnan(best_plane(f))
        error('Could not assign best plane')
    end
end

% PLOT

figure; hold on
subplot(4,1,1); hold on
imagesc(best_z_stack_relative)
colormap('jet')
ylabel('Plane')
xlabel('Frames')
title('Frames from best z-position')
set(gca,'YTick',1:nPlanes)
set(gca,'YTickLabel',0:(nPlanes-1))
xlim([0.5 length(best_plane) + .5])
ylim([0.5 nPlanes + .5])

subplot(4,1,2)
plot(best_plane)
xlim([0 length(best_plane)])
ylim([0 nPlanes])
title('Best plane')
set(gca,'YTick',1:nPlanes)

subplot(4,1,3); hold on
imagesc(best_z_stack_relative)
colormap('jet')
ylabel('Plane')
xlabel('Frames')
title('Overlay')
set(gca,'YTick',1:nPlanes)
set(gca,'YTickLabel',0:(nPlanes-1))
plot(best_plane, 'w')
xlim([0.5 length(best_plane) + .5])
ylim([0.5 nPlanes + .5])

subplot(4,1,4)
hist(best_plane, 0.5:1:(nPlanes+0.5))
set(gca,'XTick', 0.5:1:(nPlanes+0.5))
set(gca,'XTickLabel', 0:nPlanes)
title('N best frames per plane')
xlabel('Plane')
ylabel('Count')

suptitle(block.setup.block_supname)

%%

if plotFigureOnly
    return
end

%% Make Z-corrected block

disp('Performing Z-correction')

cd(block.setup.block_path)
[path,folder_name,~] = fileparts(pwd);
mkdir(path, ['Zcorrected-'  folder_name]);

filedir = dir;
filenames = {filedir(:).name};
tiffs_str = string(filenames(contains(filenames,'.ome.tif')))';
tiffs = char(tiffs_str);

%The filename format should be as follows:
%BOT_FILENAME-###_Cycle#####_Ch2_######.ome.tif
% ### is the block number
% ##### is the cycle number (used for t-series)
% ###### is the tif/frame number

%Filepath access
planeA = 13; %number of chars to start of frame number
planeB = 8;  %number of chars to end of frame number
cycleA = 23;
cycleB = 19;

filelist_lengths = strlength(deblank(tiffs_str)); %Get length of each filename
L = unique(filelist_lengths); %How many unique lengths are there
if length(L) > 1; error('Check filenames'); end

planes = double(string(tiffs(:,L-planeA:L-planeB)));
cycles = double(string(tiffs(:,L-cycleA:L-cycleB)));
unique_cycles = unique(cycles);
if length(unique_cycles) ~= length(best_plane)
    error('Number of cycles does not match Z computation')
end

moved_files = cell(length(best_plane),2);
new_timestamp = nan(length(best_plane),1);

for c = 1:length(unique_cycles)
    %select file
    current_cycle_files = tiffs_str(cycles == unique_cycles(c));
    best_plane_for_cycle = current_cycle_files(best_plane(c));
    moved_files{c,1} = best_plane_for_cycle;
    
    %Update timestamp
    file_ind = find(strcmp(tiffs_str,best_plane_for_cycle));
    new_timestamp(c) = block.timestamp(file_ind);

    %rename
    new_filename = char(best_plane_for_cycle);
    new_filename(1,L-planeA:L-planeB) = num2str(c,'%06.f');
    new_filename(1,L-cycleA:L-cycleB) = num2str(1,'%05.f'); 
    moved_files{c,2} = new_filename;
    
    %copy and rename file
    copyfile(best_plane_for_cycle, [path '/' 'Zcorrected-'  folder_name '/' new_filename])
end

%% Save new block.csv file with timestamp
cd([path '/' 'Zcorrected-'  folder_name])
block_csv_name = [new_filename(1:L-cycleB) '-botData.csv'];
textHeader = 'Timestamp';

%write header to file
fid = fopen(block_csv_name,'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);

%write data to end of file
dlmwrite(block_csv_name, new_timestamp ,'-append', 'precision', '%.7f');

%% Copy voltage recording csv to new block
copyfile([path '/' folder_name '/' block.setup.VR_filename], [path '/' 'Zcorrected-'  folder_name '/' block.setup.VR_filename])

%% Copy XML file to new block
copyfile([path '/' folder_name '/' block.setup.XML.filename], [path '/' 'Zcorrected-'  folder_name '/' block.setup.XML.filename])

%% Save old and new filenames to new block
moved_files_save = [{'Old filename', 'New filename'}; string(moved_files)];
xlswrite('BOT_file_identities', moved_files_save);

%%

disp('Done processing')

end