function AVI = tosca_create_avi_log(tl, folder)
% TOSCA_CREATE_AVI_LOG -- parse .avi.txt files into trial-based structure
% Usage: AVI = tosca_create_avi_log(tl, folder)
%
% --- Inputs ---
% tl : Tosca log output from tosca_create_log
% folder: folder containing .avi data
%

[~, fstem] = fileparts(tl.filename);

avilist = dir(fullfile(folder, [fstem '.*.avi.txt']));

frameInAVI = [];
fileNum = [];
tframe = [];

trTime = [];
trNum = [];
firstFrameInTrial = [];
overallFrameNum = [];
n0 = -1;

for k = 1:length(avilist)
   fp = fopen(fullfile(folder, avilist(k).name), 'rt');
   data = textscan(fp, '%d\t%d\t%f\t%f');
   fclose(fp);
   
   vidFile = fullfile(folder, strrep(avilist(k).name, '.txt', ''));
   vid = VideoReader(vidFile);

   n = data{2};
   if k == 1, n0 = n(1); end
   
   frameInAVI = [frameInAVI; (1:vid.NumFrames)'];
   fileNum = [fileNum; (k-1) * ones(vid.NumFrames, 1)];
   overallFrameNum = [overallFrameNum; n(n>0) - n0 + 1];
   
   t = data{3};

   npad = [n; max(n(n>0))+1];
   itr = find(n < 0);
   firstFrameInTrial = [firstFrameInTrial; npad(itr+1)-n0 + 1];
   
   % Tosca trials indicated by n < 0. Corresponding time stamps are on the
   % Tosca clock. Make the trial time the next frame immediately following.
%    trTime = [trTime; t(find(n < 0) + 1)];
%    
%    % Tosca trial number is the absolute value of 'n'
%    trNum = [trNum; abs(n(n<0))];
%    
%    % actual frame times
    ft = t(n > 0); 
%    
%    nfr = sum(n > 0);
%    if nfr ~= vid.NumFrames
%       % There is a bug where the .avi.txt file stops accumulating before
%       % the .avi file fills up. Fill in the missing data here
%       fprintf('%s: number of frames in .avi does not match log file.\n', avilist(k).name); 
%       
%       nadd = vid.NumFrames - nfr;
%       dt = mean(diff(ft)); % mean frame rate, computed from existing frames
%       ft = [ft; ft(end) + ((1:nadd)')*dt];
%    end
    tframe = [tframe; ft];
end

% overallFrameNum = cumsum(frameInAVI);
% 
% tosca_trial_times = NaN(length(tl.trials), 1);
% for k = 1:length(tl.trials)
%    tosca_trial_times(k) = tl.trials{k}.start;
% end
% 
% % Align first Tosca trial time with first .avi trial time
% tosca_trial_times = tosca_trial_times - tosca_trial_times(1) + trTime(1);

% Get the Tosca trial times in frame-referenced time
% if length(trNum) ~= length(tosca_trial_times)
%    trTime = interp1(tosca_trial_times(trNum), trTime, tosca_trial_times, 'linear', 'extrap');
% else
%    trTime = tosca_trial_times(trNum);
% end

% Create a structure array with one entry per Tosca trial
%trTime(end+1) = max(tframe);
firstFrameInTrial = [firstFrameInTrial; overallFrameNum(end)+1];
for k = 1:length(firstFrameInTrial)-1
%    AVI(k).start = trTime(k);
%    AVI(k).stop = trTime(k+1);
%    
%    idx = tframe >= trTime(k) & tframe < trTime(k+1);
   idx = overallFrameNum >= firstFrameInTrial(k) & overallFrameNum < firstFrameInTrial(k+1);
   
   AVI(k).frameInAVI = frameInAVI(idx);
   AVI(k).aviNum = fileNum(idx);
   AVI(k).frames = overallFrameNum(idx);
    AVI(k).tframe = tframe(idx);
end
