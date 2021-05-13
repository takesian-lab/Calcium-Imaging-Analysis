function TL = tosca_create_log(FN, varargin)
% TOSCA_CREATE_LOG -- add timing details to Tosca trial data.
% Usage: TL = tosca_create_log(FN)
% Usage: TL = tosca_create_log(FN, 'aviFolder', aviFolder)
%
% --- Inputs ---
% FN : main Tosca data file
% aviFolder : if specified, adds video frame information
%

aviFolder = '';
ParseArgin(varargin{:});

[D, P] = tosca_read_run(FN);

TL.filename = P.Info.Filename;
TL.version = P.Info.Version;
TL.params = P;

for k = 1:length(D)
   TR{k} = D{k};
   
   s = tosca_read_trial(P, D, D{k}.N);

   % trial timing information
   TR{k}.start = s.Time_s(1);
   TR{k}.stop = max(s.Time_s);
   TR{k}.duration = TR{k}.stop - TR{k}.start;

   % individual state timing information
   tState = [s.Time_s(1) s.Time_s(find(diff(s.State_Change)>0) + 1) max(s.Time_s)];
   names = s.History(1:2:end);
   if length(names) ~= length(tState)-1
      names = names(1:length(tState)-1);
%       warning('mismatch between state names and times'); 
   end
   
   for ks = 1:length(names)
      stateData = struct('name', strtrim(names{ks}), ...
         'start', tState(ks), 'stop', tState(ks+1), 'duration', tState(ks+1) - tState(ks));
      TR{k}.states(ks) = stateData;
   end
   
end

TL.trials = TR;

% merge .avi data, if requested
if ~isempty(aviFolder)
   aviLog = tosca_create_avi_log(TL, aviFolder);
   TL = tosca_merge_avi_log(TL, aviLog);
   TL.aviFolder = aviFolder;
end

TL = tosca_create_loco_log(TL);

% save to .mat file
save(strrep(TL.filename, '.txt', '.log.mat'), 'TL');

