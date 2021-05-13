function Data = tosca_reconstruct_di(P, Data)
% TOSCA_RECONSTRUCT_DI -- reparse the DI data to ensure correct alignment.
% 

% Backup data
disp('Tosca: backing up DI data.');
[folder, fn] = fileparts(P.Info.Filename);
if isempty(folder)
   folder = pwd;
end
backupFolder = fullfile(folder, 'backup');
if ~exist(backupFolder, 'dir')
   mkdir(backupFolder);
end

flist = dir(fullfile(folder, [fn, '-*.di.txt']));
for k = 1:length(flist)
    movefile(fullfile(folder, flist(k).name), fullfile(backupFolder, flist(k).name));
end

Trace = tosca_read_trace_data(strrep(P.Info.Filename, '.txt', '.trace.txt'));

P.Info.Filename = fullfile(backupFolder, fn);

idi_in = 1;
idi_out = 1;

idata = 1;
itrace = 1;
while true
   if idata > length(Data), break; end
   if length(Data{idata}.History)==1 && strcmp(Data{idata}.History{1}(1:5), 'ERROR')
      idata = idata + 1;
      continue;
   end
  
   % If no state enqueued (Event=5) --> duplicate error message in Tosca
   % resulting in trial with no AO and hence no digital data.
   if ~any(Trace(itrace).Event == 5)
      fprintf('Tosca: trial #%d is empty\n', idata);
      Data{idata}.N = NaN;
      Data{idata}.Result = 'empty';
%       idi_out = idi_out + 1;
      idata = idata + 1;
      itrace = itrace + 1;
      continue;
   end

   Data{idata}.N = idi_in;
   di_fn = fullfile(backupFolder, sprintf('%s-Trial%02d.di.txt', fn, idi_in));
   if ~exist(di_fn, 'file') && idata == length(Data)
      break;
   end
   
   s = tosca_read_trial(P, Data, idata);     

   tEnd = Inf;
   if itrace < length(Trace)
      tEnd = Trace(itrace+1).Time(1);
   end
   tRep = s.Time_s(diff(s.Rep_Trigger)>0.5);
   if s.Rep_Trigger(1) > 0.5, tRep = [s.Time_s(1) tRep]; end
%    missingTrialMarker = any(tRep > tEnd);
   missingTrialMarker = tRep(1) > tEnd;
   if missingTrialMarker
      fprintf('Tosca: trial #%d is empty\n', idata);
      Data{idata}.N = NaN;
      Data{idata}.Result = 'empty';
      idata = idata + 1;
      itrace = itrace + 1;
      continue;
   end
   
   
   % Does number of state changes match history length? If yes, copy the
   % file back with (potentially) new index.
   numStateChange = sum(diff(s.State_Change) < 0); % use falling edge so we count them all
   if numStateChange == length(s.History) / 2
      fn_in = sprintf('%s-Trial%02d.di.txt', fn, idi_in);
      fn_out = sprintf('%s-Trial%02d.di.txt', fn, idi_out);
      
      if ~strcmpi(fn_in, fn_out)
         fprintf('Tosca: %s --> %s\n', fn_in, fn_out);
      end
      
      copyfile(fullfile(backupFolder, fn_in), fullfile(folder, fn_out));
      
      Data{idata}.N = idi_out;

      idi_in = idi_in + 1;
      idi_out = idi_out + 1;
      idata = idata + 1;
      itrace = itrace + 1;
      continue;
   end
   
   while true
      fn_out = fullfile(folder, sprintf('%s-Trial%02d.di.txt', fn, idi_out));
      [s, repaired] = tosca_repair_di(s, Data, Trace, idata, itrace, fn_out);
   
      if bitand(repaired, 1)
         fprintf('Tosca: repaired trial #%d\n', idata);
         Data{idata}.Result = 'suspect di';
      elseif bitand(repaired, 2)
         fprintf('Tosca: marked null state in trial #%d\n', idata);
      end
      Data{idata}.N = idi_out;
      idi_out = idi_out + 1;
      idata = idata + 1;
      itrace = itrace + 1;
      
      if isempty(s), break; end
      
      disp('Tosca: SPLIT (missing trial marker)');
   end
   idi_in = idi_in + 1;
end