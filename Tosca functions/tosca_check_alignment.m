function Data = tosca_check_alignment(P, Data)
% TOSCA_CHECK_ALIGNMENT -- examines the data for a run and reports any
% instance in which the history length does not match the number of state
% changes.
%

traceFile = strrep(P.Info.Filename, '.txt', '.trace.txt');
if exist(traceFile, 'file')
   tr = tosca_read_trace_data(strrep(P.Info.Filename, '.txt', '.trace.txt'));
else
   tr = [];
end

idi = 1;

ok = true;
cause = '';

for k = 1:length(Data)
   if length(Data{k}.History)==1 && strcmp(Data{k}.History{1}(1:5), 'ERROR')
      Data{k}.N = -1;
      fprintf('Error (no trial): %d\n', k);
   elseif length(Data{k}.History)==1 || (~isempty(tr) && ~any(tr(idi).Event == 5))
      % If no state enqueued (Event=5) --> duplicate error message in Tosca
      % resulting in trial with no AO and hence no digital data.
      Data{k}.N = NaN;
      if ~isempty(tr) && any(tr(idi).Event == 6)
         warning('Weird case! Double check with Ken that everything is OK.');
         tr = tr([1:idi-1 idi+1:end]);
      end
      
   else
      Data{k}.N = idi;
      idi = idi + 1;

      [s, di_fn] = tosca_read_trial(P, Data, k);
      if isempty(s), continue; end

      % Check time stamps. The trial markers in the trace file (Event=1)
      % are sent by Tosca when it enters the flowchart block and are 
      % unambiguously correct.

      if ~isempty(tr)
         if idi < length(tr)
            tnext = tr(idi+1).Time(1);
         else
            tnext = Inf;
         end
         
         % Next trial should begin within the current .di.txt file
         if tnext < s.Time_s(1)
            ok = false;
            cause = sprintf('Trace/di mismatch: trial %d', k);
            break;
         end
      end   
      
      % Does number of state changes match history length?
      numStateChange = sum(diff(s.State_Change) < 0); % use falling edge so we count them all
      if numStateChange ~= length(s.History) / 2 && ~any(contains(s.History, 'Error'))
         ok = tosca_restore_di_markers(di_fn, s, tr(idi-1));
         if ok
            fprintf('Restored missing markers: trial %d\n', k);
         else
            cause = sprintf('State change/history mismatch: trial %d', k);
            break;
         end
      end

   end
end

if ~ok
   disp(cause);
   if ~isempty(tr)
      disp('Tosca: data alignment issue, reparsing DI.');
      Data = tosca_reconstruct_di(P, Data);
   else
     disp('Tosca: data alignment issue, no trace file.');
     for k = 1:length(Data)
        Data{k}.N = k;
     end
%       error('Tosca: data alignment issue, no trace file. Talk to Ken. Bring cash.');
   end
end

% Rewrite Run .txt file
movefile(P.Info.Filename, strrep(P.Info.Filename, '.txt', '.orig.txt'));

fpOld = fopen(strrep(P.Info.Filename, '.txt', '.orig.txt'), 'rt');
fpNew = fopen(P.Info.Filename, 'wt');

nd = 0;

while ~feof(fpOld)
   line = fgets(fpOld);
   
   if length(line)>6 && strcmp(line(1:6), 'Result')
      fprintf(fpNew, 'Result=%s\n', Data{nd}.Result);
   else
      fprintf(fpNew, '%s', line);
   end
   
   if length(line)>7 && strcmp(line(1:7), 'Version')
      fprintf(fpNew, 'AlignmentChecked=TRUE\n');
   end
   
   if length(line)>6 && (strcmp(line(1:6), '[Block') || strcmp(line(1:6), '[Track'))
      nd = nd + 1;
      fprintf(fpNew, 'N=%d\n', Data{nd}.N);
   end
   
end
fclose(fpOld);
fclose(fpNew);
