function repaired = tosca_restore_di_markers(fn, s, tr)

repaired = false;

if ~any(s.Trial_Change > 0)
   return;
end

t = tr.Time - s.Time_s(1);
di_time = s.Time_s - s.Time_s(1);

tRep = di_time(diff([0 s.Rep_Trigger]) > 0.5);
if any(tRep > t(end))
   return;
end

idxWrite1 = find(tr.Event == 5 & tr.Data == 1);

offset = 0;
if t(idxWrite1(1)) < 0
   offset = 1;
end

for k = 1:length(idxWrite1)
   % Writing the 1st buffer to the hardware occurs while the 0th buffer is
   % playing. Thus, the nearest frame boundary before the 1st buffer write
   % should correspond to a rep marker.
   tRepMarker = (floor(t(idxWrite1(k)) / s.frameSize) + offset) * s.frameSize; 

   % Look in a + 15-ms window around the nominal time
%    idi_time = di_time>=tRepMarker & di_time<tRepMarker + 15e-3;
   
%    if ~any(idi_time)
   if ~any(tRep >= tRepMarker - 5e-3 & tRep < tRepMarker + 20e-3)
      % No rep marker found. Insert it.
      
      insert_at = find(di_time < tRepMarker, 1, 'last');
      pulse_dur = min(20e-3, di_time(insert_at+1)-tRepMarker);
      
      di_time = [di_time(1:insert_at) tRepMarker tRepMarker+pulse_dur di_time(insert_at+1:end)];

      % Need to expand all of the di fields
      for kf = 1:length(s.DigitalNames)
         x = s.(s.DigitalNames{kf});
         if insert_at < length(x)
            s.(s.DigitalNames{kf}) = x([1:insert_at insert_at insert_at+1 insert_at+1:length(x)]);
         else
            s.(s.DigitalNames{kf}) = x([1:insert_at insert_at insert_at]);
         end
      end

      % Set inserted time values to edges of 20-ms pulse
      s.Time_s(insert_at + 1) = s.Time_s(1) + tRepMarker;
      s.Time_s(insert_at + 2) = s.Time_s(1) + tRepMarker + pulse_dur;
      
      % Set rep marker high
      s.Rep_Trigger(insert_at + (1:2)) = 0.5;
      
      % Is this also a state change? 
      tState = s.Time_s(diff(s.State_Change)>0.5);
      iStateEnqueued = find(tr.Event==2 & t<tRepMarker, 1, 'last');
      if ~any(tState>t(iStateEnqueued) & tState < tRepMarker)
         s.State_Change(insert_at + (1:2)) = 0.5;
      end
      
      repaired = true;
   end
end

folder = fileparts(fn);
if isempty(folder)
   folder = pwd;
end
backupFolder = fullfile(folder, 'backup');
if ~exist(backupFolder, 'dir')
   mkdir(backupFolder);
end
copyfile(fn, backupFolder);

resave(fn, s);

%--------------------------------------------------------------------------
function resave(fn, s)

fp = fopen(fn, 'wt');

for kf = 1:length(s.DigitalNames)
   fprintf(fp, '%s\t', s.DigitalNames{kf});
end
fprintf(fp, '\n');

for kr = 1:length(s.Time_s)
   for kf = 1:length(s.DigitalNames)
      if kf < 3
         fprintf(fp, '%.6f\t', s.(s.DigitalNames{kf})(kr));
      else
         fprintf(fp, '%g\t', s.(s.DigitalNames{kf})(kr));
      end
   end
   fprintf(fp, '\n');
end

fclose(fp);