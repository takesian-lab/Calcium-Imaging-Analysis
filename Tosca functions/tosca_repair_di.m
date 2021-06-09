function [srem, repaired] = tosca_repair_di(s, d, tr, trial, itrace, fn)

srem = [];

tEnd = Inf;
if itrace < length(tr)
   tEnd = tr(itrace+1).Time(1);
end
tStart = tr(itrace).Time(1);

tr = tr(itrace);
% d = d{trial};

if tStart>s.Time_s(1) || tEnd < s.Time_s(1)
   error('Unexpected DI reconstruction error.');
end

trimFrom = [];
repaired = 0;

t = tr.Time - s.Time_s(1);
di_time = s.Time_s - s.Time_s(1);

tRep = di_time(diff([0 s.Rep_Trigger]) > 0.5);

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

   % Look in a + 5-ms window around the nominal time
%    idi_time = abs(di_time - tRepMarker) < 5e-3;
%    idi_time = di_time>=tRepMarker & di_time<tRepMarker + 15e-3;
%    
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
      s.Loop_time_s(insert_at + 1) = s.Loop_time_s(1) + tRepMarker;
      s.Loop_time_s(insert_at + 2) = s.Loop_time_s(1) + tRepMarker + pulse_dur;
      
      % Set rep marker high
      s.Rep_Trigger(insert_at + (1:2)) = 0.5;
      
      % Is this also a state change? 
      tState = s.Time_s(diff(s.State_Change)>0.5);
      iStateEnqueued = find(tr.Event==2 & t<tRepMarker, 1, 'last');
      if ~any(tState>t(iStateEnqueued) & tState < tRepMarker)
         s.State_Change(insert_at + (1:2)) = 0.5;
      end
      
      % Is it a trial change?
      tReps = s.Time_s(diff(s.Rep_Trigger) < -0.5) - s.Time_s(1);
      if ~any(tReps < tRepMarker)
         s.Trial_Change(insert_at + (1:2)) = 0.5;
         trimFrom = insert_at + 1;
      end
      
      repaired = 1;
   end
end

if ~isempty(trimFrom)
   for kf = 1:length(s.DigitalNames)
      name = s.DigitalNames{kf};
      x = s.(name);
      s.(name) = x(trimFrom:end);
   end
end

% Is there a missing trial trigger? 'tEnd' (trace event == 1 for the
% following trial) will normally occur in the .di.txt file for the
% preceding trial. The indicator for a missing trial trigger is the
% presence of rep triggers beyond 'tEnd'.
tRep = s.Time_s(diff(s.Rep_Trigger)>0.5);
if any(tRep > tEnd)
   splitAt = find(s.Time_s >= tEnd, 1);
   srem = s;
   for kf = 1:length(s.DigitalNames)
      name = s.DigitalNames{kf};
      x = s.(name);
      srem.(name) = x(splitAt:end);
      s.(name) = x(1:splitAt-1);
   end
end

% While not an error per se, if the rep interval is set to zero as a means
% of bypassing the state, the state gets recorded in the history, but no
% corresponding markers are generated. 
% tStateEnqueued = tr.Time(tr.Event==6) - s.Time_s(1);
% frStateEnqueued = floor(tStateEnqueued / s.frameSize);
% ishort = find(diff(frStateEnqueued) < 2); % Two states enqueued in same frame
% if ~isempty(ishort)
%    tRepMarker = (frStateEnqueued(ishort(1)) + 1) * s.frameSize; 
%    [~, insert_at] = min(abs(di_time < tRepMarker));
%    [~, last_pt] = min(abs(di_time < tRepMarker + pulse_dur));
% 
%    % Set rep and state markers high
%    s.Rep_Trigger(insert_at:last_pt) = 0.5;
%    s.State_Change(insert_at:last_pt) = 0.5;
%    
%    repaired = repaired + 2;
% end


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
         fprintf(fp, '%d\t', s.(s.DigitalNames{kf})(kr));
      end
   end
   fprintf(fp, '\n');
end

fclose(fp);