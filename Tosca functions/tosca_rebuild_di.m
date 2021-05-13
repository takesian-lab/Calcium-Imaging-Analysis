function tosca_rebuild_di

fn = 'D:\Data\Tosca\DFI024\Session 2\DFI024-Session2-Run2.txt';

tr = tosca_read_trace_data(strrep(fn, '.txt', '.trace.txt'));

[d,p] = tosca_read_run(fn, false);

s = [];

idi_in = 1;
idi_out = 1;

for k = 1:length(tr),
   fprintf('Trial %d\n', k);
   
   tStart = tr(k).Time(1);
   if k < length(tr),
      tEnd = tr(k+1).Time(1);
   else
      tEnd = Inf;
   end
   
   while (isempty(s) || tEnd > s.Time_s(end)) && idi_in <=125,
      fprintf('Reading %d\n', idi_in);
      s = add_data(s, d, p, idi_in);
      idi_in = idi_in + 1;
   end

   istart = find(s.Time_s > tStart & s.State_Change > 0, 1);
   istop = find(s.Time_s < tEnd, 1, 'last');
   
   ikeep = istart:istop;
   s_save = s;
   for kf = 1:length(s.DigitalNames),
      name = s.DigitalNames{kf};
      x = s.(name);
      s_save.(name) = x(ikeep);
   end
   
   istateoff = find(diff(s_save.State_Change) < 0, 1);
   s_save.Trial_Change(1:istateoff) = 1;
   
   [~, f] = fileparts(fn);
   resave(sprintf('%s-Trial%02d.di.txt', f, idi_out), s_save);
   idi_out = idi_out + 1;
   
   for kf = 1:length(s.DigitalNames),
      name = s.DigitalNames{kf};
      x = s.(name);
      s.(name) = x(ikeep(end)+1:end);
   end

   
end

%--------------------------------------------------------------------------
function s = add_data(s, d, p, ifile)

s_new = tosca_read_trial(p, d, ifile);

if isempty(s),
   s = s_new;
else
   for kf = 1:length(s.DigitalNames),
      name = s.DigitalNames{kf};
      s.(name) = [s.(name) s_new.(name)];
   end
end

%--------------------------------------------------------------------------
function resave(fn, s)

fp = fopen(fn, 'wt');

for kf = 1:length(s.DigitalNames),
   fprintf(fp, '%s\t', s.DigitalNames{kf});
end
fprintf(fp, '\n');

for kr = 1:length(s.Time_s),
   for kf = 1:length(s.DigitalNames),
      if kf < 3,
         fprintf(fp, '%.6f\t', s.(s.DigitalNames{kf})(kr));
      else
         fprintf(fp, '%d\t', s.(s.DigitalNames{kf})(kr));
      end
   end
   fprintf(fp, '\n');
end

fclose(fp);