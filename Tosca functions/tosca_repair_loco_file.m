function data = tosca_repair_loco_file(fileName, TL)
% TOSCA_REPAIR_LOCO_FILE -- inserts missing trial markers into locomotion
% data.
%
% Usage: data = tosca_repair_loco_file(fileName, TL)
%

fp = fopen(fileName, 'rt');
C = textscan(fp, '%f%f%f');
fclose(fp);

t = C{1};
ch = C{2};
speed = C{3};

ttr = t(ch==1);

for k = 1:length(TL.trials)
   if min(abs(TL.trials{k}.start - ttr)) > 0.5
      insertAt = find(t < TL.trials{k}.start, 1, 'last');
      t = [t(1:insertAt); TL.trials{k}.start; t(insertAt+1:end)];
      ch = [ch(1:insertAt); 1; ch(insertAt+1:end)];
      speed = [speed(1:insertAt); NaN; speed(insertAt+1:end)];
      
      fprintf('Inserted marker for trial %d\n', k);
   end
end

movefile(fileName, strrep(fileName, '.txt', '.orig.txt'));

fp = fopen(fileName, 'wt');
for k = 1:length(t)
   fprintf(fp, '%.6f\t%.6f\t%.6f\n', t(k), ch(k), speed(k));
end
fclose(fp);

if nargout
   data = tosca_read_loco(fileName);
end

