function tosca_find_errors(P)

if nargin < 1,
   P = pwd;
end

flist = dir(fullfile(P, '*Run*.txt'));

for k = 1:length(flist),
   t = regexp(flist(k).name, '-Run([\d]+).txt', 'tokens');
   if ~isempty(t),
      disp(flist(k).name);
      run = str2double(t{1});
      [d, p] = tosca_read_run(fullfile(P, flist(k).name),false);
      for kd = 1:length(d),
         ierr = find(strncmp(d{kd}.History, 'Error', 5));
         if ~isempty(ierr),
            fprintf('Run %d-Trial%d: %s\n', run, d{kd}.trial, d{kd}.History{ierr});
         end
      end
      
   end
end