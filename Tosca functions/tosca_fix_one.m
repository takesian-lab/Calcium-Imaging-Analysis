function tosca_fix_one(P, Data, trial)
% TOSCA_RECONSTRUCT_DI -- reparse the DI data to ensure correct alignment.
% 

Trace = tosca_read_trace_data(strrep(P.Info.Filename, '.txt', '.trace.txt'));

s = tosca_read_trial(P, Data, trial);

[folder, fn] = fileparts(P.Info.Filename);
fn_out = fullfile(folder, sprintf('%s-Trial%02d.di.txt', fn, trial));
movefile(fn_out, [fn_out, '.orig']);
[s, repaired] = tosca_repair_di(s, Data, Trace, trial, fn_out);

if bitand(repaired, 1),
   fprintf('Tosca: repaired trial #%d\n', trial);
   Data{trial}.Result = 'suspect di';
elseif bitand(repaired, 2),
   fprintf('Tosca: marked null state in trial #%d\n', trial);
end
