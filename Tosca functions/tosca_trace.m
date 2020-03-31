function tosca_trace(fn, trial)

[d,p] = tosca_read_run(fn);
s = tosca_read_trial(p,d,trial);
tr = tosca_read_trace_data(strrep(fn, '.txt', '.trace.txt'));
tosca_plot_trial(s, tr(trial), 1 / p.DAQ.Poll_Rate_Hz);

[~,f] = fileparts(fn);
title(sprintf('%s-%d', f, trial));
