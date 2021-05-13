function tosca_trace(varargin)

if nargin == 2,
   fn = varargin{1};
   trial = varargin{2};
   [d,p] = tosca_read_run(fn);
else
   p = varargin{1};
   d = varargin{2};
   trial = varargin{3};
   fn = p.Info.Filename;
end
   
s = tosca_read_trial(p,d,trial);
tr = tosca_read_trace_data(strrep(fn, '.txt', '.trace.txt'));
tosca_plot_trial(s, tr(trial));

