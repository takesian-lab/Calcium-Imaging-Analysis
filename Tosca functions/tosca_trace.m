function TR = tosca_trace(varargin)
% TOSCA_TRACE
% Usage 1: tosca_trace(fn, trialNum)
% Usage 2: tosca_trace(p, d, trialNum);
%

if nargin == 2
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

if nargout
   TR = tr(trial);
   return;
end

tosca_plot_trial(s, tr(trial));
