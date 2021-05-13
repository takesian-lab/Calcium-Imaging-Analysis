function varargout = tosca_read_ThorSync_data(folder, varargin)

fn = fullfile(folder, 'tsd.mat');

if exist(fn, 'file'),
   load(fn);
else
   fn = fullfile(folder, 'Episode001.h5');
   s = LoadThorSyncData(fn, ...
      'vars', {'time', 'Frame_Out', 'Tosca_Rep', 'Tosca_State', 'Tosca_Trial'}, ...
      varargin{:});

   tframe = s.time(diff(s.Frame_Out) > 0);
   trep = s.time(diff(s.Tosca_Rep) > 0);
   tstate = s.time(diff(s.Tosca_State) > 0);
   ttrial = s.time(diff(s.Tosca_Trial) > 0);

   save(fullfile(folder, 'tsd.mat'), 'tframe', 'ttrial', 'trep', 'tstate');
end

if nargout > 0,
   varargout{1} = tframe;
end
if nargout > 1,
   varargout{2} = ttrial;
end
if nargout > 2,
   varargout{3} = tstate;
end
if nargout > 3,
   varargout{4} = trep;
end
