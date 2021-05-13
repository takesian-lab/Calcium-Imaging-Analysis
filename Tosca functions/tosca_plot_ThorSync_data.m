function varargout = tosca_plot_ThorSync_data(folder, varargin)

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

figure;
pos = get(gcf, 'Position');
w = 1600;
pos(1) = max(0, pos(1) - (w - pos(3))/2);
pos(3) = w;
set(gcf, 'Position', pos);

hold on;
h = stem(tframe, ones(size(tframe)), 'r', 'Marker','none', 'Color', 0.7*[1 1 1]);
yd = get(h, 'YData');
set(h, 'YData', yd * 3);

h = stem(trep, ones(size(trep)), 'r', 'Color', [1 0.4 0], 'Marker','none', 'LineWidth', 2);
yd = get(h, 'YData');
set(h, 'YData', yd * 3);
set(h, 'BaseValue', 2);

h = stem(tstate, ones(size(tstate)), 'g', 'Marker','none', 'LineWidth', 2, 'Color', [0 0.7 0]);
yd = get(h, 'YData');
set(h, 'YData', yd * 2);
set(h, 'BaseValue', 1);

h = stem(ttrial, ones(size(ttrial)), 'b', 'Marker','none', 'LineWidth', 2);
set(h, 'BaseValue', 0);

% y = 0.5*s.speaker1/max(abs(s.speaker1));
% plot(s.time, y+1.5, 'k-');

set(gca, 'YTick', 0.5:1:2.5, 'YTickLabel', {'Trial','State','Rep'});
set(gca, 'TickDir', 'out');

xlabel('Time (s)');
