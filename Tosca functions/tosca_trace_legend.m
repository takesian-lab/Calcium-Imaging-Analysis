% function tosca_trace_legend()

if ishandle(555),
   return;
end

figure(555);
set(gcf, 'NumberTitle', 'off', 'Name', 'Trace Legend');

x = [0 1];
y = [0 NaN];

hold on;
plot(x, y, 'g-', 'LineWidth', 2);
plot(x, y, 'k-', 'LineWidth', 2);
plot(x, y, 'b:', 'LineWidth', 2);
plot(x, y, 'r-', 'LineWidth', 2);
plot(x, y, 'c-', 'LineWidth', 2);
plot(x, y, 'b-', 'LineWidth', 2);
plot(x, y, 'm:', 'LineWidth', 2);
plot(x, y, 'm-', 'LineWidth', 2, 'Color', [0.85 0.6 0]);

hl = legend( ...
   'Unspecified', ...
   'Trial start', ...
   'State enqueued', ...
   'Result', ...
   'Parser end received', ...
   'State dequeued', ...
   'Timeout sent', ...
   'TTL' ...
   );
   
axis off;
legend boxoff;

set_axis_size([2 2]);
p = get(gca, 'Position');
set(hl, 'Position', p);

tight_figure;
