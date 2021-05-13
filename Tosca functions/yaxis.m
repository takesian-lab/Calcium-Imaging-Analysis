function ylimits = yaxis(ymin, ymax)
%YAXIS - Reset limits of Y axis
%usage: yaxis(ymin, ymax)

% limits = axis;
% if nargin == 0,
% 	ylimits = limits(3:4);
% elseif nargin == 1,
%    limits(3) = ymin;
%    axis(limits);
% else
% 	limits(3:4) = [ymin ymax];
% 	axis(limits)
% end
% 
set(gca, 'YLim', [ymin ymax]);