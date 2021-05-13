function xlimits=xaxis(xmin, xmax)
%XAXIS - Reset limits of X axis
%usage: xaxis(xmin, xmax)

% limits=axis;
% 
% if nargin == 0,
% 	xlimits=limits(1:2); 
% elseif nargin == 1,
%    if length(xmin) > 1,
%       limits(1:2) = [min(xmin) max(xmin)];
%    else
%       limits(1) = xmin;
%    end
%    axis(limits);
% else	
% 	limits(1:2) = [xmin xmax];
% 	axis(limits)
% end
set(gca, 'XLim', [xmin xmax]);