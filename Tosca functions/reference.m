function H = reference(whichAxis, refVal, lineSpec, varargin)
%REFERENCE -- add reference lines to a plot.
%Usage: h = reference(whichAxis, refVal, lineSpec)
%
% Inputs:
% whichAxis : 'x' or 'y'
% refVal    : vector of reference values
% lineSpec  : default = 'k-'
%
% Output:
% h         : handle of reference line objects
%
% Should be called *after* the plot limits have been set.
%

if nargin < 3, lineSpec='k-'; end

xl = get(gca, 'XLim');
if isequal(get(gca,'XScale'), 'log'),
   xl(1) = min(get(gca, 'XTick'));
end
yl = get(gca, 'YLim');

hold on;

switch lower(whichAxis),
case 'x',
   for k = 1:length(refVal),
      h(k) = plot(refVal(k)*[1 1], yl, lineSpec, varargin{:});%, 'LineWidth', 1.5);
   end
   set(gca, 'YLim', yl);
   
case 'y',
   for k = 1:length(refVal),
      h(k) = plot(xl, refVal(k)*[1 1], lineSpec, varargin{:});
   end
   
case 'x==y',
   v = [max([xl yl]) min([xl yl])];
   h = plot(v, v, lineSpec);

case 'y=x+',
   v = [max([xl yl]) min([xl yl])];
   h = plot(v, v+refVal(1), lineSpec);

case 'pt',
   h = plot([xl(1) refVal(1) refVal(1)], [refVal(2) refVal(2) yl(1)],...
      lineSpec, 'LineWidth', 1.5, varargin{:});
   
   
otherwise,
   error('reference: invalid parameter.');
   
end

for k = 1:length(h),
   set(h(k), 'Tag', ['reference,' whichAxis]);
   hg = get(h(k), 'Annotation');
   hl = get(hg, 'LegendInformation');
   set(hl, 'IconDisplayStyle', 'off');
end

if nargout, H=h; end