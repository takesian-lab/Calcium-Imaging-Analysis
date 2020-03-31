function figsize(a, SF)

if nargin == 0, a='thissize'; end
if nargin<2, SF=1; end

maximize = false;

if ischar(a)
   switch a,
      case 'portrait',
         Px = 8.5;
         Py = 11;
         orient portrait;
      case 'landscape',
         Px = 11;
         Py = 8.5;
         orient landscape;
      case 'maximize'
         screensize = get(0, 'ScreenSize');
         Px = screensize(3);
         Py = screensize(4);
         maximize = true;
         SF = 1;
      otherwise
         error('Invalid figure size specification string.');
  end
else
  Px = a(1);
  Py = a(2);
end

if length(a)==3,
   margin = a(3);
else
   margin = 0;
end

ou = get(gcf,'Units');
set(gcf,'Units','inches');
oldp = get(gcf,'Position');
top = oldp(2)+oldp(4);
left = oldp(1);
if strcmp(a,'thissize'),
   Px = oldp(3);
   Py = oldp(4);
end
if maximize || any(SF*[Px Py] > 10),
   set(gcf,'Position',[0 0 Px Py]);
else
   set(gcf,'Position',[left top-Py*SF Px*SF Py*SF]);
end
set(gcf,'Units',ou);

set(gcf,'PaperPosition',[margin margin Px-2*margin Py-2*margin]);
