function tosca_plot_raw_data(fn)

text = fileread(fn);

inl = find(text==char(13), 1);
headerRow = text(1:inl-1);

a = textscan(strrep(headerRow,' ', '_'), '%s');

varNames = a{1};
Nvar = length(varNames);

data = textscan(text(inl+1:end), '%f');
M = reshape(data{1}, Nvar, length(data{1})/Nvar)';
Npts = size(M,1);

H = tosca_read_header(fn);
dt = 1 / H.DAQ.Poll_Rate_Hz;


figure;
figsize([10 4]);
% SetFigStyle('work');

hold on;

cmap = jet(Nvar);
for k = 1:Nvar,
   y = M(:,k);
   y(y>0) = k-1 + 0.95;
   y(y==0) = k-1;
   h = stairs((0:Npts-1) * dt, y);
   set(h, 'Color', cmap(k,:));
end

fnRawDI = strrep(fn, '.txt', '.rawDI');
if exist(fnRawDI, 'file'),
   [yraw, traw] = tosca_read_raw_di(fnRawDI);
   for k = 1:size(yraw,2),
      plot(traw, 0.95*yraw(:,k)+ Nvar + k - 1);
%       set(h, 'Color', cmap(k,:));
   end
   Nvar = Nvar + size(yraw, 2);
%    varNames = [varNames 
end

set(gca, 'XLim', [0 Npts] * dt);
set(gca, 'YLim', [-0.25 Nvar]);
set(gca, 'YTick', -0.5 + 1:Nvar);
set(gca, 'YTickLabel', varNames);
set(gca, 'Color', 'none');

   