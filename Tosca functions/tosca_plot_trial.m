function tosca_plot_trial(S, TR, varargin)
% TOSCA_PLOT_TRIAL -- plots detailed data for single Tosca trial.
% Usage: tosca_plot_trial(S)
% 
% Input:
%   S : structure returned by TOSCA_READ_TRIAL
%

if nargin < 2
   TR = [];
end
frameSize = S.frameSize;
Loco = [];
LocoMax = 20;
LocoThreshold = 3;
AI = [];
VBL = [];
labelFrames = true;

ParseArgin(varargin{:});

names = S.DigitalNames(3:end);
% t = S.Loop_time_s;
t = S.Time_s;
t = t - t(1);
y = NaN(length(t), length(names));

tStateChange = [t(1) t(find(diff(S.State_Change)>0) + 1)];
for k = 1:length(names)
   y(:,k) = 0.9*S.(names{k})+(k-1);
end

numRows = length(names);
if ~isempty(Loco), numRows = numRows + 1; end
if ~isempty(VBL), numRows = numRows + 1; end
if ~isempty(AI), numRows = numRows + length(AI.names); end

figure;
figsize([15 4]);
cmap = get(gca, 'ColorOrder');

if ~isempty(TR),
   tmin = 0;
   t0 = tmin:frameSize:max(t);
   
   frameColor = 0.85*[1 1 1];
   
   xx = [0 0 1 1 0];
   yy = [0 1 1 0 0];
   for k = 2:2:length(t0),
      h = patch(frameSize*xx + t0(k), yy*(length(names)+ 1), 'c');
      set(h, 'EdgeColor', frameColor, 'FaceColor', frameColor);
   end
   hold on;
   
   tRep = [-Inf t(1) t(find(diff(S.Rep_Trigger)>0) + 1) t(end) Inf];
   krep = 1;
   kfr = 0;
   for k = 1:length(t0),
      if t0(k) >= tRep(krep+1)-0.005 && tRep(krep+1)>=0,
         kfr = 0;
         krep = krep + 1;
      end
      
      if labelFrames,
         h = text(t0(k)+frameSize/2, numRows+0.25, num2str(kfr));
         set(h, 'HorizontalAlignment', 'center', 'FontSize', 8);
      end
      kfr = kfr + 1;
   end
end

h = stairs(t, y, 'LineWidth', 2);
set(h(1), 'LineStyle', 'none', 'Marker', '.');
hold on;
tmax = max(t);
icol = mod(length(names), size(cmap, 1)) + 1;

if ~isempty(Loco)
   names{end+1} = 'Speed';  
   tspeed = Loco.t - S.Time_s(1);
   ifilt = tspeed >=0 & tspeed<=max(t);
   tspeed = tspeed(ifilt);
   speed = Loco.speed(ifilt) / (2*LocoMax) + 0.5 + length(names)-1;
   
   plot(tspeed, speed, 'k', 'Color', cmap(icol, :));
   reference('y', length(names)-0.5);
   reference('y', length(names)-0.5 + LocoThreshold/(2*LocoMax), 'k:');
   
   icol = mod(icol, size(cmap, 1)) + 1;
end

if ~isempty(VBL)
   names{end+1} = 'VBL';

   [t, isort] = sort([VBL.Ton; VBL.Toff]);
   y = [ones(size(VBL.Ton)); zeros(size(VBL.Toff))];
   y = y(isort);
   ifilt = t>=S.Time_s(1) & t<=S.Time_s(end);
   t = [S.Time_s(1); t(ifilt); S.Time_s(end)] - S.Time_s(1);
   y = [0; y(ifilt); 0];
   
   y = 0.9*y + length(names)-1;

   stairs(t, y, 'b', 'LineWidth', 2, 'Color', cmap(icol, :));
   icol = mod(icol, size(cmap, 1)) + 1;
end

if ~isempty(AI)
   for k = 1:length(AI.names)
      names{end+1} = ['AI: ' AI.names{k}];
      
      y = AI.data(:,k) / (2*max(abs(AI.data(:,k)))) + 0.5 + length(names)-1;
      
      plot(AI.t, y, 'k', 'Color', cmap(icol,:));
      icol = mod(icol, size(cmap, 1)) + 1;
      reference('y', length(names)-0.5);
   end
end

figsize([15 0.5*length(names)]);
yaxis(0, numRows);

xaxis(-0.025*tmax, 1.025*tmax);
set(gca, 'YTick', 0.5 + (0:length(names)-1));
set(gca, 'YTickLabel', names);

xlabel('Time (s)');
reference('x', tStateChange, 'k:');

if ~isempty(TR),
   t = double(TR.Time - S.Time_s(1));

   % 0: Unspecified
   if any(TR.Event == 0),
      val = TR.Data(TR.Event==0);
      h = reference('x', t(TR.Event == 0), 'g-', 'LineWidth', 2);
      set(h(val==1), 'LineStyle', ':');
   end
   % 1: Trial start
   reference('x', t(TR.Event == 1), 'k-', 'LineWidth', 2);
   % 2: State enqueued
   reference('x', t(TR.Event == 2), 'b:', 'LineWidth', 2);
   % 3: Result
   reference('x', t(TR.Event == 3), 'r-', 'LineWidth', 1);
   % 4: Parser end received
   if any(TR.Event == 4),
      reference('x', t(TR.Event == 4), 'c-', 'LineWidth', 2);
   end
   
   % 5: Audio frame sent
   if labelFrames,
      ifr = TR.Event == 5;
      tfr = t(ifr);
      fr_num = TR.Data(ifr);
      for k = 1:length(tfr),
         h = text(tfr(k), numRows+0.75, num2str(fr_num(k)));
         set(h, 'HorizontalAlignment', 'center');
      end
   end
   
   % 6: State queue received by AO thread
   if any(TR.Event == 6),
      reference('x', t(TR.Event == 6), 'b-', 'LineWidth', 2);
   end
   % 7: Time out sent by AO thread
   if any(TR.Event == 7),
      reference('x', t(TR.Event == 7), 'm:', 'LineWidth', 2);
   end
   % 8: TTL change
   if any(TR.Event == 9),
      reference('x', t(TR.Event == 9), 'm-', 'LineWidth', 2, 'Color', [0.85 0.6 0]);
   end
   
end

yaxis(-0.5, numRows+1);
set(gca, 'TickDir', 'out', 'TickLen', 0.005*[1 1]);

box off;

tStateChange(end+1) = max(t);
for k = 1:2:length(S.History),
   ks = (k+1)/2;
   if ks < length(tStateChange),
      h = text(mean(tStateChange(ks:ks+1)), -0.25, strrep(S.History{k},'_','\_'));
      set(h, 'HorizontalAlignment', 'center');
   end
end
for k = 2:2:length(S.History)-1,
   ks = k/2 + 1;
   if ks < length(tStateChange),
      h = text(tStateChange(ks), length(names)+0.55, S.History{k});
      set(h, 'HorizontalAlignment', 'center');
   end
end

set(gcf, 'NumberTitle', 'off', 'Name', S.name);

set_axis_size([15 0.5*numRows]);
tight_figure();

zoom xon;