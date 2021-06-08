function [x,fitted_yes,fitted_dprime,threshold,slope] = ...
    plotPsych_dprime(results,dprimemat,options,plotOptions,zHR)
%[x,fitted_yes,fitted_dprime,threshold,slope] = ...
%    plotPsych_dprime(results,dprimemat,options,plotOptions,zFA)
%
%This function plots psychometric fits after transforming to dprime space.
%Input variables:
%   results: results structure created by psignifit 4
%   dprimemat: M x 2 matrix arranged as [stimulus (dB re: 100%), dprime]
%   options: structure created by setOptions function for
%            generating fits with psignifit
%   plotOptions: structure created by setOptions function for setting plot
%             features
%   zFA: scored FA rate
%
%
%Written by ML Caras Dec 5 2016





%Calculate threshold and slope
threshold = getThreshold(results,options.threshPC,1) %scaled threshold
slope = getSlopePC(results,options.threshPC,1); %scaled slope


%Establish x values
xlength = max(results.data(:,1))-min(results.data(:,1));
xLow = min(results.data(:,1))- plotOptions.extrapolLength*xlength;
xHigh = max(results.data(:,1))+ plotOptions.extrapolLength*xlength;
x  = linspace(xLow,xHigh,1000);

%Define percent "yes" responses (y values) for fit
fitted_yes = (1-results.Fit(3)-results.Fit(4))*arrayfun(@(x)...
    results.options.sigmoidHandle(x,results.Fit(1),...
    results.Fit(2)),x)+results.Fit(4);

%Define dprime (y values) for fit - This is the CGS edit 6/2/21
fitted_dprime  = zHR - (sqrt(2)*erfinv(2*fitted_yes-1));


%Plot raw dprime values
plot(dprimemat(:,1),dprimemat(:,2),'.',...
    'MarkerSize',plotOptions.dataSize,...
    'Color',plotOptions.dataColor);
hold on;

%Plot dprime fit
plot(x,fitted_dprime,'Color', plotOptions.lineColor,...
    'LineWidth',plotOptions.lineWidth);


%Plot threshold lines
ylimits = get(gca,'ylim');
ymin = ylimits(1);
x_thresh = [threshold,threshold];
y = [ymin,options.dprimeThresh];
plot(x_thresh,y,'-','Color',plotOptions.lineColor);

% xlimits = get(gca,'xlim');
% xmin = xlimits(1);
% x = [xmin,threshold];
% y = [options.dprimeThresh,options.dprimeThresh];
% plot(x_thresh,y,'-','Color',plotOptions.lineColor);


%Format plot
axis tight
set(gca,'FontSize',plotOptions.fontSize)

ylabel('d''','FontName',plotOptions.fontName,...
    'FontSize', plotOptions.labelSize);
xlabel(plotOptions.xLabel,'FontName',plotOptions.fontName,...
    'FontSize', plotOptions.labelSize);

ylim([0,3.5]);
set(gca,'TickDir','out')
box off

end