function [options, plotOptions] = setOptions
%[options, plotOptions] = setOptions
%
%This function sets the options for fitting behavioral data with psignifit,
%and also sets options for the plotting features.
%
%Written by ML Caras Dec 5 2016

%Set up fit options
options.sigmoidName = 'norm'; %use cumulative gaussian fit
options.expType = 'YesNo';
options.confP  = 0.95; %confidence level for confidence intervals

%Set up plotOptions structure
plotOptions.dataColor      = [0 0 0];
plotOptions.plotData       = 1;
plotOptions.lineColor      = plotOptions.dataColor;
plotOptions.lineWidth      = 2;
plotOptions.xLabel         = '\Delta Frequency in octaves from target';
plotOptions.yLabel         = '%go';
plotOptions.labelSize      = 15;
plotOptions.fontSize       = 10;
plotOptions.fontName       = 'Arial';
plotOptions.tufteAxis      = false;
plotOptions.plotAsymptote  = false;
plotOptions.plotThresh     = false;
plotOptions.aspectRatio    = false;
plotOptions.extrapolLength = .2;
plotOptions.CIthresh       = false;
plotOptions.dataSize       = 25;
plotOptions.dataColor = 'k';
plotOptions.lineColor = plotOptions.dataColor;
plotOptions.plotThresh     = true;               

end
