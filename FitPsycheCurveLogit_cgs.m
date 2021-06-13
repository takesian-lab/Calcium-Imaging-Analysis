% http://matlaboratory.blogspot.co.uk/2015/04/introduction-to-psychometric-curves-and.html

function [coeffs, stats, curve, threshold] = ...
    FitPsycheCurveLogit(xAxis, yData, weights, targets)

% Transpose if necessary
if size(xAxis,1)<size(xAxis,2)
    xAxis=xAxis';
end
if size(yData,1)<size(yData,2)
    yData=yData';
end
if size(weights,1)<size(weights,2)
    weights=weights';
end

% Perform fit
[coeffs, ~, stats] = ...
    glmfit(xAxis, [yData, weights], 'binomial','link','logit');

% Create a new xAxis with higher resolution
fineX = linspace(min(xAxis),max(xAxis),numel(xAxis)*50);
% Generate curve from fit
curve = glmval(coeffs, fineX, 'logit');
if max(weights)<=1
    % Assume yData was proportional
    curve = [fineX', curve];
else
    % Assume yData was % or actual number of trials
    curve = [fineX', curve*100];
end

% If targets (y) supplied, find threshold (x), else find 25, 50 and 75%
% values
if nargin==4
else
    targets = [0.25, 0.5, 0.75];
end
% Calculate
threshold = (log(targets./(1-targets))-coeffs(1))/coeffs(2);