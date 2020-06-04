%% CLEAN AND CLEAR

clear
close all 
clc

%% NOTES

% This script is meant to demonstrate by example(s) the use of functions we
% developed for the analysis of transients in time series data. 
% In this case, we will apply these functions to synthetic velocity data. 
% 
% However, these functions are general and could be applied to any time series, 
% most notably they could be directly applied to the analysis of transients
% in calcium imaging data (e.g. GCaMP dF/F time series).

% This script is loosely inspired by the methods for transient detection outlined in 
% Beaulieu-Laroche & Harnett (2019).
% 
% CITATION:
% Beaulieu-Laroche L, Toloza EHS, Brown NJ, Harnett MT (2019). Widespread and highly 
% correlated somato-dendritic activity in cortical layer 5 neurons. Neuron 103(2):235-241

% We will be using two thresholds (A and B) to define our transients
% in time.  See the methods section from the above reference for the meaning of these thresholds.

% Real velocity data is sampled at 10 Hz and will need to be upsampled to
% 30 Hz.  Synthetic velocity data is generated at a sampling frequency of 30 Hz.

%% Create and filter a synthetic noisy transient

% This is the variance of the noise we will add to the synthetic velocity
% time series data (chosen arbitrarily)
noiseStd = 0.5; % cm/s

% This is the sample rate for the synthetic velocity time series
% Note: Real velocity data is sampled at 10 Hz and will need to be upsampled
samplingFrequency = 30; % Hz

% Generate a random velocity transient
% CreateRandomComplexTransient() creates a transient of random duration 
% comprised of a mixture of a random number of gaussian peaks 
% with random peak locations in time with random variance.
% See the documentation for CreateRandomComplexTransient() for more details
[time,transient] = createRandomComplexTransient(samplingFrequency,noiseStd);

%% Segment the synthetic velocity data

thresholdA = 10;  % cm/s
thresholdB = 2.5; % cm/s

motorData = segmentMotorData(time,transient,samplingFrequency,'Fixed',thresholdA,thresholdB);

%% Visualize the data and segmentation

figure;
% plot the original velocity datas
plot(motorData.velocityTimeStamps,motorData.velocity,'LineWidth',2,'Color','Black')
hold on 
% loop over all motor bouts in the trace
for ithMotorBout = 1:size(motorData.motorBouts,1)
    % plot the full motor bout segment
    plot(motorData.motorBoutTimeStamps{ithMotorBout,:},motorData.motorBouts{ithMotorBout,:},'LineWidth',3,'Color','Blue')
    % marker for the full motor bout segment
    plot(motorData.motorBoutTimeStamps{ithMotorBout,:},-2*ones(1,length(motorData.motorBouts{ithMotorBout,:})),'LineWidth',3,'Color','Blue')
    % plot the rise event segment
    plot(motorData.riseEventTimeStamps{ithMotorBout,:},motorData.riseEvents{ithMotorBout,:},'LineWidth',3,'Color','Red')
    % marker for the rise event segment
    plot(motorData.riseEventTimeStamps{ithMotorBout,:},-1*ones(1,length(motorData.riseEvents{ithMotorBout,:})),'LineWidth',3,'Color','Red')
end
% Zero velocity axis
yline(0);
% threshold A
yline(thresholdA);
% threshold B
yline(thresholdB);
% x-axis plotting limit 
xlim([motorData.velocityTimeStamps(1) motorData.velocityTimeStamps(end)])
% Labels & titles
title('Velocity Data Segmentation', 'FontSize', 30)
xlabel('Time (seconds)', 'FontSize', 24)
ylabel('Speed (cm/second)', 'FontSize', 24)
ax = gca;
ax.FontSize = 16; 
% legend('Location','south')
% 'DisplayName','Example',
