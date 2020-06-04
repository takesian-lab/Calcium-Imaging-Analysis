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

% Plot the original, noisy synthetic velocity transient
figure(1)
plot(time,transient,'Color','Black','LineWidth',1)
hold on
% Zero velocity axis
yline(0);
% x-axis plotting limit 
xlim([time(1) time(end)])

%% Design filters for the velocity and filter the original signal

% Velocity filter cutoff frequency
velocityCutoffFrequency = 2; % Hz

% Lowpass filter order (chosen arbitrarily)
LPFilterOrder = 12;

% We will now compute the filter coefficients for velocity filtering.
% We will use a butterworth filter since they are maximally flat in the passband 
% and stopband for a given filter order.
% NOTE: According to matlab conventions:
% [b, a] are the coefficients in the denominater and numerator of our filter's 
% transfer function, respectively.
[bVelocity, aVelocity] = butter(LPFilterOrder,velocityCutoffFrequency/(samplingFrequency/2));

% This is the filtered velocity signal
% NOTE: filtfilt() gives us a zero phase filter that is double the order of our
% butterworth filter (designed above).
transientFiltered = filtfilt(bVelocity,aVelocity,transient);

% Continue plotting on figure 1 (started above)
% Plot the filtered signal
plot(time,transientFiltered,'Blue','LineWidth',2)

%% Threshold Detection Parameters 

% We will be using standard deviation based thresholds, however
% detectThreshold() also supports 'STD' threshold definitions
thresholdType = 'Fixed'; %'STD';
% See Beaulieu-Laroche et al.
% Threshold A is always larger than Threshold B
thresholdA = 10.0;
thresholdB = 2.5;

% Calculate the standard deviation of the original velocity trace
velocityStd = std(transient);

%% Detect Threshold B (See Beaulieu-Laroche & Harnett (2019) for details)

% These are crossing indexes
[thresholdBCrossingIndices] = detectThreshold(transientFiltered, thresholdType, thresholdB);
% Translate indices into time stamps
thresholdBTimeStamps = time([thresholdBCrossingIndices]);

% Continue plotting on figure 1 (started above)
% Plot horizontal threshold crossings
if strcmp(thresholdType,'STD')
    yline(thresholdB*velocityStd,'Color','Black');
elseif strcmp(thresholdType,'Fixed')
    yline(thresholdB,'Color','Black');
end

% Continue plotting on figure 1 (started above)
% Plot vertical threshold crossings and crossing points
% Loop over transients
for ithTransient = 1:size(thresholdBTimeStamps,1)
    % For a given transient loop over the time stamps (threshold crossings)
    for ithTS = thresholdBTimeStamps(ithTransient,:)
        xline(ithTS,'-k');
        % Before plotting check which threshold type we are using
        if strcmp(thresholdType,'STD')
            plot(ithTS,thresholdB*velocityStd,'ko','MarkerFaceColor','Black')
        elseif strcmp(thresholdType,'Fixed')
            plot(ithTS,thresholdB,'ko','MarkerFaceColor','Black')
        end
    end
end

%% Detect Threshold A (See Beaulieu-Laroche & Harnett (2019) for details)

% These are crossing indexes
[thresholdACrossingIndices] = detectThreshold(transientFiltered, thresholdType, thresholdA);
% reshape zero crossings indices (flatten)
thresholdACrossingIndices = reshape(thresholdACrossingIndices,1,[]);
% Translate indices into time stamps
thresholdATimeStamps = time([thresholdACrossingIndices]);

% Continue plotting on figure 1 (started above)
% Plot horizontal threshold crossings
% Before plotting check which threshold type we are using
if strcmp(thresholdType,'STD')
    yline(thresholdA*velocityStd,'Color','Red');
elseif strcmp(thresholdType,'Fixed')
    yline(thresholdA,'Color','Red');
end

% Continue plotting on figure 1 (started above)
% Plot vertical threshold crossings and crossing points
% Loop over time stamps (threshold crossings)
for ithTS = thresholdATimeStamps
    xline(ithTS,'-r');
    % Before plotting check which threshold type we are using
    if strcmp(thresholdType,'STD')
        plot(ithTS,thresholdA*velocityStd,'ro','MarkerFaceColor','Red')
    elseif strcmp(thresholdType,'Fixed')
        plot(ithTS,thresholdA,'ro','MarkerFaceColor','Red')
    end
end

%% Compute the acceleration (the derivative of velocity) and filter it

% NOTE: This is a little different than the methods outlined in
% Beaulieu-Laroche & Harnett (2019)
% 
% For calcium transients the authors compute the derivative on the
% filtered signal and then lowpass the derivative after the fact (with a
% cutoff frequency of 0.5 Hz).  
% 
% For locomotion: velocity was low-pass filtered at 0.05 Hz, and 2.5 cm/s 
% was used as a 'Fixed' threshold to detect moving epochs. They did not
% calculate acceleration. Nor did they segment motor bouts in the same
% way they did calcium transients.
% 
% In contrast, here we compute the acceleration on the original velocity signal, 
% then use the same cutoff frequency for filtering both the velocity and acceleration (2 Hz).
% 
% Differences in filter cutoff frequencies are likely due to differences in sampling 
% frequency.  The sampling rate for running speed used by Beaulieu-Laroche et al. is unknown.

% Acceleration filter cutoff frequency
accelerationCutoffFrequency = velocityCutoffFrequency; % Hz

% We will now compute the filter coefficients for acceleration filtering.
% We will use a butterworth filter since they are maximally flat in the passband 
% and stopband for a given filter order.
% NOTE: According to matlab conventions:
% [b, a] are the coefficients in the denominater and numerator of our filter's 
% transfer function, respectively.
[bAcceleration, aAcceleration] = butter(LPFilterOrder,accelerationCutoffFrequency/(samplingFrequency/2));

% Compute the derivative of the velocity, the acceleration.
acceleration = gradient(transientFiltered); % Alternatively: gradient(transient);
% Then we filter the acceleration
transientFilteredAcceleration = filtfilt(bAcceleration,aAcceleration,acceleration); % Alternatively: acceleration;

% Continue plotting on figure 1 (started above)
% Plot the acceleration
plot(time,transientFilteredAcceleration,'LineWidth',2,'Color','Green')

%% Calculate acceleration zero crossing time stamps

% These are zero crossing indices
[accelerationZeroCrossingIndices] = detectThreshold(transientFilteredAcceleration, thresholdType, 0);
% reshape zero crossings indices (flatten)
accelerationZeroCrossingIndices = reshape(accelerationZeroCrossingIndices,1,[]);
% Translate indices into time stamps
accelerationZeroCrossingTimeStamps = time([accelerationZeroCrossingIndices]);

% Continue plotting on figure 1 (started above)
% Plot vertical threshold crossings and crossing points
% Loop over time stamps (zero crossings)
for ithTS = accelerationZeroCrossingTimeStamps
    xline(ithTS,'-g');
    % Before plotting check which threshold type we are using
    if strcmp(thresholdType,'STD')
        plot(ithTS,0,'go','MarkerFaceColor','Green')
    elseif strcmp(thresholdType,'Fixed')
        plot(ithTS,0,'go','MarkerFaceColor','Green')
    end
end

%% Temporal Segmentation of the velocity transient
% Here is where we define our "motor bouts" and the "motor onset transients" 
% We will assign the time stamps associated with each, and segment the
% original velocity trace using them.

clc

thresholdATimeStamps
thresholdBTimeStamps
accelerationZeroCrossingTimeStamps

% This is the number of locomotion events that might qualify as a motor bout
numberOfCandidateTransients = size(thresholdBTimeStamps,1);
