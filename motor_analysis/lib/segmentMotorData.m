function [motorData] = segmentMotorData(timeStamps,velocityTimeSeries,samplingFrequency,thresholdType,thresholdA,thresholdB)
% segmentMotorData(timeStamps,velocityTimeSeries,samplingFrequency,thresholdType,thresholdA,thresholdB) 
% 
% segmentMotorData extracts motor activity bouts and onsets and the associated time stamps 
% from velocity time series data.
% 
% TODO: Documentation in progress
% 
% Written by Wisam Reid - June 2020 - wisam@g.harvard.edu
%
% ARGUMENTS:       
% 
%          timeStamps:  (double) vector
%                       Time stamps for the velocity time series.
% 
%  velocityTimeSeries:  (double) vector
%                       A velocity time series.
%                       Note: This time series should not be smoothed
% 
%   samplingFrequency:  (double) 
% 
%       thresholdType:  (string) 
%                       options: {'STD','Fixed'} 
% 
%          thresholdA:  (double) 
%                       A multiplier of the standard deviation of input signal or 
%                       a fixed threshold value
% 
%          thresholdB:  (double) 
%                       A multiplier of the standard deviation of input signal or 
%                       a fixed threshold value
% 
% RETURNS: 
% 
%           motorData:  (struct)
%                                   
%                       motorData has the following fields:
% 
%                       1) motorData.velocity 
%                       The original velocity trace 
% 
%                       2) motorData.motorBoutTimeStamps
%                       Entire motor bouts (time stamps between a positive and negative pair of 
%                       crossings through threshold B, conditioned on a threshold A crossing within 
%                       the interval) 
% 
%                       3) motorData.riseEventTimeStamps
%                       Motor rise events (time stamps between a positive threshold B crossing 
%                       and the first zero crossing of the derivative (i.e., the acceleration is zero), 
%                       after the first threshold A crossing) 
% 
%                       4) motorData.nonmotorBoutTimeStamps
%                       All non-motor bouts (all other time stamps that do not meet the above 
%                       criteria for a locomotor bout) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% NOTES:
% 
% This script is loosely inspired by the methods for transient detection outlined in 
% Beaulieu-Laroche & Harnett (2019).
% 
% CITATION:
% Beaulieu-Laroche L, Toloza EHS, Brown NJ, Harnett MT (2019). Widespread and highly 
% correlated somato-dendritic activity in cortical layer 5 neurons. Neuron 103(2):235-241
% 
% We will be using two thresholds (A and B) to define our transients
% in time.  See the methods section from the above reference for the meaning of these thresholds.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example Use:
% 
% % Create synthetic transient
% noiseStd = 0; % cm/s
% samplingFrequency = 30; % Hz
% [time,transient] = createRandomComplexTransient(samplingFrequency,noiseStd);
% 
% % Segment synthetic transient
% thresholdA = 10;  % cm/s
% thresholdB = 2.5; % cm/s
% motorData = segmentMotorData(time,transient,samplingFrequency,'Fixed',thresholdA,thresholdB);
% 
% % Plot Segmentation
% figure;
% % plot the original velocity datas
% plot(motorData.velocityTimeStamps,motorData.velocity,'LineWidth',2,'Color','Black')
% hold on 
% % loop over all motor bouts in the trace
% for ithMotorBout = 1:size(motorData.motorBouts,1)
%     % plot the full motor bout segment
%     plot(motorData.motorBoutTimeStamps{ithMotorBout,:},motorData.motorBouts{ithMotorBout,:},'LineWidth',3,'Color','Blue')
%     % marker for the full motor bout segment
%     plot(motorData.motorBoutTimeStamps{ithMotorBout,:},-2*ones(1,length(motorData.motorBouts{ithMotorBout,:})),'LineWidth',3,'Color','Blue')
%     % plot the rise event segment
%     plot(motorData.riseEventTimeStamps{ithMotorBout,:},motorData.riseEvents{ithMotorBout,:},'LineWidth',3,'Color','Red')
%     % marker for the rise event segment
%     plot(motorData.riseEventTimeStamps{ithMotorBout,:},-1*ones(1,length(motorData.riseEvents{ithMotorBout,:})),'LineWidth',3,'Color','Red')
% end
% % Zero velocity axis
% yline(0);
% % threshold A
% yline(thresholdA);
% % threshold B
% yline(thresholdB);
% % x-axis plotting limit 
% xlim([motorData.velocityTimeStamps(1) motorData.velocityTimeStamps(end)])
% % Labels & titles
% title('Velocity Data Segmentation', 'FontSize', 30)
% xlabel('Time (seconds)', 'FontSize', 24)
% ylabel('Speed (cm/second)', 'FontSize', 24)
% ax = gca;
% ax.FontSize = 16; 
%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEBUGGING FLAGS

% Plot the results from the threshold detection
debug_transient_detection = 0;

%% parse arguments

%% Save the original velocity data and time stamps

motorData.velocity = velocityTimeSeries;
motorData.velocityTimeStamps = timeStamps;

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
transientFiltered = filtfilt(bVelocity,aVelocity,velocityTimeSeries);

%% Threshold Detection Parameters 

% % Threshold A is always larger than Threshold B (See Beaulieu-Laroche & Harnett (2019) for details)
% thresholdA = 10.0;
% thresholdB = 2.5;

% % Calculate the standard deviation of the original velocity trace
% velocityStd = std(transient);

%% Detect Threshold B (See Beaulieu-Laroche & Harnett (2019) for details)

% These are crossing indexes
[thresholdBCrossingIndices] = detectThreshold(transientFiltered, thresholdType, thresholdB);
% Translate indices into time stamps
thresholdBTimeStamps = timeStamps([thresholdBCrossingIndices]);

%% Detect Threshold A (See Beaulieu-Laroche & Harnett (2019) for details)

% These are crossing indexes
[thresholdACrossingIndices] = detectThreshold(transientFiltered, thresholdType, thresholdA);
% reshape zero crossings indices (flatten)
thresholdACrossingIndices = reshape(thresholdACrossingIndices,1,[]);
% Translate indices into time stamps
thresholdATimeStamps = timeStamps([thresholdACrossingIndices]);

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
accelerationZeroTSs = gradient(transientFiltered); % Alternatively: gradient(transient);
% Then we filter the acceleration
transientFilteredAcceleration = filtfilt(bAcceleration,aAcceleration,accelerationZeroTSs); % Alternatively: acceleration;

%% Calculate acceleration zero crossing time stamps

% These are zero crossing indices
[accelerationZeroCrossingIndices] = detectThreshold(transientFilteredAcceleration, thresholdType, 0);
% reshape zero crossings indices (flatten)
accelerationZeroCrossingIndices = reshape(accelerationZeroCrossingIndices,1,[]);
% Translate indices into time stamps
accelerationZeroCrossingTimeStamps = timeStamps([accelerationZeroCrossingIndices]);

%% Temporal Segmentation of the velocity transients
% Here is where we define our "motor bouts" and the "motor onset transients" 
% We will assign the time stamps associated with each, and segment the
% original velocity trace using them.

% This is the number of locomotion events that might qualify as a motor bout
numberOfCandidateMotorBouts = size(thresholdBTimeStamps,1);

% Loop over candidate motor bouts
for ithMotorBout = 1:numberOfCandidateMotorBouts
    % Extract the starting and ending time stamps (TSs) for the current candidate motor bout.
    % These are the start (positive threshold B crossing) and the end
    % (negative threshold B crossing) of the motor bout.
    % There should always be 2 TSs here.
    thresholdBTSs = thresholdBTimeStamps(ithMotorBout,:);
    % Grab the time stamps for all threshold A crossings within the current
    % candidate motor bout time interval.
    thresholdATSs = thresholdATimeStamps(thresholdATimeStamps > thresholdBTSs(1) & thresholdATimeStamps < thresholdBTSs(2));
    % Is this a motor bout? 
    if ~isempty(thresholdATSs)
        % This is now considered to be a motor bout and 
        % This is the starting TS for the current motor bout
        motorBoutStartTS = thresholdBTSs(1);
        % This is the starting index for the current motor bout
        motorBoutStartIdx = find(timeStamps == motorBoutStartTS);
        % This is the ending TS for the current motor bout
        motorBoutEndTS = thresholdBTSs(2);
        % This is the ending index for the current motor bout
        motorBoutEndIdx = find(timeStamps == motorBoutEndTS);
        % Now we can also define the beginning TS of the rise event for the current motor bout
        riseEventStartTS = thresholdBTSs(1);
        % This is the starting index for the rise event for the current motor bout
        riseEventStartIdx = find(timeStamps == riseEventStartTS);
        % Grabbing all TSs where the acceleration is zero and the velocity is above threshold A 
        accelerationZeroTSs = sort(accelerationZeroCrossingTimeStamps(accelerationZeroCrossingTimeStamps >= thresholdATSs(1) & accelerationZeroCrossingTimeStamps <= thresholdATSs(2)),'ascend');
        % Now we can define the end of the rise event for this bout as the
        % first time the acceleration is zero and the velocity is above threshold A
        % (this will be in between the threshold A crossings).
        riseEventEndTS = accelerationZeroTSs(1);
        % This is the starting index for the rise event for this bout
        riseEventEndIdx = find(timeStamps == riseEventEndTS);
        % Now we grab all the TSs for the entire motor bout
        motorData.motorBoutTimeStamps{ithMotorBout,:} = timeStamps(timeStamps >= motorBoutStartTS & timeStamps <= motorBoutEndTS);
        motorData.motorBouts{ithMotorBout,:} = velocityTimeSeries(motorBoutStartIdx:motorBoutEndIdx);
        % Now we grab all the TSs for the current rise event
        motorData.riseEventTimeStamps{ithMotorBout,:} = timeStamps(timeStamps >= riseEventStartTS & timeStamps <= riseEventEndTS);
        motorData.riseEvents{ithMotorBout,:} = velocityTimeSeries(riseEventStartIdx:riseEventEndIdx);
 
    end
end
    
%% RETURN

%% VISUALIZE/DEBUG

if debug_transient_detection
    
    motorData.thresholdATimeStamps = thresholdATimeStamps;
    motorData.thresholdBTimeStamps = thresholdBTimeStamps;
    motorData.accelerationZeroCrossingTimeStamps = accelerationZeroCrossingTimeStamps;
    
    % Calculate the standard deviation of the original velocity trace
    velocityStd = std(motorData.velocity);
    
    % Plot the original, noisy synthetic velocity transient
    figure;
    plot(motorData.velocityTimeStamps,motorData.velocity,'Color','Black','LineWidth',1)
    hold on
    % Zero velocity axis
    yline(0);
    % x-axis plotting limit 
    xlim([motorData.velocityTimeStamps(1) motorData.velocityTimeStamps(end)])
    % Plot the filtered signal
    plot(motorData.velocityTimeStamps,transientFiltered,'Blue','LineWidth',2)
    % Plot horizontal threshold crossings
    if strcmp(thresholdType,'STD')
        yline(thresholdB*velocityStd,'Color','Black');
    elseif strcmp(thresholdType,'Fixed')
        yline(thresholdB,'Color','Black');
    end
    % Plot vertical threshold crossings and crossing points
    % Loop over transients
    for ithTransient = 1:size(motorData.thresholdBTimeStamps,1)
        % For a given transient loop over the time stamps (threshold crossings)
        for ithTS = motorData.thresholdBTimeStamps(ithTransient,:)
            xline(ithTS,'-k');
            % Before plotting check which threshold type we are using
            if strcmp(thresholdType,'STD')
                plot(ithTS,thresholdB*velocityStd,'ko','MarkerFaceColor','Black')
            elseif strcmp(thresholdType,'Fixed')
                plot(ithTS,thresholdB,'ko','MarkerFaceColor','Black')
            end
        end
    end
    % Plot horizontal threshold crossings
    % Before plotting check which threshold type we are using
    if strcmp(thresholdType,'STD')
        yline(thresholdA*velocityStd,'Color','Red');
    elseif strcmp(thresholdType,'Fixed')
        yline(thresholdA,'Color','Red');
    end
    % Plot vertical threshold crossings and crossing points
    % Loop over time stamps (threshold crossings)
    for ithTS = motorData.thresholdATimeStamps
        xline(ithTS,'-r');
        % Before plotting check which threshold type we are using
        if strcmp(thresholdType,'STD')
            plot(ithTS,thresholdA*velocityStd,'ro','MarkerFaceColor','Red')
        elseif strcmp(thresholdType,'Fixed')
            plot(ithTS,thresholdA,'ro','MarkerFaceColor','Red')
        end
    end
    % Plot the acceleration
    plot(motorData.velocityTimeStamps,transientFilteredAcceleration,'LineWidth',2,'Color','Green')
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
    % Labels & titles
    title('Velocity Threshold Detection', 'FontSize', 30)
    xlabel('Time (seconds)', 'FontSize', 24)
    ylabel('Speed (cm/second)', 'FontSize', 24)
    ax = gca;
    ax.FontSize = 16; 
    
end

end