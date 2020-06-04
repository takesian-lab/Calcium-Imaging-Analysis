%% CLEAN AND CLEAR

clear
close all 
clc

%% Test noisy transient

noise_std = 0.5; % cm/s

fs = 30; % Hz
% Generate a random transient
[time,transient] = createRandomComplexTransient(fs,noise_std);

% Velocity filter cutoff frequency
fc1 = 2; % Hz
% Acceleration filter cutoff frequency
fc2 = 0.5; % Hz
% LPF order
LP_filter_order = 8;
% return filter coefficients for velocity filtering
[b1,a1] = butter(LP_filter_order,fc1/(fs/2));
% return filter coefficients for acceleration filtering
[b2,a2] = butter(LP_filter_order,fc2/(fs/2));

transient_filtered = filtfilt(b1,a1,transient);
transient_poorly_filtered = filter(b1,a1,transient);

plot(time,transient,'Color','Black','LineWidth',1)
hold on 
% The derivative of filtered velocity is the filtered acceleration
transient_derivative_pre_filtered = gradient(transient_filtered);
% The derivative of velocity is the acceleration
transient_derivative = gradient(transient);
% Low pass filtered acceleration (Not Good)
transient_derivative_post_filtered = filtfilt(b2,a2,transient_derivative);

plot(time,transient_derivative_pre_filtered,'LineWidth',2,'Color','Green')
% plot(time,transient_derivative_post_filtered,'LineWidth',1,'Color','Green')
plot(time,transient_filtered,'Blue','LineWidth',2)
plot(time,transient_poorly_filtered,'--r','LineWidth',2)

% Horizontal Thresholds
yline(0);
yline(0.5*std(transient));
yline(1.5*std(transient));
xlim([time(1) time(end)])
