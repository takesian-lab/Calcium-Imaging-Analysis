function [velocity_activity_time_stamps,synthetic_velocity_activity] = createRandomComplexTransient(fs,noise_std)
%CreateRandomVelocityTransient Summary of this function goes here
%   Detailed explanation goes here

disp('Generating a random transient...')

%% Sampling

% sampling period
dt = 1/fs;

%% Generate random transient duration

% How much baseline activity before and after the transient
% TODO: this could be separated and/or randomized
synthetic_baseline_duration = 0.5; % seconds

% Create upper and lower bounds for the transient duration
min_duration = 1;   % seconds
max_duration = 5;   % seconds
% Generate a random duration for the transient (seconds)
random_transient_duration = (max_duration-min_duration).*rand() + min_duration; % seconds
% Total bevahior trace duration (seconds)
synthetic_velocity_activity_duration = 2*synthetic_baseline_duration + random_transient_duration; % seconds

%% Create behavioral time stamps

random_step_height = 1;

baseline_activity = zeros(1,length(0:dt:synthetic_baseline_duration)-1);
transient_step_activity = random_step_height*ones(1,length(0:dt:random_transient_duration));

% create time stamps
velocity_activity_time_stamps = 0:dt:synthetic_velocity_activity_duration;
% extract transient time stamps
velocity_transient_time_stamps = velocity_activity_time_stamps(numel(baseline_activity)+1:numel(baseline_activity)+numel(transient_step_activity));

baseline_activity = baseline_activity + normrnd(0, noise_std, 1,length(baseline_activity));
% create a tapered to baseline transient that is normalized to amplitude one
baseline_activity = baseline_activity - min(baseline_activity);

num_transient_time_steps = numel(velocity_transient_time_stamps);
disp(['Number of transient time steps: ' num2str(num_transient_time_steps)])

% plot(velocity_activity_time_stamps,synthetic_baseline_step_activity)

%% Generate a random number of gaussians to mix
% We will use this to create smoothly but randomly varying synthetic
% velocity transients

min_num_gaussians = 5;  % A finite number of gaussians
max_num_gaussians = 100;% A finite number of gaussians
% This is how many gaussian we will use
random_number_gaussians = min(randi([min_num_gaussians max_num_gaussians],1),floor(num_transient_time_steps/2)); % A finite number of gaussians
disp(['Number of random gaussians: ' num2str(random_number_gaussians)])

%% Generate a velocity transient event

% Largest possible temporal pad
max_pad = 5; % number of time stamps
% how many time stamps to pad the transient with?
pad = min(floor(num_transient_time_steps/2)-1,max_pad); % time stamps
velocity_transient_time_stamps_padded = velocity_transient_time_stamps(pad+1:num_transient_time_steps-pad);
num_velocity_transient_time_stamps_padded = numel(velocity_transient_time_stamps_padded);
disp(['Number of padded transient time steps: ' num2str(num_velocity_transient_time_stamps_padded)])

% Randomly choose means (time stamps) during the transient
% We are creating random transient events peaks (gaussians centered at
% these times)
random_indices = randi(num_velocity_transient_time_stamps_padded,1,random_number_gaussians);
mus = velocity_transient_time_stamps_padded(random_indices);

% Create upper and lower bounds for the transient standard deviations (in time)
min_std = 0.1;   % Unit: std
max_std = 0.5;   % Unit: std

% pre-allocate transient activity trace
synthetic_velocity_activity = zeros(1,num_transient_time_steps);
% Loop over all transient events
for mu = mus
    % Generate random standard deviation
    random_std = (max_std-min_std).*rand() + min_std; % Unit: std
    % sum gaussian transient peaks
    synthetic_velocity_activity = synthetic_velocity_activity + normpdf(velocity_transient_time_stamps,mu,random_std).*rand();
end

% create a tapered to baseline transient that is normalized to amplitude one
synthetic_velocity_activity = synthetic_velocity_activity - min(synthetic_velocity_activity);
% To smooth out the transition:
% We point-wise multiply the transient with a normal distribution with a mean
% centered in the transient (in time) with std 1
synthetic_velocity_activity = synthetic_velocity_activity.*normpdf(velocity_transient_time_stamps,velocity_transient_time_stamps(ceil(num_transient_time_steps/2)),1);
synthetic_velocity_activity = synthetic_velocity_activity/max(synthetic_velocity_activity);

% Create upper and lower bounds for the velocity peak amplitude (m/s)
min_peak_velocity = 25;   % Unit: m/s
max_peak_velocity = 35;   % Unit: m/s
% velocity peak amplitude (m/s)
random_peak_velocity = (max_peak_velocity-min_peak_velocity).*rand() + min_peak_velocity; % Unit: m/s
synthetic_velocity_activity = random_peak_velocity.*synthetic_velocity_activity;

% Add baseline (zero) velocity activity before and after the transient
synthetic_velocity_activity = [baseline_activity synthetic_velocity_activity baseline_activity];

synthetic_velocity_activity = abs(synthetic_velocity_activity + normrnd(0, noise_std, 1,length(synthetic_velocity_activity)));
% % Makes sure the velocity is always positive
% synthetic_velocity_activity = synthetic_velocity_activity - min(synthetic_velocity_activity);

end

