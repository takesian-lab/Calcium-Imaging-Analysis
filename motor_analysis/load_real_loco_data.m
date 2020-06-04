%% CLEAN AND CLEAR

clear
close all 
clc

%%

% Load Carolyn's data
% Locomotor_data_2 is the velocity trace
%   Column1 = timestamps 
%   Column3 = activity (cm/s)
load locomotor_data_carolyn

% Load Maryse's data
% There's 8 traces in the same format as Carolyn's
load locomotor_data_maryse

%% Visualize velocity trace

fs = 10; % Hz

time = Locomotor_data_2(:,1);
velocity = Locomotor_data_2(:,3);

figure; plot(time,velocity)

% Grab an example transient 
start = 30; % seconds
duration = 5; % seconds

figure; plot(time((start:start+duration))*fs,velocity((start:start+duration)))

%% Do we have a time stamp issue?

df = diff(locodata.trace3(:,1))
uniqueVals = uniquetol(df)
countmember(uniqueVals,df)
size(uniqueVals)
