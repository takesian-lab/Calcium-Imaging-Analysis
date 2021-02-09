function visualize_zcorr(block)
% DOCUMENTATION IN PROGRESS
%
% This function allows you to preview the output of ops.zcorr from a single block
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   
% 
% Notes:
%
%
% TODO: 
% Search 'TODO'

%% Setup

zcorr = block.ops.zcorr(1:100,:);
[best_z_val, best_z] = max(zcorr,[],1);
recording_depth = mode(best_z);
timestamp = block.timestamp;

%locomotor activity
loco_time = block.loco_times;
loco_speed = block.loco_activity;
        
%% PLOT

figure;
subplot(3,1,1)
imagesc(zcorr)
ylabel('Z stack position')
xlabel('Frames')

subplot(3,1,2)
plot(best_z)
hline(recording_depth, 'r')
ylabel('Best Z stack position')
xlabel('Frames')
xlim([1, length(best_z)])

subplot(3,1,3)
plot(loco_time, loco_speed)
xlim([0 timestamp(end)]) %loco will be shorter because tosca ends before PV

end