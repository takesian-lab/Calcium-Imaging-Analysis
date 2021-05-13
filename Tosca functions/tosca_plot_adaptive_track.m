function tosca_plot_adaptive_track(d, track)
% TOSCA_PLOT_ADAPTIVE_TRACK -
%

if nargin < 2, track = 1; end

L = [];
trialType = [];
result = [];

for k = 1:length(d),
   if d{k}.block ~= track, continue; end

   L(end+1) = d{k}.cue.Speaker_1.Level.dB_SPL; 
   trialType(end+1) = strcmpi(d{k}.Type, 'adapt');
   result(end+1) = strcmpi(d{k}.Result, 'go');
   
end

figure;
hold on;
idx = find(trialType == 1);
plot(idx, L(idx), 'bo-');
idx = find(trialType==1 & result == 1);
plot(idx, L(idx), 'bo', 'MarkerFaceColor', 'b');
