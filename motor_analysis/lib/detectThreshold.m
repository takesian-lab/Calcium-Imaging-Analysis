function [crossings] = detectThreshold(signal, thresholdType, threshold)
% detectThreshold finds the domain index where the input signal crosses a user defined threshold.
% 
% Written by Wisam Reid - June 2020 - wisam@g.harvard.edu
%  
% This function locates both upward and downward crossings by looking for a 
% sign change in the threshold shifted signal (thresholdShiftedSignal). 
% Where, thresholdShiftedSignal = signal - threshold 
% 
% This function is vectorized to optimize speed.
%
% Arguments:       
%                  signal: (double) vector
%           thresholdType: (string) options: {'STD','Fixed'} 
%               threshold: (double) multiplier of the std of signal or a fixed threshold
%
% Returns:      
%               crossings: (double) 2 column array where the first column contains 
%                          upward crossings and the second contains downward crossings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example Use:
%              
% dt = 0.001; % (1 ms in seconds)
% time = 0:dt:20; % (seconds)
% triangle = [0:dt:10 flip(0:dt:10-dt)];
% threshold = 3.01;
% crossings = detectThreshold(triangle, 'Fixed', threshold);
% figure; % plot demo
% plot(time,triangle,'-b','DisplayName','Signal','LineWidth',2)
% hold on
% yline(threshold,'-g','DisplayName','Threshold','LineWidth',2);
% xline(time(crossings(1)),'-r','DisplayName','Positive Crossing','LineWidth',2);
% xline(time(crossings(2)),'-k','DisplayName','Negative Crossing','LineWidth',2);
% p = plot(time(crossings(1)),threshold,'ro','MarkerFaceColor','Red');
% p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% p = plot(time(crossings(2)),threshold,'ko','MarkerFaceColor','Black');
% p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% xlim([time(1) time(end)]);
% title('detectThreshold Demo')
% xlabel('Time (seconds)') 
% ylabel('Amplitude (generic units)')
% legend('Location','south')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Are we using a standard deviation based or fixed threshold?
switch thresholdType
    case 'STD'
        % The threshold as a factor of the std of the input signal
        threshold = threshold*std(signal);
    case 'Fixed'
        threshold = threshold;
end

% Compute the thresholded shifted signal. 
% We will use this to look for sign changes.
thresholdShiftedSignal = signal - threshold;

% Now we will locate the crossings by looking for where consecutive elements of
% signal change sign 
% i.e. their product is negative
crossings = find(thresholdShiftedSignal.*circshift(thresholdShiftedSignal,1) < 0);

% Now we will pair the upward and downward crossing into a cell array:

% 1) Reshape the array into a 2 column array where the first column
% contains upward crossings and the second contains downward crossings
crossings = reshape(crossings,2,[])';

% % % % 2) Combine the upward and downward crossings into each element of the cell array 
% % % crossings = num2cell(crossings,2);

end