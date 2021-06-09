function [sp, sp_by_int, sp_by_freq, h] = compute_tuning_sparseness_v2(RF, plotFigure, units, freqs, ints, cellNumber)
% Calculate lifetime sparseness
%
% Isaacson papers: Lin et al 2019 PNAS Arousal Regulates Frequency
% Tuning... and Kato, Asinof, Isaacson 2017 Neuron Network-Level Control...
%
% Lifetime sparseness (Rolls and Tovee, 1995; Willmore and Tolhurst, 2001), which was calculated as:
% (1 - {[sum of Rj/N]^2 / [sum of Rj^2/N]}) / (1 - 1/N)
% where Rj was the response peak amplitude of the cell to tone j, and N was the total number of tones.
% 1 – Sp provides a measure of how much the response probability of a neuron was distributed equally among all tones
% (non-selective: 1 - Sp = 1) versus attributable entirely to one tone (highly selective: 1- Sp = 0).
% Since calculation of 1 - Sp does not rely on the TRF shape, it can also be used for quantifying the selectivity of
% inhibitory responses, which often do not have clear V-shape.
%
% Argument(s): 
%   RF_orig (X x Y matrix) - RF matrix of cell's response to frequency x intensity stim combinations
%   plotFigure - 0 or 1 to plot figure
%   units - 'df/f' or 'spikes' for plot labels
%   freqs (X vector) - list of frequencies from low to high corresponding to RF stim
%   ints (Y vector - list of intensities from high to low corresponding to RF stim
%   cellNumber (int) - suite2p cell number
%
% Returns:
%   sp - sparseness for full RF
%   sp_by_int - frequency sparseness by intensity
%   sp_by_freq - intensity sparseness by frequency
%   h - figure handle
% 
% Notes:
%
% TODO: 
% Search 'TODO'

%% Calculate sparseness of full trace not including NaNs
R = RF(~isnan(RF));
N = numel(R); %N = number of stimuli
E1 = sum((R/N))^2;
E2 = sum(R.^2/N);
sp = (1 - E1/E2)/(1 - 1/N);


%% Calculate freq. sparseness at each dB intensity
sp_by_int = nan(size(RF,1),1);
for i = 1:size(RF,1)
    R = RF(i,:);
    R = R(~isnan(R));
    if isempty(R)
        continue;
    end
    N = length(R); %N = number of tones
    E1 = sum((R/N))^2;
    E2 = sum(R.^2/N);
    sp_by_int(i) = (1 - E1/E2)/(1 - 1/N);
end

%% Calculate int. sparseness at each freq
sp_by_freq = nan(size(RF,2),1);
for i = 1:size(RF,2)
    R = RF(:,i);
    R = R(~isnan(R));
    if isempty(R)
        continue;
    end
    N = length(R); %N = number of tones
    E1 = sum((R/N))^2;
    E2 = sum(R.^2/N);
    sp_by_freq(i) = (1 - E1/E2)/(1 - 1/N);
end

%% Figure

if plotFigure
    
    h = figure; hold on
    subplot(2,2,3)
    imagesc(RF_orig); hold on
    ylabel('Intensity')
    xlabel('Frequency')
    set(gca,'YTick',1:length(ints))
    set(gca,'YTickLabel', ints)
    set(gca,'XTick',1:length(freqs))
    set(gca,'XTickLabel', freqs)
    
    subplot(2,2,4)
    plot(sp_by_int); hold on
    scatter(1:length(ints),sp_by_int)
    camroll(-90)
    xlabel('Intensity')
    ylabel('Freq. sparseness by Intensity')
    ylim([0 1])
    xlim([0.5 length(ints)+0.5])
    set(gca,'XTick',1:length(ints))
    set(gca,'XTickLabel', ints)
    
    subplot(2,2,1)
    plot(sp_by_freq); hold on
    scatter(1:length(freqs),sp_by_freq)
    xlabel('Frequency')
    ylabel('Int. sparseness by Frequency')
    ylim([0 1])
    xlim([0.5 length(freqs)+0.5])
    set(gca,'XTick',1:length(freqs))
    set(gca,'XTickLabel', freqs)
    
    suptitle(['Total sparseness = ' num2str(sp) ' (' units ') Cell - ' num2str(cellNumber)])

else
    h = [];
end

%% USE THIS CODE TO TEST FUNCTION (run in separate script)
% 
% plotFigure = 1;
% units = 'A.U.';
% freqs = [4, 5.7, 8, 11.3, 16, 22.6, 32, 45.2];
% ints = [80 70 60 50 40 30 20 10];
% cellNumber = 0;
% 
% RF = zeros(8); %Divide by zero error: nans
% RF = ones(8); %Sparseness = 0
% RF = -ones(8); %Same as above
% RF = eye(8); %Sparseness = 1 at all freqs and ints
% RF = flipud(eye(8)); %Same as above
% RF = -eye(8); %Same as above
% RF = rand(8);
% RF = -rand(8);
% RF = zeros(8); RF(randi(64,10,1)) = 1;
% 
% [sp, sp_by_int, sp_by_freq, h] = compute_tuning_sparseness_v2(RF, plotFigure, units, freqs, ints, cellNumber);

end