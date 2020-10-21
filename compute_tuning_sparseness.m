function [sparseness, sparseness_by_intensity] = compute_tuning_sparseness(RF, plot_graph)
%Isaacson code uses response peak amplitude

RF(RF < 0) = 0; %Not sure how this performs with negative numbers

%sparseness of full trace
N = size(RF,1)*size(RF,2); %N = number of tones
r = reshape(RF,[N,1]);
E1 = sum((r/N))^2;
E2 = sum(r.^2/N);
sparseness = (1 - E1/E2)/(1 - 1/N);


%sparseness at dB intensity
sparseness_by_intensity = nan(size(RF,1),1);
for i = 1:size(RF,1)
    r = RF(i,:);
    N = length(r); %N = number of tones
    E1 = sum((r/N))^2;
    E2 = sum(r.^2/N);
    sparseness_by_intensity(i) = (1 - E1/E2)/(1 - 1/N);
end

if plot_graph == 1
    
    figure; hold on
    subplot(1,2,1)
    imagesc(RF)
    ylabel('Intensity')
    xlabel('Frequency')
    
    subplot(1,2,2)
    plot(sparseness_by_intensity)
    camroll(-90)
    xlabel('Intensity')
    ylabel('Sparseness')
    %ylim([0 1])
    
    suptitle(['Total sparseness = ' num2str(sparseness)])
end