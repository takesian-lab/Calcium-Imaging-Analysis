function [sp, sp_by_int, sp_by_freq] = compute_tuning_sparseness(RF, plot_graph)
%Isaacson code uses response peak amplitude

RF_orig = RF;

RF(RF < 0) = 0; %Not sure how this performs with negative numbers
RF(isnan(RF)) = 0;

%sparseness of full trace
N = size(RF,1)*size(RF,2); %N = number of tones
r = reshape(RF,[N,1]);
E1 = sum((r/N))^2;
E2 = sum(r.^2/N);
sp = (1 - E1/E2)/(1 - 1/N);


%sparseness at dB intensity
sp_by_int = nan(size(RF,1),1);
for i = 1:size(RF,1)
    r = RF(i,:);
    N = length(r); %N = number of tones
    E1 = sum((r/N))^2;
    E2 = sum(r.^2/N);
    sp_by_int(i) = (1 - E1/E2)/(1 - 1/N);
end

%sparseness at freq
sp_by_freq = nan(size(RF,2),1);
for i = 1:size(RF,2)
    r = RF(:,i);
    N = length(r); %N = number of tones
    E1 = sum((r/N))^2;
    E2 = sum(r.^2/N);
    sp_by_freq(i) = (1 - E1/E2)/(1 - 1/N);
end

if plot_graph == 1
    
    figure; hold on
    subplot(2,2,3)
    imagesc(RF_orig)
    ylabel('Intensity')
    xlabel('Frequency')
    
    subplot(2,2,4)
    plot(sp_by_int)
    camroll(-90)
    xlabel('Intensity')
    ylabel('Sparseness per Intensity')
    %ylim([0 1])
    
    subplot(2,2,1)
    plot(sp_by_freq)
    xlabel('Frequency')
    ylabel('Sparseness per Frequency')
    %ylim([0 1]) 
    
    suptitle(['Total sparseness = ' num2str(sp)])
end