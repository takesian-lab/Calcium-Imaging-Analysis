function [FRA, data, fig1, fig2, fig3, fig4, fig5] = compute_frequency_tuning_v2(responsive, RF, plot_tuning, freqs, ints, cell)
% Compute tuning properties of a single cell and plot figures
%
% Argument(s): 
%   responsive (X x Y matrix) - matrix of 0s and 1s indicating whether the
%   cell is significantly responsive to that stim combination
%   RF (X x Y matrix) - matrix of cell's response to that stim combination
%   plot_tuning - 0 or 1 to plot figures
%   freqs (X vector) - list of frequencies from low to high corresponding to RF stim
%   ints (Y vector - list of intensities from high to low corresponding to RF stim
%   cell (int) - cell number
%
% Returns:
%   FRA
%   data (1x12 vector) containing:
%     - BF   = frequency with the strongest response
%     - BF_I = intensity corresponding to BF
%     - CF   = frequency with the strongest response at threshold
%     - CF_I = threshold/intensity corresponding to CF
%     - BestInt = (measured at BF) best intensity/peak of the gaussian intensity tuning curve
%     - Int_RMS = (measured at BF) RMS width of the gaussian intensity tuning curve at half the maximum
%     - Int_halfwidth = (measured at BF) width of the gaussian intensity tuning curve at half the maximum 
%     - BW_20_RMS = freq bandwidth 20dB above threshold (RMS)
%     - BW_20_halfwidth = freq bandwidth 20dB above threshold (width)
%     - BW_BF_RMS = freq bandwith at BF (RMS)
%     - BW_BF_halfwidth = freq bandwith at BF (width)
%     - ISI  = Intensity Selectivity Index (from Mesik et al)
%     - type = 1. simple (1 peak) or 2. complex (multipeaked) RF
%   fig1 - Intensity tuning
%   fig2 - CF
%   fig3 - BW20
%   fig4 - Responses at BF
%   fig5 - FRA
% 
% Notes:
%
%
% TODO: 
% Search 'TODO'
%% Initial parameters

%disp(['Computing frequency tuning for cell ' num2str(cell)]);

%In case RF is not 8x8
if numel(RF) > length(freqs)*length(ints)
    RF = RF(1:length(ints), 1:length(freqs));
    responsive = responsive(1:length(ints), 1:length(freqs));
end

%Check RF for NaNs:
% if any(isnan(RF(:)))
%     warning('RF contains NaNs, fits will remove to run')
%     %RF will contain NaNs for conditions where there were no trials (for
%     %example due to loco activity)
% end

%Figures will return as [] if plot_tuning is 0
[fig1, fig2, fig3, fig4, fig5] = deal([]);

%Variables will return as NaN if fits cannot be made
[BF, BF_I, CF, CF_I, BestInt, Int_RMS, Int_halfwidth, BW_20_RMS, BW_20_halfwidth, BW_BF_RMS, BW_BF_halfwidth, ISI, typeNumber] = deal(nan);

%convert freqs to octaves 
octaves = log2(freqs./freqs(1));
x_freq = octaves(1):0.1:octaves(length(octaves));

%flip RF such that lowest intensities on top
if ~issorted(ints,'descend')
    error('Intensities should be descending')
end
flip_RF = flip(RF,1);
flip_ints = flip(ints);
flip_responsive = flip(responsive,1);


%% Find threshold and BF

%find CF_I (threshold)
[row,col] = find(flip_responsive==1);
CF_row = min(row);
CF_I = flip_ints(CF_row);

%find BF
[M, I] = max(flip_RF(:));
[BF_row, BF_col] = ind2sub(size(flip_RF),I);
BF_I = flip_ints(BF_row);
BF = freqs(BF_col);

%% Compute intensity tuning (ISI)

%applies Mesik et al approach to compute intensity tuning
%'To generate an intensity tuning curve, spike rates at the CF and two 
%neighboring frequencies were averaged for each tone intensity. Intensity 
%selectivity index (ISI) was calculated as 1 minus the ratio between the 
%spike count at 30 dB above the best intensity (i.e., the intensity that 
%produced maximum spike count) or the highest intensity tested and 
%that at the best intensity.'

%1 - (response at highest intensity tested/response at BF_I)
%Result will be between 0 and 1: intensity selective = 1, monotonic = 0

%normalize RF_F 
RF_norm = RF./M;
flip_RF_norm = flip(RF_norm);

best = flip_RF_norm(BF_row,BF_col);
if best <0; best = 0; end
highest = RF_norm(1,BF_col);

if highest <0; highest = 0; end
ISI = 1-(highest/best); 

%% Fit #1 - Intensity Tuning
% determine intensity tuning using peak of gaussian fit across intensities at BF
options = fitoptions('gauss1');
options.Lower = [0 1 0];
options.Upper = [inf 80 80];
x = flip_ints;
y = flip_RF(:,BF_col);

if length(y(~isnan(y))) >= 3
    gauss_fit = fit(x(~isnan(y)), y(~isnan(y)), 'gauss1',options);
    x_ints = flip_ints(1):0.1:flip_ints(length(flip_ints));
    gauss_curve = gauss_fit(x_ints);

    Int_RMS = gauss_fit.c1*2; %gauss_RMS_width
    Int_halfwidth = gauss_fit.c1*2.3548; %full width of the gaussian curve at half the maximum 

    [BI_A,BI] = max(gauss_fit(x_ints)); 
    BestInt = x_ints(BI);
    
    if plot_tuning
        fig1 = figure; %plot gaussian fits
        plot(flip_ints,flip_RF(:,BF_col));
        hold on; plot(x_ints, gauss_curve,'r');
        sz_half_BW_int = Int_halfwidth./2;
        y = [BI_A/2 BI_A/2];
        x = [BestInt-(sz_half_BW_int) BestInt+(sz_half_BW_int)];
        hold on; plot(x,y,'c');
        title(['Intensity Tuning Curve for Cell = ',num2str(cell)]);
        xlabel('Intensities (dB SPL)');
        ylabel('Amplitude (df/f)');
    end
else
    warning('Not enough points to fit intensity tuning')
end

%% Fit #2 - CF
%determine CF using peak of gaussian fit at threshold intensity
options = fitoptions('gauss1');
options.Lower = [0 1 0];
options.Upper = [inf 4 2];
x = octaves;
y = flip_RF(CF_row,:)';

if length(y(~isnan(y))) >= 3
    gauss_fit = fit(x(~isnan(y)), y(~isnan(y)), 'gauss1',options);
    gauss_curve = gauss_fit(x_freq);
    [CF_A,CF_F] = max(gauss_fit(x_freq)); 

    if plot_tuning
        fig2 = figure; %plot gaussian fits
        plot(octaves,flip_RF(CF_row,:));
        hold on; plot(x_freq, gauss_curve,'r');
        title(['Threshold Responses for Cell =  ',num2str(cell)]);
        xlabel('Frequencies (octaves)');
        ylabel('Amplitude (df/f)');
    end
    CF = x_freq(CF_F);
else
    warning('Not enough points to fit CF')
end

%% Fit #3 - BW20 
%find bandwidth 20dB above threshold
if CF_row + 2 < length(ints) %Check if BW20 exists
    options = fitoptions('gauss1');
    options.Lower = [0 1 0];
    options.Upper = [inf 4 2];
    x = octaves;
    y = flip_RF(CF_row+2,:)';
    
    if length(y(~isnan(y))) >= 3
        gauss_fit = fit(x(~isnan(y)), y(~isnan(y)), 'gauss1',options);
        [A,F] = max(gauss_fit(x_freq));
        F = x_freq(F);
        gauss_curve = gauss_fit(x_freq);
        BW_20_RMS = gauss_fit.c1*2; %gauss_RMS_width
        BW_20_halfwidth = gauss_fit.c1*2.3548; %full width of the gaussian curve at half the maximum 

        if plot_tuning
             fig3 = figure; %plot gaussian fits
             plot(octaves,flip_RF(CF_row+2,:));
             hold on; plot(x_freq, gauss_curve,'r');
             sz_half_BW = BW_20_halfwidth./2;
             y = [A/2 A/2];
             x = [F-(sz_half_BW) F+(sz_half_BW)];
             hold on; plot(x,y,'c');
             title(['20dB above Threshold Responses for Cell =  ',num2str(cell)]);
             xlabel('Frequencies (octaves)');
             ylabel('Amplitude (df/f)');
        end
    else
       warning('Not enough points to fit BW20')
    end   
end

%% Fit #4 - bandwidth at BF
%find bandwidth at BF
options = fitoptions('gauss1');
options.Lower = [0 1 0];
options.Upper = [inf 4 2];
x = octaves;
y = flip_RF(BF_row,:)';

if length(y(~isnan(y))) >= 3
    gauss_fit = fit(x(~isnan(y)), y(~isnan(y)), 'gauss1',options);

    [A2,F2] = max(gauss_fit(x_freq));
    gauss_curve = gauss_fit(x_freq);
    BW_BF_RMS = gauss_fit.c1*2; %gauss_RMS_width
    BW_BF_halfwidth = gauss_fit.c1*2.3548; %full width of the gaussian curve at half the maximum 

    F2 = x_freq(F2);

    if plot_tuning
         fig4 = figure; %plot gaussian fits
         plot(octaves,flip_RF(BF_row,:));
         hold on; plot(x_freq, gauss_curve,'r');
         sz_half_BW = BW_BF_halfwidth./2;
         y = [A2/2 A2/2];
         x = [F2-(sz_half_BW) F2+(sz_half_BW)];
         hold on; plot(x,y,'c');
         title(['Responses at BF for Cell =  ',num2str(cell)]);
         xlabel('Frequencies (octaves)');
         ylabel('Amplitude (df/f)');
    end
else
   warning('Not enough points to fit bandwidth at BF')
end    

%% Fit # 5 - Gaussian at highest intensity (not used?)
%fit gaussian at highest intensity           
options = fitoptions('gauss1');
options.Lower = [0 1 0];
options.Upper = [inf 4 2];
x = octaves;
y = flip_RF(size(RF,1),:)';

if length(y(~isnan(y))) >= 3
    gauss_fit = fit(x(~isnan(y)), y(~isnan(y)), 'gauss1',options);
    [A3,F3] = max(gauss_fit(x_freq));
else
   warning('Not enough points to fit gaussian at highest intensity')
end   

%% FRA - Frequency Receptive Area
%make frequency receptive area (FRA) contour
%determine type, single peak if continuous, dual peak, or comples (multi-peak)   
[FRA,L] = bwboundaries(responsive,'noholes');

if length(FRA)==1
    type = 'single peak';
    typeNumber = 1;
else
    type = 'complex peak';
    typeNumber = 2;
end

if plot_tuning
    fig5 = figure; %plot FRA contour
    h=imagesc(RF);
    %imshow(label2rgb(L, @jet, [.5 .5 .5]))
    hold on;
    %FRA_flip = flip(FRA,1); % flip back for normal FRA
    FRA_white = ones(size(L))*0.3;
    FRA_white(L>0) = 1;
     for k = 1:length(FRA)
         boundary = FRA{k};
         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth',2);
         hold on;
     end
     set(h, 'AlphaData', FRA_white);

     ylabel('dB');
     xlabel('Frequencies');
     set(gca,'XTickLabel',freqs);
     set(gca, 'YTickLabel', ints);
     title(['Frequency Receptive Areas for cell =  ',num2str(cell)]);
end

%% Store Data
data = [BF, BF_I, CF, CF_I, BestInt, Int_RMS, Int_halfwidth, BW_20_RMS, BW_20_halfwidth, BW_BF_RMS, BW_BF_halfwidth, ISI, typeNumber];

end %end function