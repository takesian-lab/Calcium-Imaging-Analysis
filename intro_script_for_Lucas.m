%Introduction to Takesian Lab Calcium Imaging Analysis

%Load a compiled block (this one is FM sweep)
cd('\\apollo\research\ENT\Takesian Lab\Maryse\2p analysis\CompiledBlocks_VIPvsNDNF')
load('Compiled_NxDD070420F2_FOV1_2020-09-16_Session01_Run01_Block001_FM.mat')

%Preview the locomotor activity, raw data, and FOV
visualize_block(block)

%This script automatically pulls out sound-responsive (either positively or negatively
%responding) cells and plots their average activity. This code is still very much a
%work in progress, so we are still working on the best automatic classification of active cells
visualize_active_cells(block)

%Look at all cells regardless of responsiveness (warning, this makes a lot of graphs)
visualize_cell(block, block.cell_number)

%For a synchronization analysis, you'll probably want to use spontaneous data.
%For these we record ~2 minutes of spontaneous activity (not in silence because the
%2p rig makes ~60dB noise with the equipment on) with 20 sham trial markers embedded
load('Compiled_NxDD070420F2_FOV1_2020-09-16_Session01_Run04_Block004_spontaneous.mat')

%Compute correlations between cells using either df/f or deconvolved spikes:

%Pull out raw signal and correct for neuropil
F = block.F; %Raw fluorescence trace for all cells (n cells * n frames)
Fneu = block.Fneu; %Neuropil trace for all cells
F7 = F - block.setup.constant.neucoeff*Fneu; %Neuropil-corrected trace

%Compute df/f
mean_F = mean(F7,2);
df_f = (F7 - mean_F)./mean_F; %(total-mean)/mean

%Deconvolved spike train
spks = block.spks;

%INSERT AMAZING SYNCHRONIZATION COMPUTATIONS HERE :D
neuron1 = df_f(1,:);
neuron2 = df_f(2,:);
framerate = block.setup.framerate; %frames per second
maxlag = 2*framerate; %Max lag to compute in frames for xcorr, this is 2 seconds

%Compute cross-correlogram
[cc, lags] = xcorr(neuron1, neuron2, maxlag, 'coeff'); %Are these the right parameters for 2p data?

figure
plot(lags, cc)
ylabel('Cross correlation')
xlabel('Time lags (frames)')

%Variables to pull out from cc for each neuron pair:
% CC peak
% CC lag (time at peak)
% CC width @ half of peak
% Z-score (statistical measure of synchronization to say whether pair is
% significantly correlated or not)

%REFERENCES:
%- Eggermont JJ. 1992. Neural interaction in cat primary auditory cortex.
%Dependence on recording depth, electrode separation, and age. J Neurophysiol. 68:1216–1228.
%- Brosch M, Schreiner CE. 1999. Correlations between neural discharges are related to receptive field
%properties in cat primary auditory cortex. Eur J Neurosci. 11:3517–3530.

%CAVEATS & CONSIDERATIONS:
%   - How do you decide which cells to include? Although we manually curate
%   all Suite2p ROIs to establish our final set of neurons, it is possible
%   that some ROIs are not actually neurons, or the cells are out of plane
%   for that block, or the cells are not spontaneously active. We have
%   conversations about this with Anne all the time, and still don't have
%   the perfect answer

%   - Is the brain activity synchronized with locomotor activity? If so, is
%   this just a motion artifact or is it true synchronization? You may want
%   to exclude locomotor bouts or analyze them separately
    activity = block.loco_activity;
    activity(activity < block.setup.constant.locoThresh) = 0; %correct for baseline noise
    inactive_time = activity == 0;
    
%   For SOME fields of view, we recorded Z-stacks and were able to compute
%   the degree of brain movement and correlate this with loco and neural
%   activity, see here:
    visualize_zcorr(block)
    
%   - Degree of synchronization between cells is related to both distance and CF
%   You'll definitely want to pull out information about distance between
%   cells (we haven't done this yet) and CF would be a bonus.
%   Info about the cells' X/Y position is located in block.stat

%   To get the CF data, you could use a Receptive Field block recorded in
%   the same field of view and check out extract_data_v3 for how
%   to pull out CF for each cell. However, in some cases matching cells
%   across blocks is not trivial (but Carolyn and I are working on this)
