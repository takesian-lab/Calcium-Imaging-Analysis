% List of ops to check

% This script saves an ops.mat file with a list of Suite2p ops that the 
% user would like to check against their compiled blocks. This check will
% happen during define_suite2p_singleblock.

% Each user should set their own ops file and save it in the same folder as Info

% Steps:
% 1. Save script as ops_user
% 2. Modify filepath and filename
% 3. Modify list of ops [comment those you don't want to check]
% 4. Run script

% The list of ops is based on the Suite2p GitHub version from May 2020

%% Save filepath

filepath = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
filename = 'Maryse_ops2.mat';

%% File paths
%ops.look_one_level_down
%ops.save_path0
%ops.fast_disk
ops.bruker = 0;

%% Main settings
ops.nplanes = 1;
ops.nchannels = 1;
ops.functional_chan = 1;
ops.tau = 0.7;
ops.fs = 30;
ops.do_bidiphase = 1;
ops.bidiphase = 0;

%% Output settings
ops.preclassify = 0;
%ops.save_mat = 1;
ops.combined = 1;
%ops.reg_tif
%ops.reg_tif_chan2
ops.aspect = 1;
%ops.delete_bin
%ops.move_bin

%% Registration
ops.do_registration = 1;
ops.align_by_chan = 1;
ops.nimg_init = 300;
ops.batch_size = 500;
ops.smooth_sigma = 1.15;
ops.smooth_sigma_time = 0;
ops.maxregshift = 0.1;
ops.th_badframes = 1;
%ops.keep_movie_raw
%ops.two_step_registration

%% Nonrigid
ops.nonrigid = 1;
ops.block_size = [128, 128];
ops.snr_thresh = 1.2;
ops.maxregshiftNR = 5.0;

%% 1P
ops.spatial_hp = 50;
ops.pre_smooth = 2;
ops.spatial_taper = 50;

%% ROI detection
ops.roidetect = 1;
ops.sparse_mode = 1;
%ops.diameter = 10;
ops.spatial_scale = 0;
ops.connected = 1;
ops.threshold_scaling = 1;
ops.max_overlap = 0.75;
ops.max_iterations = 20;
ops.high_pass = 100;

%% Extraction\Neuropil
ops.allow_overlap = 0;
ops.inner_neuropil_radius = 3;
ops.min_neuropil_pixels = 350;

%% Deconvolution
ops.spikedetect = 1;
ops.win_baseline = 60;
ops.sig_baseline = 10;
ops.neucoeff = 0.7;

%% Save ops
cd(filepath)
save(filename, 'ops', '-mat');

