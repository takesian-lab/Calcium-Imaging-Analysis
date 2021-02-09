function ops2 = get_abridged_ops(ops)
% This function extracts the user-defined variables in Fall.ops (defined
% in the suite2p ops window) to prepare them to be stored in block. This is
% helpful because otherwise Fall.ops is too big to save.
% 
% Argument(s): 
%   Fall.ops
% 
% Returns:
%   ops2 (an abridged version of Fall.ops)
% 
% Notes:

%% File paths
ops2.look_one_level_down = ops.look_one_level_down;
ops2.save_path0 = ops.save_path0;
ops2.fast_disk = ops.fast_disk;
ops2.bruker = ops.bruker;

%% Main settings
ops2.nplanes = ops.nplanes;
ops2.nchannels = ops.nchannels;
ops2.functional_chan = ops.functional_chan;
ops2.tau = ops.tau;
ops2.fs = ops.fs;
ops2.do_bidiphase = ops.do_bidiphase;
ops2.bidiphase = ops.bidiphase;

%% Output settings
ops2.preclassify = ops.preclassify;
ops2.save_mat = ops.save_mat;
ops2.combined = ops.combined;
ops2.reg_tif = ops.reg_tif;
ops2.reg_tif_chan2 = ops.reg_tif_chan2;
ops2.aspect = ops.aspect;
ops2.delete_bin = ops.delete_bin;
%ops2.move_bin = ops.move_bin;

%% Registration
ops2.do_registration = ops.do_registration;
ops2.align_by_chan = ops.align_by_chan;
ops2.nimg_init = ops.nimg_init;
ops2.batch_size = ops.batch_size;
ops2.smooth_sigma = ops.smooth_sigma;
% ops2.smooth_sigma_time = ops.smooth_sigma_time;
ops2.maxregshift = ops.maxregshift;
ops2.th_badframes = ops.th_badframes;
ops2.keep_movie_raw = ops.keep_movie_raw;
% ops2.two_step_registration = ops.two_step_registration;

%% Nonrigid
ops2.nonrigid = ops.nonrigid;
ops2.block_size = ops.block_size;
ops2.snr_thresh = ops.snr_thresh;
ops2.maxregshiftNR = ops.maxregshiftNR;

%% 1P
ops2.spatial_hp = ops.spatial_hp;
ops2.pre_smooth = ops.pre_smooth;
ops2.spatial_taper = ops.spatial_taper;

%% ROI detection
ops2.roidetect = ops.roidetect;
ops2.sparse_mode = ops.sparse_mode;
ops2.diameter = ops.diameter;
ops2.spatial_scale = ops.spatial_scale;
ops2.connected = ops.connected;
ops2.threshold_scaling = ops.threshold_scaling;
ops2.max_overlap = ops.max_overlap;
ops2.max_iterations = ops.max_iterations;
ops2.high_pass = ops.high_pass;

%% Extraction\Neuropil
ops2.allow_overlap = ops.allow_overlap;
ops2.inner_neuropil_radius = ops.inner_neuropil_radius;
ops2.min_neuropil_pixels = ops.min_neuropil_pixels;

%% Deconvolution
ops2.win_baseline = ops.win_baseline;
ops2.sig_baseline = ops.sig_baseline;
ops2.neucoeff = ops.neucoeff;

%% Zcorr
if isfield(ops,'zcorr')
    ops2.zcorr = ops.zcorr;
end
end

