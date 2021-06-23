# Calcium-Imaging-Analysis

Authors: Anne Takesian, Carolyn Sweeney, Maryse Thomas, and Wisam Reid

##### This is a work in progress


The stim_protocol code is:
- noiseburst      = 1
- ReceptiveField  = 2
- FM sweep        = 3
- widefield       = 4
- SAM             = 5
- SAM freq        = 6
- behavior        = 7
- behavior        = 8
- random h20      = 9
- noiseburst_ITI  = 10
- random air      = 11
- spontaneous     = 12
- behavior MT     = 13

Output block structure:

block.setup 		= list of variables created from Info spreadsheet
block.setup.constant 	= user-defined constants from the top of compile_blocks_from_info script
block.start_time 	= raw Tosca timestamp for the start of each trial
block.Tosca_times 	= raw Tosca timestamps for every sample in each trial
block.errors 		= list of error trials (which have been removed from every variable saved in block)
block.New_sound_times 	= sound time in seconds relative to each trial start
block.lick_time		= lick raster (0 or 1) for every sample in each trial
block.concat_times 	= concatenated trials in seconds (time-corrected)
block.concat_licks  	= concatenated lick raster corresponding to time in continuous_trials
block.outcome		= CS+/CS- outcomes: 1 = hit, 0 = miss, 3 = withhold, 4 = false alarm; nan = other
block.trialType		= CS+/CS- identities: 1 = CS+, 0 = CS-
block.TargetFreq	= Stim frequency for behavior blocks
block.parameters	= stimulus parameters for every trial (e.g. freq. and dB)
block.loco_data		= raw loco data if block does not have Bruker data
block.loco_activity	= concatenated locomotor activity trace for the full block (time-corrected)
block.loco_times	= timestamp in seconds corresponding to loco_activity
block.loc_Trial_times 	= locomotor activity for each trial (time-corrected)
block.loc_Trial_activity= timestamp in seconds corresponding to loc_Trial_times
block.rxn_time 		= reaction time in ms
block.locIDX		= indices of loco_activity/loco_times that correspond to each trial

-----

# Running the code
## Command line

1. Compile blocks
- Make Info excel spreadsheet outside of MATLAB. Look at the template spreadsheet and "Data Structure for Info" PDF within this folder for help.
- Update paths at the top of compile_blocks_from_info_v2
- Run compile_blocks_from_info_v2
- Visualize the output of a single block with the function visualize_block
- Visualize the output of a single cell from a block with the function visualize_cell or visualize_active_cells

2. Analyze
- You can use either the same Info spreadsheet as before or create a dedicated spreadsheet for each type of analysis you want to do.
- Update paths at the top of analysis_2P_v3
- Run analysis_2P_v3

# Code Structure

-----

# Development

Search **TODO** in the source code files


| TODO: | Code Location              | Call Location              | Task                                                             |
|:-----:|:--------------------------:|:--------------------------:|:----------------------------------------------------------------:|
| 1.    | [example.m]      | [main.m] | add a feature                                 |


## Miscellaneous Notes

# Dependencies
