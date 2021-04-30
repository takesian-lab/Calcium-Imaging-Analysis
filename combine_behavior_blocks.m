function combinedBlock = combine_behavior_blocks(block1, block2)
%Maryse March 2021
   
combinedBlock = struct;
combinedBlock.setup = block1.setup;
if ~isfield(block1, 'setup2')
    combinedBlock.setup2 = block2.setup;
elseif ~isfield(block1, 'setup3')
    combinedBlock.setup2 = block1.setup2
    combinedBlock.setup3 = block2.setup;
else
    error('More than 2 blocks found, combine manually')
end
combinedBlock.Combined = true;
combinedBlock.start_time = [block1.start_time, block2.start_time];
combinedBlock.Tosca_times = [block1.Tosca_times, block2.Tosca_times];
combinedBlock.errors = [block1.errors, block2.errors];
combinedBlock.New_sound_times = [block1.New_sound_times, block2.New_sound_times];
combinedBlock.New_sound_idx = [block1.New_sound_idx, block2.New_sound_idx];
combinedBlock.lick_time = [block1.lick_time; block2.lick_time];
combinedBlock.concat_times = [block1.concat_times; block2.concat_times + block1.concat_times(end)];
combinedBlock.concat_licks = [block1.concat_licks; block2.concat_licks];
combinedBlock.Outcome = [block1.Outcome, block2.Outcome];
combinedBlock.trialType = [block1.trialType, block2.trialType];
combinedBlock.TargetFreq = [block1.TargetFreq, block2.TargetFreq];
combinedBlock.parameters.variable1 = [block1.parameters.variable1, block2.parameters.variable1];
combinedBlock.parameters.variable2 = [block1.parameters.variable2, block2.parameters.variable2];
%skip combinedBlock.loco_data
combinedBlock.loco_activity = [block1.loco_activity; block2.loco_activity];
combinedBlock.loco_times = [block1.loco_times; block2.loco_times + block1.loco_times(end)];
combinedBlock.loc_Trial_times = [block1.loc_Trial_times, block2.loc_Trial_times];
combinedBlock.loc_Trial_activity = [block1.loc_Trial_activity, block2.loc_Trial_activity];
combinedBlock.rxn_time = [block1.rxn_time, block2.rxn_time];
combinedBlock.locIDX = [block1.locIDX, block2.locIDX];
combinedBlock.holdingPeriod = [block1.holdingPeriod, block2.holdingPeriod];
combinedBlock.stim_level = unique([block1.stim_level, block2.stim_level]);
