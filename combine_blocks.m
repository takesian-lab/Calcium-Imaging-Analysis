function [combined_block, filename] = combine_blocks(block1, block2)
% DOCUMENTATION IN PROGRESS
% 
% This function allows you to combine the data from two blocks, given that
% they have the same ROIs. This can be used if you recorded two blocks back
% to back, for example
% 
% Argument(s): 
%   block2 and block2 (structs)
% 
% Returns:
%   combined_block
% 
% Notes:
%
%
% TODO: Verify that data are combined in a way that we can still use them.
% For example, in some cases it might not sense to concatenate without
% changing the time information (usu. first column)
% Aligned to stim data should be constrained so that they have the same z
% dimension (n frames)
% Search 'TODO'

%% Create new block struct
combined_block = struct;
combined_block.files = {block1.setup.block_filename, block2.setup.block_filename};
combined_block.setup = block1.setup;
combined_block.setup2 = block2.setup;
filename = strcat('Combined_', block1.setup.block_filename, '.mat');

%% Combine Tosca data
combined_block = hcat(combined_block, block1, block2, 'New_sound_times');
combined_block = hcat(combined_block, block1, block2, 'start_time');
combined_block = vcat(combined_block, block1, block2, 'lick_time');
combined_block = hcat(combined_block, block1, block2, 'Tosca_times');
combined_block = hcat(combined_block, block1, block2, 'Outcome');
combined_block = hcat(combined_block, block1, block2, 'trialType');
combined_block = check_equal(combined_block, block1, block2, 'TargetFreq');
parameters = struct;
combined_block.parameters = hcat(parameters, block1.parameters, block2.parameters, 'variable1');
combined_block.parameters = hcat(combined_block.parameters, block1.parameters, block2.parameters, 'variable2');
combined_block = vcat(combined_block, block1, block2, 'loco_data');
combined_block = vcat(combined_block, block1, block2, 'loco_activity');
combined_block = vcat(combined_block, block1, block2, 'loco_times');
combined_block = vcat(combined_block, block1, block2, 'locomotion_data');
combined_block = vcat(combined_block, block1, block2, 'active_time');
combined_block = hcat(combined_block, block1, block2, 'Sound_Time');
combined_block = vcat(combined_block, block1, block2, 'locTime');

%% Combine Bruker data
combined_block = vcat(combined_block, block1, block2, 'timestamp');
combined_block = hcat(combined_block, block1, block2, 'isLoco');

%% Combine Suite2p data
%combined_block.ops = check_equal(combined_block, block1, block2, 'ops'); %Does not need to be equal
combined_block = check_equal_struct(combined_block, block1, block2, 'img');
combined_block = check_equal_struct(combined_block, block1, block2, 'iscell');
combined_block = check_equal_struct(combined_block, block1, block2, 'cell_number');
combined_block = check_equal_struct(combined_block, block1, block2, 'stat');
combined_block = hcat(combined_block, block1, block2, 'F');
combined_block = hcat(combined_block, block1, block2, 'Fneu');
combined_block = hcat(combined_block, block1, block2, 'spks');
combined_block = check_equal_struct(combined_block, block1, block2, 'redcell');

%% Combine aligned to stim data
try
    combined_block = vcat(combined_block, block1, block2, 'aligned_to_stim');
catch
    warning('aligned_to_stim mats do not match on all dimensions. Not included in combined data')
end

end

%% Helper functions

function combined_block = vcat(combined_block, block1, block2, fieldname)

if isfield(block1, fieldname) && isfield(block2, fieldname)
    combined_block.(fieldname) = [block1.(fieldname); block2.(fieldname)];
end

end

function combined_block = hcat(combined_block, block1, block2, fieldname)

if isfield(block1, fieldname) && isfield(block2, fieldname)
    combined_block.(fieldname) = [block1.(fieldname), block2.(fieldname)];
end

end

function combined_block = check_equal(combined_block, block1, block2, fieldname)

if isfield(block1, fieldname) && isfield(block2, fieldname)
    if isnan(block1.(fieldname)) && isnan(block2.(fieldname))
        combined_block.(fieldname) = nan;
    elseif isnan(block1.(fieldname)) || isnan(block2.(fieldname))
        error('Blocks do not match')
    elseif ~isequal(block1.(fieldname), block2.(fieldname))
        error('Blocks do not match')
    else %blocks match
        combined_block.(fieldname) = block1.(fieldname);
    end
end

end


function combined_block = check_equal_struct(combined_block, block1, block2, fieldname)

if isfield(block1, fieldname) && isfield(block2, fieldname)
    if ~isequal(block1.(fieldname), block2.(fieldname))
        error('Blocks do not match')
    else %blocks match
        combined_block.(fieldname) = block1.(fieldname);
    end
end

end