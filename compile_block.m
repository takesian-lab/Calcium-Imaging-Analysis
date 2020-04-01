% Pull out and compile:
%   - Tosca data (locomotion, trial types)
%   - Bruker-derived timestamps
%   - 2p data from Fall.mat

function block = compile_block(setup)

    disp('Processing...');
    disp(setup.block_name);

    block = struct;
    block.setup = setup;

    %% Behavior, locomotion, sound, and Suite2p data

    %pull out the Tosca-derived, behaviorally relevant data
    [block] = behavior_RF_singleblock(block);

    %pull out the Bruker-derived timestamps from BOTs and Voltage Recordings
    [block] = define_sound_singleblock(block);

    %determine which trials are considered "active (locomotor)"
    % This might not be necessary to do here, but leaving in for now.
    [block] = define_loco_singleblock(block);
    % [loco_activity,isLocoSound,data] = isLoco(setup,data);

    %pull out block-specific data from Fall.mat
    [block] = define_suite2p_singleblock(block);

end