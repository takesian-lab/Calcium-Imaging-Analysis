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


    tempFrame_set = {};
    for k = 1:size(currentInfo_R,1)
        path = strcat(currentInfo_R{k,P}, '/', currentInfo_R{k,M}, '/', currentInfo_R{k,AP});
        cd(path);
        load('Fall.mat', 'ops');
        Imaging_Block = currentInfo_R{k,IS};
        showTable = 0;
        if k == 1
            display(currentMouse);
            showTable = 1;
        end
        tempFrame_set{1,k} = get_frames_from_Fall(ops,Imaging_Block,showTable);
    end
    
        setup.Frame_set         =   tempFrame_set;


    block = struct;
    block.setup = setup;

end