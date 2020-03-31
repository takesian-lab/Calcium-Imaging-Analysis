function block = compile_block(setup, 2p_path, Tosca_path)

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
