function    [data]=concatBlocks_aligned(data);
%
% DOCUMENTATION IN PROGRESS
%
% What does this function do?
% concatBlocks takes individual blocks (generated in
% compile_blocks_from_info) and concatenates them based on stim-type. This
% funcion is for the 'aligned_to_stim' data.
%
% Argument(s):
%   data(struct)
%
% Returns:
%   data(struct)
%
% Notes:
%
%
% TODO: Magic numbers, does this work with multiple ROIs?
% Search 'TODO'

setup = data.setup;

for a=1:size(setup.mousename,1) %Mice
    for b=1:size(setup.mousename,2) %ROIs
        
        if isempty(setup.mousename{a,b})
            continue;
        end
        
        mouseID=setup.mousename{a,b};
        Imaging_Block=setup.Imaging_sets{a,b};
 
        for i=1:length(Imaging_Block)
            
            unique_block_name = setup.unique_block_names{a,b}(i);
            block = data.([mouseID]).([unique_block_name]);
            
            if ismissing(block.setup.suite2p_path)
                disp('Skipping Suite2p data for block...');
                disp(unique_block_name);
                return
            end
            
            if setup.run_redcell ==1
              data.([mouseID]).parameters.red_idx  = block.redcell;
              data.([mouseID]).parameters.green_idx = ~block.redcell;
            end
            
            
%             timestamp =  block.timestamp;
%             Sound_Time = block.Sound_Time;
            isLoco = block.isLoco;
            
            
            % correct for the third dimension of data.
            if i == 1
                F_cat = block.aligned_stim.F_stim;
                F7_cat = block.aligned_stim.F7_stim;
                spks_cat = block.aligned_stim.spks_stim;
                neu_cat = block.aligned_stim.Fneu_stim;
                V1_cat = block.parameters.variable1;
                V2_cat = block.parameters.variable2;
                isLoco_cat = block.isLoco;
            else i>1
                % check to make sure that the size (frames) of F, F7,
                % spks,and neu is the same across blocks. The third dim of
                % each of these matricies should be the same, so I only
                % used F_cat to test for size issues:
                if size(F_cat,3) == size(block.aligned_stim.F_stim,3)
                    F_cat = cat(2,F_cat,block.aligned_stim.F_stim);
                    F7_cat = cat(2,F7_cat,block.aligned_stim.F7_stim);
                    spks_cat = cat(2,spks_cat,block.aligned_stim.spks_stim);
                    neu_cat = cat(2,neu_cat,block.aligned_stim.Fneu_stim);
                    
                elseif size(F_cat,3) > size(block.aligned_stim.F_stim)
                    % make new array dim3 same size as F_cat dim3
                    diff_size = size(F_cat,3) - size(block.aligned_stim.F_stim,3);
                    add_size = NaN(size(F_cat,1),size(F_cat,2),size(diff_size));
                    
                    % pad the arrays with Nans
                    F_nan = cat(3,block.aligned_stim.F_stim,add_size);
                    F7_nan =  cat(3,block.aligned_stim.F7_stim,add_size);
                    spks_nan = cat(3,block.aligned_stim.spks_stim,add_size);
                    neu_nan =  cat(3,block.aligned_stim.Fneu_stim,add_size);
                    
                    % concat the arrays
                    F_cat = cat(2, F_cat,F_nan);
                    F7_cat = cat(2,F7_cat,F7_nan);
                    spks_cat = cat(2,spks_cat,spks_nan);
                    neu_cat = cat(2,neu_cat,neu_nan);
                    
                else size(F_cat,3) < size(block.aligned_stim.F_stim)
                    %make dim3 larger for F_cat...
                    diff_size = size(block.aligned_stim.F_stim,3)- size(F_cat,3); 
                    add_size = NaN(size(F_cat,1),size(F_cat,2),size(diff_size));
                    F_cat = cat(3,F_cat,add_size);
                    F7_cat =  cat(3,F7_cat,add_size);
                    spks_cat = cat(3,spks_cat,add_size);
                    neu_cat =  cat(3,neu_cat,add_size);
                    
                    % then concat the arrays...
                    F_cat = cat(2, F_cat,block.aligned_stim.F_stim);
                    F7_cat = cat(2,F7_cat,block.aligned_stim.F7_stim);
                    spks_cat = cat(2,spks_cat,block.aligned_stim.spks_stim);
                    neu_cat = cat(2,neu_cat,block.aligned_stim.Fneu_stim);
                end
                % cat is Loco and variables...
                V1_cat = cat(2,V1_cat,block.parameters.variable1);
                V2_cat = cat(2,V2_cat,block.parameters.variable2);
                isLoco_cat = cat(2,isLoco_cat,block.isLoco);
            end     
        end
    end

    data.([mouseID]).cat.F_cat = F_cat;
    data.([mouseID]).cat.F7_cat = F7_cat;
    data.([mouseID]).cat.spks_cat = spks_cat;
    data.([mouseID]).cat.neu_cat = neu_cat;
    data.([mouseID]).cat.V1_cat = V1_cat;
    data.([mouseID]).cat.V2_cat = V2_cat;
    data.([mouseID]).cat.isLoco_cat = isLoco_cat;
    

end
end

