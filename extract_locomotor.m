% use extractData to analyze motor activity
% CGS 10/19/20
%%
PC_name = getenv('computername');
'RD0332' %Carolyn
        cellList_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020';
        blocks_path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\CompiledBlocks';
        save_path = 'Z:\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS\Inactive';
        cellList_filename = 'Responsive cells v2';

%% 
dataType = 'FM';
MouseType = 'VIP';

%%
cd(cellList_path)
cellList = importfile(cellList_filename);
cellList(1,:) = []; %Remove header

cd(blocks_path)

%Look at only one stim type at a time
dataTypes = [cellList{:,4}]';
MouseTypes = [cellList{:,1}]';


if ~isempty(dataType)
    cellList = cellList(strcmpi(dataTypes, dataType),:);
end


%Loop through all cells in each block
blocks = [cellList{:,5}]';
uniqueBlocks = unique(blocks);

for b = 1:length(uniqueBlocks)
    currentBlock = uniqueBlocks{b};
    block_cellList = cellList(strcmpi(blocks, currentBlock),:);
    load(currentBlock);
    [loco_cor_cell] = visualize_locomotor(block);
    DataLoco.loco_cor_cell{b} = loco_cor_cell;
    
end
