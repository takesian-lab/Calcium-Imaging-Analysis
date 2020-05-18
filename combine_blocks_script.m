% combine_blocks_script
% user-friendly way to run combine_blocks

%% For now you can only combine 2 blocks. Later we could make this longer.

block_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocksDeleteLater';
cd(block_path);

disp('Load first block');
[FileName,PathName] = uigetfile('.mat');
load([PathName '/' FileName])
block1 = block;

disp('Load second block');
[FileName,PathName] = uigetfile('.mat');
load([PathName '/' FileName])
block2 = block;

[block, filename] = combine_blocks(block1, block2);

save(filename, 'block');