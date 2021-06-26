function convert_bidirectional_tiffs(path)
% convert_bidirectional_tiffs

% Convert piezo z-stack recordings acquired in bidirectional mode to one-directional (1D) for suite2p processing
% The script does this by saving a new block with the Tiffs renamed in one-directional order
% Suite2p expects the Tiff/Frame number to match the plane number:
%
% BOT NUMBERS: 1 2 3 1 2 3
% 1D PLANES:   1 2 3 1 2 3
% 2D PLANES:   1 2 3 3 2 1
%
%%
%path = '\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\VxDG011121F1\2021-04-12\TSeries_VxDG011121F1_FMsweep_bidirectional-003';
cd(path)
[path,folder_name,~] = fileparts(pwd);
mkdir(path, ['1D-'  folder_name]);

filedir = dir;
filenames = {filedir(:).name};
tiffs_str = string(filenames(contains(filenames,'.ome.tif')))';
tiffs = char(tiffs_str);

%The filename format should be as follows:
%BOT_FILENAME-###_Cycle#####_Ch2_######.ome.tif
% ### is the block number
% ##### is the cycle number (used for t-series)
% ###### is the tif/frame number

%Filepath access
planeA = 13; %number of chars to start of frame number
planeB = 8;  %number of chars to end of frame number
cycleA = 23; %number of chars to start of cycle number
cycleB = 19; %number of chars to end of cycle number
chanA  = 15; %number of chars to channel number

filelist_lengths = strlength(deblank(tiffs_str)); %Get length of each filename
L = unique(filelist_lengths); %How many unique lengths are there
if length(L) > 1; error('Check filenames'); end

planes = double(string(tiffs(:,L-planeA:L-planeB)));
cycles = double(string(tiffs(:,L-cycleA:L-cycleB)));
chans = double(string(tiffs(:,L-chanA)));
unique_cycles = unique(cycles);
unique_chans = unique(chans);

%Move each tiff one by one
%Odd cycles are already sequential so filenames don't need to be changed
%The numbering for even cycles has to be reversed
%The last sequence could have < or > nPlanes

moved_files = cell(length(tiffs_str),2); %Record old and new filenames
count = 1;

for c = 1:length(unique_chans)
    
    cc = unique_chans(c);

    tiffs_ch = tiffs_str(chans == cc);
    cycles_ch = cycles(chans == cc);
    planes_ch = planes(chans == cc);
    nFrames = length(tiffs_ch);
    nPlanes = length(tiffs_ch(cycles_ch == 1));

    for f = 1:nFrames
        current_filename = tiffs_ch(f);
        current_cycle = cycles_ch(f);
        current_plane = planes_ch(f);
        cycleFiles = tiffs_ch(cycles_ch == current_cycle,1);
        nPlanesInCycle = length(cycleFiles);

        if current_cycle ~= cycles(end) || nPlanesInCycle <= nPlanes
            if mod(current_cycle,2) == 1 %Odd
                new_filename = current_filename;
            else %Even
                inverted_files = flipud(cycleFiles);
                new_filename = inverted_files(current_plane);
            end
        else 
            if mod(current_cycle,2) == 1 %Odd
                regular_files = cycleFiles(1:nPlanes,:);
                inverted_files = flipud(cycleFiles(nPlanes+1:end,:));
                combined_files = [regular_files; inverted_files];
            else %Even
                inverted_files = flipud(cycleFiles(1:nPlanes,:));
                regular_files = cycleFiles(nPlanes+1:end,:);
                combined_files = [inverted_files; regular_files];
            end
            new_filename = combined_files(current_plane);
        end

        moved_files{count,1} = current_filename;
        moved_files{count,2} = new_filename;

        %copy and rename file
        copyfile(current_filename, [path '/' '1D-'  folder_name '/' char(new_filename)])
        
        count = count + 1;
    end
end

%% Save excel file with old and new filenames to new block
cd([path '/' '1D-' folder_name])
moved_files_save = [{'Old filename', 'New filename'}; string(moved_files)];
xlswrite('BOT_file_identities', moved_files_save);

%% Copy BOT and voltage recording CSVs to new block
csv_str = string(filenames(contains(filenames,'.csv')))';
for c = 1:length(csv_str)
    copyfile([path '/' folder_name '/' char(csv_str(c))], [path '/' '1D-'  folder_name '/' char(csv_str(c))])
end

%% Copy XML files to new block
xml_str = string(filenames(contains(filenames,'.xml')))';
for x = 1:length(xml_str)
    copyfile([path '/' folder_name '/' char(xml_str(x))], [path '/' '1D-'  folder_name '/' char(xml_str(x))])
end

%%

disp('Done processing')
