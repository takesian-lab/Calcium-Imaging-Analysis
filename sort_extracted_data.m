% load extracted data
stimTypes = {'FM','RF','SAM','SAMfreq','water','air'};
data_path = ('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020');
matchFile = ('Matching cells');
cd(data_path)

for i=1:length(stimTypes)
    s = stimTypes{1,i}
    name = strcat('extractedData_', s, '.mat');
    load(name)
    Data.([s])=ExtractedData;
end



%load matching cells
matchList = importfile(matchFile);
matchList(1,:) = []; %Remove header
Match.info = matchList;

%% find the cells that are responsive and are matched
mouseID = [matchList{:,1}];
uniquemice = unique(mouseID);

count=1
for i = 1:length(uniquemice);
    %find FOVs
    for j = 1:length(matchList);
        if matchList{j,1} == uniquemice(i);
            FOV(count) = matchList{j,2};
            count = count+1;
        end
        FOVlist = unique(FOV);
        % run through each mouse/FOV to build the cell ID matrix
        for k = 1:length(FOVlist);
            for m = 1:length(stimTypes);
                stimcolumn = m+2;
                %find each cell number on the sheet....
                if matchList{j,1} == uniquemice(i);  if matchList{j,2}==FOVlist(k);
                        cellcheck = matchList{j,stimcolumn};
                    end
                end
                
                %check if the cell is in the responsive cell dataset
                nomdat = Data.([stimTypes{m}]).NominalData;
                autodat = Data.([stimTypes{m}]).AutoActivity;
                numdat = Data.([stimTypes{m}]).NumericalData;
                caRast = Data.([stimTypes{m}]).Calcium_Raster;
                sRast = Data.([stimTypes{m}]).Spikes_Raster;
                for n = 1:length(nomdat);
                    if nomdat{n,2} == uniquemice(i);
                        if nomdat{n,3}==FOVlist(k);
                            if autodat{n,1}~('none');
                                rcell = nomdat{n,6};
                                try
                                if cellcheck == rcell;
                                    Match.([stimTypes{m}]).numerical{j,1} = numdat(n,:);
                                    Match.([stimTypes{m}]).Calcium_Raster{j,1} = caRast(n,:);
                                    Match.([stimTypes{m}]).Spikes_Raster{j,1} = sRast(n,:);
                                end
                                catch
                                if cellcheck=='NaN';
                                    Match.([stimTypes{m}]).numerical{j,1} = 'NaN';
                                    Match.([stimTypes{m}]).Calcium_Raster{j,1} = 'NaN';
                                    Match.([stimTypes{m}]).Spikes_Raster{j,1} = 'NaN';
                                end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
