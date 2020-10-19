% load extracted data
% stimTypes = {'FM','RF','SAM','SAMfreq','NoiseITI','water','air'};
stimTypes = {'FM','RF','SAM','SAMfreq','NoiseITI'};
data_path = ('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS\Inactive');
save_path = ('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS\Inactive');
matchFile = ('Matched Cells v2');
cd(data_path)

for i=1:length(stimTypes)
    s = stimTypes{1,i};
    name = strcat('extractedData_', s, '.mat');
    load(name)
    Data.([s])=ExtractedData;
end



%load matching cells
matchList = importfile(matchFile);
Match.info = matchList;
matchList(1,:) = []; %Remove header


%% pull out unique mice, FOV, and cell numbers
mouseID = [matchList{:,1}]';
uniquemice = unique(mouseID);
count=1;
for i = 1: length(uniquemice)
    mr(i,:) = strcmp(mouseID,uniquemice(i))';
    rows =  find(mr(i,:)==1);
    mouserow{i} = rows;
    fov = [matchList{rows,2}];
    FOVlist{i} = unique(fov);
end

for i = 1:length(uniquemice)
    for k = 1:length(FOVlist{1,i}(:))
        count=1;
        numf = num2str(FOVlist{1,i}(k));
        for j = 1:length(mouserow{1,i}(:))
            rr = mouserow{1,i}(1,j);
            if matchList{rr,2} == FOVlist{1,i}(k)
                Rowlist.([uniquemice(i)]).(['FOV' numf])(count) = rr;
                count = count+1;
            end
        end
    end
end

%% find the cells that were included in the receptive ExtractData set

%look at the rows that correspond to a single mouse
tic
for i =1:length(uniquemice)%loop through mice\
   toc
    for k = 1:length(FOVlist{1,i}(:))%loop through FOV
        numf = num2str(FOVlist{1,i}(k));
        rw =  Rowlist.([uniquemice(i)]).(['FOV' numf])(:);
        for q = 1:length(rw) %only look at rows that = mouse and FOV of interest
            for s = 1:length(stimTypes) %one stimtype at a time
                ss = s+2; %the first two columns are mouse id and FOV
                
                %pull out data that matches the stim type
                nomdat = Data.([stimTypes{s}]).NominalData;
                autodat = Data.([stimTypes{s}]).AutoActivity;
                numdat = Data.([stimTypes{s}]).NumericalData;
                caRast = Data.([stimTypes{s}]).Calcium_Raster;
                sRast = Data.([stimTypes{s}]).Spikes_Raster;
                
                % find the cell in this row...
                cellnum = matchList{rw(q),ss};
          
                % find see if this cell is in the Extracted data
                for n = 1:length(nomdat)
                    if nomdat{n,2} == uniquemice(i); 
                        if nomdat{n,3} == FOVlist{1,i}(k);
%                             if autodat{n,1}~('none');
                                rcell = nomdat{n,6};
%                             end
                        end
                    end
%                     A = exist('rcell');
%                     if A ==1

                        try
                            if cellnum ==rcell
                                Match.([stimTypes{s}]).numerical{rw(q),1} = numdat(n,:);
                                Match.([stimTypes{s}]).Calcium_Raster{rw(q),1} = caRast(n,:);
                                Match.([stimTypes{s}]).Spikes_Raster{rw(q),1} = sRast(n,:);
                            end
                        catch
                            if cellnum == 'NaN'
                                Match.([stimTypes{s}]).numerical{rw(q),1} = 'NaN';
                                Match.([stimTypes{s}]).Calcium_Raster{rw(q),1} = 'NaN';
                                Match.([stimTypes{s}]).Spikes_Raster{rw(q),1} = 'NaN';
                            end
                        end
%                     end
                    clear rcell
                end
            end
            
        end
    end
end
display('done sorting...')




%% plot data - in progress

% make single matrix of avg numerical data
nacheck = 'NaN';
for i = 1:length(stimTypes)
    for j=1:length(Match.([stimTypes{i}]).numerical)
        x = isempty(Match.([stimTypes{i}] ).numerical{j,1});
        if x==1
            meanval(j,i) = -1000;
        elseif isequal(Match.([stimTypes{i}]).numerical{j,1},nacheck);
            meanval(j,i) = NaN;
        else
            meanval(j,i)=nanmean(Match.([stimTypes{i}]).numerical{j,1}(:));
        end
    end
end

% CGS needs to fix this,but it will work for the purposes of APAN
z = find(meanval==-1000);
meanval(z)=0;

%graph - 
data = meanval;

[nr,nc] = size(data);
% ax(1) = subplot(2,1,1); 
% imagesc(data); 
% ax(2) = subplot(2,1,2); 
pcolor([data nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(ax, 'clim', [0 1]);

save_path = ('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS\Inactive');
cd(save_path)
save(['MatchedData_all.mat'], 'Match');