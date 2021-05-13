function S = LoadThorSyncData(filepath, varargin)

varsToRead = {};

idx = find(strcmp(varargin, 'vars'));
if ~isempty(idx),
   varsToRead = varargin{idx+1};
   varargin = varargin(setdiff(1:length(varargin), idx+[0 1]));
end

% Load params from XML:
clockRate = 20000000;
try
   sampleRate = LoadSyncXML(fileparts(filepath));
catch
   sampleRate = 30000;
end

% Start loading HDF5:
info = h5info(filepath);

% Parse input:
props = {'start','length','interval'};
% data = {[1 1],[1 Inf],[1 1]};
data = {[1 1],[1 Inf],[1 1]};

if ~isempty(varargin),
    assert(rem(length(varargin),2)==0 && iscellstr(varargin(1:2:end)), 'Inputs failed to conform to expected string-value pair format, eg> ''start'',1');
    %foundProps = intersect(varargin(1:2:end), props); 
    IdxCell = cellfun(@(x) strcmpi(x,props),varargin(1:2:end),'UniformOutput',false);
    val = double(cell2mat(varargin(2:2:end)))*sampleRate;
    for i=1:length(val)
        data{cell2mat(IdxCell(i))>0} = [1 val(i)];
    end
end

% Read HDF5:

for j = 1:length(info.Groups),
    for k = 1:length(info.Groups(j).Datasets)
        datasetPath = strcat(info.Groups(j).Name,'/',info.Groups(j).Datasets(k).Name);
        datasetName = info.Groups(j).Datasets(k).Name;
        datasetName(isspace(datasetName))='_';   
        if strcmp(info.Groups(j).Name,'/Global'),
            datasetName = 'time';
        end
        
        if ~isempty(varsToRead) && ~any(strcmp(varsToRead, datasetName)),
           continue;
        end
        
        datasetValue = h5read(filepath, datasetPath, data{1}, data{2}, data{3})';
        
        % load digital line in binary:
        if strcmp(info.Groups(j).Name,'/DI'),
            datasetValue(datasetValue>0) = 1;
        end
        
        % create time variable out of gCtr,
        % account for 20MHz sample rate:
        if strcmp(info.Groups(j).Name,'/Global'),
            datasetValue = double(datasetValue)./clockRate;
        end
        
        S.(datasetName) = datasetValue;
    end
end

