function [Data, Params] = tosca_read_run(fn, check)
% TOSCA_READ_RUN -- reads summary data for Tosca run.
% Usage: [Data, Params] = tosca_read_run(ExpID, Session, Run)
% 
% Inputs:
%   fn:      run data file, e.g. 'XYZ999-Session1-Run12.txt'
%
% Outputs:
%   Data:    structure containing summary data
%   Params:  structure containing exhaustive description of Tosca setup
%

if nargin < 2, check = true; end

if ~exist(fn, 'file')
   error(['Cannot find data file: ' fn]);
end

% Read parameter sections
Params = parse_ini_config(fn);
Params.Info.Filename = fn;

% Read data section
fp = fopen(fn, 'rt');
while 1
   line = fgetl(fp);
   if isequal(line, '[DATA]') || feof(fp)
      break;
   end
end

if feof(fp)
   Data = {};
   warning(['Could not find [DATA] section in file: ' fn]);
end

if Params.Info.Version < 2777 %originally 2734, CGS changed to 2776 to get the air to run
   isAdapt = isfield(Params.Schedule, 'Mode') && strcmpi(Params.Schedule.Mode, 'adapt');
else
   isAdapt = strcmpi(Params.Tosca.Schedule_Mode, 'adapt');
end

if Params.Info.Version < 1988
   % Parse header line
   s = fgetl(fp);
   c = textscan(s, '%s', 'delimiter', '\t');
   names = c{1};
   for k = 1:length(names)
      c = textscan(names{k}, '%s');
      names{k} = c{1}{1};
      if isequal(names{k}, '---')
         names{k} = 'Y';
      end
   end
   
   resultColumn = find(strcmp('Result', names), 1);
   if isempty(resultColumn)
      error('Could not find "Result" column in data file: %s', fn);
   end
   
   % Read data
   nd = 1;
   while 1
      s = fgetl(fp);
      if isempty(s) || ~ischar(s), break; end
      
      x = textscan(s, '%s', 'delimiter', '\t');
      for k = 1:resultColumn-1
         Data(nd).(strrep(names{k},'.','')) = str2double(x{1}{k});
      end
      Data(nd).Result = x{1}{resultColumn};
      c = textscan(x{1}{resultColumn+1}, '%s', 'delimiter', ',');
      Data(nd).History = c{1};
      nd = nd + 1;
   end
   d = Data;
   Data = {};
   for k = 1:length(d)
      Data{k} = d(k);
   end
   
else
   Data = {};
   while 1
      d = read_one_trial(fp, isAdapt);
      if isempty(d)
         break;
      else
         Data{end+1} = d;
      end
   end
end
fclose(fp);

% Match trial with .di.txt file
if check && ~isempty(Data) && (~isfield(Params.Info, 'AlignmentChecked') || Params.Info.Version < 1988),
   Data = tosca_check_alignment(Params, Data);
end

if check
   ikeep = [];
   for k = 1:length(Data)
      if Data{k}.N >=0
         ikeep(end+1) = k;
      end
   end
   Data = Data(ikeep);
end

%--------------------------------------------------------------------------
function T = read_one_trial(fp, isAdapt)

T = [];

while 1
   s = fgetl(fp);
   if feof(fp), return; end
   if ~isempty(s) && s(1) == '[', break; end
end

if isAdapt
   n = deal(sscanf(s, '[Track%d Trial%d]'));
else
   n = deal(sscanf(s, '[Block%d Trial%d]'));
end
T.block = n(1);
T.trial = n(2);

while 1
   s = fgetl(fp);
   if isempty(s), break; end
   
   idx = find(s == '=');
   if isempty(idx), break; end
   
   curKey = create_valid_varname(s(1:idx-1));
   
   val = parse_value(s(idx+1:end));
   
   if ischar(val)
      eval(['T.' curKey '=''' val ''';']);
   else
      eval(['T.' curKey '=val;']);
   end
   
   if feof(fp), break; end
end

if isfield(T, 'History')
   c = textscan(T.History, '%s', 'delimiter', ',');
   T.History = c{1};
end

if isfield(T, 'N') && isequal(T.N, 'NaN')
   T.N = NaN;
end

%--------------------------------------------------------------------------
function Val = parse_value(valStr)

% Parse value
valStr = strtrim(valStr);
if isempty(valStr)
   Val = '';
   
elseif isequal(lower(valStr), 'true')
   Val = true;
elseif isequal(lower(valStr), 'false')
   Val = false;
elseif isequal(valStr(1), '<')
   idx = find(valStr == '>');
   dim = sscanf(valStr(2:idx-1), '%d');
   dim = dim(:)';
   if length(dim) == 1
      dim = [1 dim];
   end
   valStr = valStr(idx+1:end);
   Val = sscanf(valStr, '%g,');
   if ~isempty(Val)
      Val = reshape(Val, fliplr(dim));
   elseif prod(dim) > 0
      Val = dim;
%       IsStructArray = true;
   end

elseif valStr(1)=='-' || (valStr(1)>='0' && valStr(1)<='9')
   Val = str2double(valStr);
   if isnan(Val)
      Val = valStr;
   end
   
else
   Val = valStr;
   Val = strrep(Val, 'line*~|.feed', '\n');
end

%--------------------------------------------------------------------------
% END OF tosca_read_run.m
%--------------------------------------------------------------------------
