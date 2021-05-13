function S = parse_ini_config(FN, stopAt)
% PARSE_INI_CONFIG -- reconstruct data from INI format.
% Usage: S = parse_ini_config(FN)
%

if nargin < 2, stopAt = '[]'; end

fp = fopen(FN, 'rt');
if fp == -1,
   error('Cannot open file for reading: %s', FN);
end

S = [];
line = fgetl(fp);
while 1,
   if ~isempty(line) && line(1) == '[',
      Section = strrep(line(2:end-1), ' ', '_');
      Section = create_valid_varname(Section);
      [S.(Section), line] = parse_section(fp);
   else
      line = fgetl(fp);
   end
   
   if isequal(line, '[RAW_DATA]') || isequal(line, '[DATA]') || isequal(line,stopAt) || feof(fp),
      break;
   end
end

fclose(fp);

%--------------------------------------------------------------------------
function [P, line] = parse_section(fp)
% PARSE_SECTION
%
P = [];
IsStructArray = false;
while 1,
   if ~IsStructArray,
      line = fgetl(fp);
   end
   IsStructArray = false;
   if (~isempty(line) && isequal(line(1), '[')) || feof(fp),
      break;
   end

   idx = find(line == '=');
   if isempty(idx), continue; end
   
   curKey = create_valid_varname(line(1:idx-1));
   
   [val, IsStructArray] = parse_value(line(idx+1:end));
   if IsStructArray,
      [val, line] = parse_struct_array(fp, line(1:idx-1), val);  
   end
   
   if ischar(val),
      eval(['P.' curKey '=''' val ''';']);
   else
      eval(['P.' curKey '=val;']);
   end
end

names = fields(P);
if length(names) == 1,
   P = P.(names{1});
end

%--------------------------------------------------------------------------
function [Val, IsStructArray] = parse_value(valStr)

IsStructArray = false;

% Parse value
valStr = strtrim(valStr);
if isempty(valStr),
   Val = '';
   
elseif isequal(lower(valStr), 'true'),
   Val = true;
elseif isequal(lower(valStr), 'false'),
   Val = false;
elseif isequal(valStr(1), '<'),
   idx = find(valStr == '>');
   dim = sscanf(valStr(2:idx-1), '%d');
   if isempty(dim),
      Val = valStr;
   else
      dim = dim(:)';
      if length(dim) == 1,
         dim = [1 dim];
      end
      valStr = valStr(idx+1:end);
      Val = sscanf(valStr, '%g,');
      if ~isempty(Val),
         Val = reshape(Val, fliplr(dim));
      elseif prod(dim) > 0,
         Val = dim;
         IsStructArray = true;
      end
   end
   
elseif valStr(1)=='-' || (valStr(1)>='0' && valStr(1)<='9'),
   Val = str2double(valStr);
   if isnan(Val),
      Val = valStr;
   end
   
else
   Val = valStr;
   Val = strrep(Val, 'line*~|.feed', '\n');
end

%--------------------------------------------------------------------------
function [S, line] = parse_struct_array(fp, name, dim)

line = '';

while 1,
   if isempty(line),
      line = fgetl(fp);
   end
   line = strrep(line, '..', '.');
   idx = strfind(line, name);
   ieq = find(line == '=');
   if isempty(idx) || idx(1) > ieq(1), break; end
   line = line(idx:end);
   
   idot = find(line == '.');
   idot = idot(idot > length(name));
   
   ieq = find(line == '=');
   if isempty(idot) || idot(1)>ieq,
      index = 1 + sscanf(line(length(name)+1:ieq(1)-1), '%d');
      val = parse_value(line(ieq(1)+1:end));
      if isscalar(index),
         S{index} = val;
      elseif length(index) == 2,
         S(index(1), index(2)) = val;
      else
         error('Error parsing array.');
      end
      line = '';
      continue;
   end
   
   index = 1 + sscanf(line(length(name)+1:idot(1)-1), '%d');

   curKey = create_valid_varname(line(idot(1)+1:ieq(1)-1));
   
   parentVars = strtrim(strsplit(name, '.'));
   curVars = strsplit(curKey, '.');
   if strcmp(parentVars{end}, curVars{1}) || strcmp(parentVars{end}, [curVars{1} 's']) ...
         || (strcmp(parentVars{end}, 'StimChans') && strcmp(curVars{1}, 'Stimulus')) ...
         || (strcmp(parentVars{end}, 'Families') && strcmp(curVars{1}, 'Family')) ...
         || (strcmp(parentVars{end}, 'Vars') && strcmp(curVars{1}, 'Element')) ...
         || (strcmp(parentVars{end}, 'Flowchart') && strcmp(curVars{1}, 'Flow_Element')) 
      curKey = '';
      for k = 2:length(curVars),
         curKey = [curKey '.' curVars{k}];
      end
   else
      curKey = ['.' curKey];
   end

   
   [val, IsStructArray] = parse_value(line(ieq(1)+1:end));
   if IsStructArray,
      [val, line] = parse_struct_array(fp, line(idot(end)+1:ieq(1)-1));
   else
      line = '';
   end
   
   if ischar(val),
      eval(['S(index)' curKey '=''' val ''';']);
   else
      eval(['S(index)' curKey '=val;']);
   end
   
end

%--------------------------------------------------------------------------
% END OF PARSE_INI_CONFIG.M
%--------------------------------------------------------------------------
