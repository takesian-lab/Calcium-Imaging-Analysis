function varName = create_valid_varname(str)

varName = deblank(str);
varName = strrep(varName, ' ', '_');
varName = strrep(varName, '---.', '');
varName = strrep(varName, '-', '_minus');
varName = strrep(varName, '(', '');
varName = strrep(varName, ')', '');
varName = strrep(varName, '/', '_');
varName = strrep(varName, '\', '_');
varName = strrep(varName, '%', '');
varName = strrep(varName, '?', '');
varName = strrep(varName, '+', '_plus');
varName = strrep(varName, '#', 'num');
varName = strrep(varName, ':', '.');
varName = strrep(varName, '<', '');
varName = strrep(varName, '>', '');
varName = strrep(varName, char(181), 'u');
if varName(1)>='0' && varName(1)<='9',
    varName = ['Key_' varName];
end

varName = strrep(varName, '._', '.');

varName = strrep(varName, '0x20', '');
varName = strrep(varName, '__', '_');

if varName(1) == '_',
    varName = varName(2:end);
end

