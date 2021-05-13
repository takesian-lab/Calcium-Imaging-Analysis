function H = tosca_read_header(fn)

idx = strfind(fn, '-Trial');
if ~isempty(idx),
   fn = [fn(1:idx-1) '.txt'];
end

% Read parameter sections
H = parse_ini_config(fn);
