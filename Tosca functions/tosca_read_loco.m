function L = tosca_read_loco(fn)
% TOSCA_READ_LOCO -- reads Tosca locomotion data.
% Usage: L = tosca_read_loco(fn)
%
% *** Inputs ***
% fn : data file path
%
% *** Outputs ***
% L       
% L.t     : time, s
% L.ch    : channel (0=speed, 1=trial marker)
% L.speed : speed, cm/s
%

if ~contains(fn, '.loco')
   fn = strrep(fn, '.txt', '.loco.txt');
end

if ~exist(fn, 'file')
   error('Tosca locomotion file not found.');
end

fp = fopen(fn, 'rt');
C = textscan(fp, '%f%f%f');
fclose(fp);

L.t = C{1};
L.ch = C{2};
L.speed = C{3};


