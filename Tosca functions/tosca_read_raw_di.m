function [Y, T] = tosca_read_raw_di(fn)

H = tosca_read_header(fn);

ptsPerRead = H.DAQ.DI_Sampling_Rate_Hz / H.DAQ.Poll_Rate_Hz;

% ndi = length(H.DigitalInputs);
ndi = 3;

fp = fopen(fn, 'rb');

Y = [];

while true,
   y = fread(fp, ndi * ptsPerRead, 'uint8');

   if isempty(y), break; end
   
   Y = [Y; reshape(y, ndi, ptsPerRead)'];
end

fclose(fp);

T = (0:size(Y,1)-1) / H.DAQ.DI_Sampling_Rate_Hz;
