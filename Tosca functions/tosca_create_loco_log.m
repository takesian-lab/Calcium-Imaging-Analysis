function TL = tosca_create_loco_log(tl)
% TOSCA_CREATE_LOCO_LOG -- parse .loco.txt file into trial-based structure
% Usage: TL = tosca_create_loco_log(tl)
%
% --- Inputs ---
% tl : Tosca log output from tosca_create_log
%
locoFile = strrep(tl.filename, '.txt', '.loco.txt');
if ~exist(locoFile, 'file')
   TL = tl;
   return;
end

loco = tosca_read_loco(locoFile);

itr = find(loco.ch == 1);
if length(itr) ~= length(tl.trials)
   fprintf('Locomotion data has %d/%d trials. Repairing data...\n', length(itr), length(tl.trials));
   loco = tosca_repair_loco_file(locoFile, tl);
end

tloco = loco.t(loco.ch==0);
speed = loco.speed(loco.ch==0);

for k = 1:length(tl.trials)
   tl.trials{k}.loco.t = tloco(tloco>=tl.trials{k}.start & tloco<=tl.trials{k}.stop);
   tl.trials{k}.loco.speed = speed(tloco>=tl.trials{k}.start & tloco<=tl.trials{k}.stop);
   
   for ks = 1:length(tl.trials{k}.states)
      tl.trials{k}.states(ks).loco.t = tloco(tloco>=tl.trials{k}.states(ks).start & tloco<=tl.trials{k}.states(ks).stop);
      tl.trials{k}.states(ks).loco.speed = speed(tloco>=tl.trials{k}.states(ks).start & tloco<=tl.trials{k}.states(ks).stop);
   end
end

TL = tl;

