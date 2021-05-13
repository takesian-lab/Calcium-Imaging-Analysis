function tosca_zone_pref(id, session, run)
% TOSCA_ZONE_PREF -- computes "zone preference" across trials in a run.
% Usage: tosca_zone_pref(id, session, run)
%

% Read parameters and summary data
[D, P] = tosca_read_run(id, session, run);

Ntr = length(D);
Zone1 = zeros(Ntr, 1);
Zone2 = zeros(Ntr, 1);
Ntotal = zeros(Ntr, 1);

for ktr = 1:Ntr,
   T = tosca_read_trial(P, D, ktr);
   
   Ntotal(ktr) = length(T.Zone_1);
   Zone1(ktr) = sum(T.Zone_1);
   Zone2(ktr) = sum(T.Zone_2);
end

figure;
plot(1:Ntr, 100*(Zone1./Ntotal), 'r', 1:Ntr, 100*(Zone2./Ntotal), 'b');
set(gca, 'XLim', [0 Ntr+1], 'YLim', [-5 105]);
xlabel('Trial');
ylabel('Percent time in zone');

title(sprintf('%s-Session%d-Run%d', id, session, run));
leg{1} = sprintf('Zone1: %d%%', round(100*sum(Zone1)/sum(Ntotal)));
leg{2} = sprintf('Zone2: %d%%\n\n', round(100*sum(Zone2)/sum(Ntotal)));
legend(leg);

fprintf('\n%s-Session%d-Run%d\n', id, session, run);
fprintf('Zone1: %d%%\n', round(100*sum(Zone1)/sum(Ntotal)));
fprintf('Zone2: %d%%\n\n', round(100*sum(Zone2)/sum(Ntotal)));
