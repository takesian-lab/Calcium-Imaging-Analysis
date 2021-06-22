function [isActiveComparedToBlanks, pValue] = compare_to_blank_trials(trials, blank_trials)

%Return NaNs if function cannot be completed for any reason
isActiveComparedToBlanks = nan;
pValue = nan;

%Number of blank trials should exceed number of actual trials to bootstrap
nTrials = size(trials,1);
nBlankTrials = size(blank_trials,1);

if nBlankTrials <=
N = min(nTrials, nBlankTrials);



end