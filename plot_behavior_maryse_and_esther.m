%plot_behavior_maryse_and_esther

%% Psychometric curve

%Block data
mousename = block.setup.mousename;
expt_date = block.setup.expt_date;

%Stimulus data
allFrequencies = block.TargetFreq; %get frequencies from block
uniqueFrequencies = unique(allFrequencies); %use the function unique to extract each frequency only once and in numerical order
nRepsPerFrequency = length(allFrequencies)/length(uniqueFrequencies); %how many times each frequency was presented

%Behavioral data
outcomes = block.Outcome; %get outcomes (hit, miss, FP, withholds) from block
misses = outcomes == 0;
hits = outcomes == 1;
withholds = outcomes == 3;
FPs = outcomes == 4;

%calculate hit rate per frequency
hitsPerFrequency = nan(1,length(uniqueFrequencies)); %make empty vector to fill with data
for i = 1:length(uniqueFrequencies)
    freqResponse = hits(allFrequencies == uniqueFrequencies(i));
    hitsPerFrequency(i) = sum(freqResponse);
end
hitRatePerFrequency = hitsPerFrequency/nRepsPerFrequency;

%set up axes for psychometric curve
x = 1:length(hitRatePerFrequency);
y = hitRatePerFrequency;

%make figure
figure
plot(x,y)
ylabel('Hit Rate')
xlabel('Alternating frequency')
set(gca, 'XTick', x)
set(gca, 'XTickLabel', uniqueFrequencies)
title([mousename ' ' expt_date])



