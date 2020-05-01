function [parameters] = indexStimuli(parameters,setup)
  FreqList=unique(parameters.frequency_list);%what is Variable#1
    LevList=unique(parameters.level_list);%what is Variable#
   parameters.frequencies=FreqList;%store in a list for later
  parameters.levels=LevList;%store in a list for later
    n1=parameters.frequency_list;%pull out variable#1 in order presented in experiment
    n2=parameters.level_list;%pull out variable#2 in order presented in experiment
    for m=1:length(FreqList)%loop through variable1 (frequency for TRF
        p=find(n1==FreqList(m));%pull out a particular stimulus (Var #1) (i.e. 4kHz)
        for q=1:length(LevList)
            r=find(n2==LevList(q));%pull out a particular stimulus (Var #2) (i.e. 60dB)
            [s]=intersect(p,r); %find specific stim types (i.e. 4khz, 60dB)
            parameters.stimIDX(m,q)={s};%stim index (Var#1xVar#2, i.e. freq x level)
        end
    end

end