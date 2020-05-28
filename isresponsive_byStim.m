function [data] = isresponsive_byStim(data,std_level)
%
% DOCUMENTATION IN PROGRESS
%
% What does this function do?
% This function determines if a cell is responsive to a particular stimulus
% Argument(s): 
%   data(struct)
%   std_level
%   
% Returns:
%   data(struct)
% 
% Notes:
%
%
% TODO: Magic numbers
% Search 'TODO'
display('...determine responsivity to each stimulus....')
data.setup.std_level_byStim = std_level;
setup = data.setup;

for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)};
    stimdata = data.([mouseID]).stim_df_f;
    
    Var1List=unique(data.([mouseID]).parameters.variable1);%what is Variable#1
    Var2List=unique(data.([mouseID]).parameters.variable2);%what is Variable#2
    data.([mouseID]).parameters.Var1List=Var1List;%store in a list for later
    data.([mouseID]).parameters.Var2List=Var2List;%store in a list for later
    n1=data.([mouseID]).parameters.variable1;%pull out variable#1 in order presented in experiment
    n2=data.([mouseID]).parameters.variable2;%pull out variable#2 in order presented in experiment
    for m=1:length(Var1List)%loop through variable1 (frequency for TRF)
        p=find(n1==Var1List(m));%pull out a particular stimulus (Var #1) (i.e. 4kHz)
        for q=1:length(Var2List)
            r=find(n2==Var2List(q));%pull out a particular stimulus (Var #2) (i.e. 60dB)
            [s]=intersect(p,r); %find specific stim types (i.e. 4khz, 60dB)
            data.([mouseID]).parameters.stimIDX(m,q)={s};%stim index (Var#1xVar#2, i.e. freq x level)
            
           for i=1:size(stimdata.F7_df_f,1);
        all_average(i,:) = squeeze(mean(stimdata.F7_df_f(i,s,:), 2)); %mean of all trials for each cell
        SEM_trace(i,:) = std(stimdata.F7_df_f(i,s,:))./sqrt(size(stimdata.F7_df_f(i,s,:),2));
        
        window_avg(i,:) = mean(stimdata.window_trace(i,s), 2); %average means responses (sound to 2s post sound) across trials for each cell
        peak_avg(i,:) = mean(stimdata.peak_val(i,s), 2); %average peak response
        maxpeak_avg(i,:) = mean(stimdata.max_peak(i,s), 2);%average around max peak
        minpeak_avg(i,:) = mean(stimdata.min_peak(i,s), 2); %average around negative peak
        
        %determine whether each cell is responsive (defined as average response
        %more than defined STDS above baseline)
        base_mean(i,:) = mean(stimdata.baseline(i,s,:),3);
        base_std(i) = std(base_mean(i));
%         f = mean(data.([mouseID]).std_baseline(i,:), 2); %average baseline STD across trials for each cell
         data.([mouseID]).response.isRespPosStim(m,q,i) = maxpeak_avg(i,:) > std_level*mean(base_std(i)) & window_avg(i,:)>0; %will be 0 or 1
         data.([mouseID]).response.isRespNegStim(m,q,i) = minpeak_avg(i,:) < -std_level*mean(base_std(i)) & window_avg(i,:)<0;
    end
            

          end
 
        end
        
    end
end
