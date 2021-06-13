%% Find blocks of interest/day
info_path = 'Z:\Carolyn\Behavior\SSRI_mice\Info_sheets';
compiled_blocks_path = 'Z:\Carolyn\Behavior\SSRI_mice\compiled blocks\Cn0012621F1';
save_path = 'Z:\Carolyn\Behavior\SSRI_mice\Analyzed data';
info_filename = 'Info_Cn0012621F1';

stim_protocol = 7;
Hit_threshold = 0.7;
remove_dailyPrep = 1;

%right now all of our analysis is by date; however, I am putting this in as
%a way to potentially change the way we look at our data
by_date = 1;

% add in Maryse's method for loading data
cd(info_path)
Info = importfile(info_filename);
[data] = fillSetupFromInfoTable_v3(Info, compiled_blocks_path,stim_protocol);
[options_fit, plotOptions] = setOptions;
options_dp = options_fit;

%the save_path is where the data are
% block_path = 'Z:\Carolyn\Behavior\SSRI_mice\compiled blocks\Cn0012621M3';
% cd(block_path)
% allfiles=dir('*Compiled*');
%% Go through the days
mouselist = data.setup.mousename;
mice = unique(mouselist);
date = string(data.setup.expt_date);
uniquedays = unique(date);
sessionlist = string(data.setup.Session);
sessions = unique(sessionlist);
Outcome = [];
Freq = [];

%%
 %loop through unique mice
for i2 = 1:length(mice)
    mouseIDX = find(strcmp(mouselist, mice{i2}));
    count = 0;
%loop through unique days
    for j2 = 1:length(uniquedays) 
        dayIDX = find(strcmp(date, uniquedays{j2}));
        
        RunsOfInterest = intersect(mouseIDX,dayIDX);
        Outcome = [];
        Freq = [];
        for m2=1:length(RunsOfInterest)
            runnum = num2str(data.setup.Tosca_Runs{(RunsOfInterest(m2))});
            sessionname = strcat('Session', sessionlist{RunsOfInterest(m2)}, '_Run', runnum);
            blockname = strcat('data.', mouselist{RunsOfInterest(m2)}, '.Session', sessionlist{RunsOfInterest(m2)}, '_Run', runnum);
            block = eval(blockname);
            
            % only look at blocks that meet hit threshold and exclude the
            % daily prep trials.
            if block.HitRate>Hit_threshold  & block.prepTrial ==0
                Outcome = [Outcome,block.Outcome];
                Freq = [Freq, block.parameters.variable1];
                target = block.TargetFreq;
                
            end
        end
        
        if ~isempty(Outcome)
            count = count+1;
        end
        
       %% 
        % convert frequencies to difference in octaves
        for i=1:length(Freq)
            if Freq(i)>0
                freq1=Freq(1,i);
                log=log10(freq1);
                logtarg=log10(target);
                diffFreq=abs(log-logtarg);
                percFreq=diffFreq/0.3010;
                percFreq=round(percFreq,1);
                Freq(1,i)= percFreq;
            end
        end
        
        %older scripts had -1000 khz tones - this will remove those if needed
        Freq(Freq(:,1) < 0, :) = [];
        stim_list=unique(Freq(1,:));
        
        
        %%
        %create data structure for psignifit
        datamat(:,1) = stim_list';
        dprimemat(:,1) = stim_list(2:end)';
        


        
             %loop through the stimuli
        for f=1:length(stim_list) 
            stim_percent(f)=nanmean(Outcome(find(Freq(:)==stim_list(f)))); % percent 'go'
            go_n = Outcome(find(Freq(:)==stim_list(f)));
            go_nan = isnan(go_n); %correct for occasional NaN
            go_n(go_nan) = [];%correct for occasional NaN
            stim_n(f,:)=sum(~isnan(find(Freq(1,:)==stim_list(f)))); % n trials per stim
            
            % FA's are identified as a 4 and His are 1 - this makes the math the same for the two conditions
            if f>1 
                stim_percent(f)=stim_percent(f)-3;
                go_n = go_n -3;
            end
           
            
            corrects = find(go_n==1);
            ypsych_size(f) = length(corrects); % number of 'go'
            datamat(f,2) = length(corrects);
            datamat(f,3) = stim_n(f);
            
            % hit rate is rate_response(1) and rest are FA rates
            responseRate(f) = sum(go_n)/stim_n(f);
            
            % dprime measure doesnt allow for hit=1 or =0. This will
            % correct for those such that we dont have inf values
            %Correct floor
            if  responseRate(f) <0.05
                responseRate(f) = 0.05;
            end
            
            %Correct ceiling
            if responseRate(f) >0.95
                responseRate(f) = 0.95;
            end
            
            % correct number of FA/hit trials to match the rate
            response_n(f) = responseRate(f)*stim_n(f);
            datamat(f,2) = response_n(f);
            
            if f ==1 %hits
            z_hit = sqrt(2)*erfinv(2*(responseRate(f))-1);
            else % FAs
                z_FA(f-1) = sqrt(2)*erfinv(2*(responseRate(f))-1);
                dprime(f-1) = z_hit - z_FA(f-1);
            end
        end
        
        
        datamat(:,1) = datamat(:,1)*-1; % make x values negative to correct for fitting 
        dprimemat(:,1) = datamat(2:end,1);
        dprimemat(:,2) = dprime';
        
        
        % generate results for psychometric fit. This will give a threshold
        % at 0.5% go. 
        results_fit = psignifit(datamat,options_fit);
        [threshold_fit,~] = getThreshold(results_fit, results_fit.options.threshPC, 1);
        
         name =strcat(' Day#', num2str(count));
        
        figure;
        subplot(2,1,1)
         suptitle(name)
        h = plotPsych(results_fit,plotOptions);
        
        options_dp.dprimeThresh =1;
        [options_dp,results_dp,zHR] = find_threshPC(datamat,options_dp);
        options_dp = results_dp.options;
          
     

  subplot(2,1,2)
        [x,fitted_yes,fitted_dprime,threshold,slope] = ...
    plotPsych_dprime_cgs(results_dp,dprimemat,results_dp.options,plotOptions,zHR);
       
   
        
      
        
        %store data
        datenumber = strcat('day',num2str(datenum(uniquedays{j2})));
        FrequencyDisc.([mice{i2}]).([datenumber]).Outcome = Outcome;
        FrequencyDisc.([mice{i2}]).([datenumber]).Frequencies = Freq;
        FrequencyDisc.([mice{i2}]).([datenumber]).targetFrequency = target;
        FrequencyDisc.([mice{i2}]).([datenumber]).threshold =threshold_fit;
        
        
    end
end



%% pull out thresholds

for i = 1:length(mice)
    daycount = 0;
    for j = 1:length(uniquedays)
        daynumber = strcat('day',num2str(datenum(uniquedays{j})));
        try
            th = abs(FrequencyDisc.([mice{i}]).([daynumber]).threshold(1));
            daycount = daycount+1;
            FrequencyDisc.allThresh(i,daycount) = th;
        catch
            warning('no thresholds for date')
        end
    end
    
    mousetitle = mice{(i)};
    figure;
    plot (FrequencyDisc.allThresh);
    title([mousetitle]);
end

% 


%%
% Save data
% cd(save_path);
% save([folder '\' filename]);









