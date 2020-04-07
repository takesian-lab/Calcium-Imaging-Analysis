clear all; close all;
Tosca_file_name = (['VK11_21-Session1-Run2.txt']); %find all the .txt file
 [Data, Params] = tosca_read_run(Tosca_file_name); %run Ken's read_run program to pull out parameters
 
 for i = 1:length(Data)
    nums1(i) = Data{i}.cue.CurrentSource.Level.dB_re_1_Vrms;
end
currents = unique(nums1);

allTrials = {};
for i = 1:length(currents)
    allTrials{i} = NaN(length(Data),3);
end

for i = 1:length(Data)
    
    %Find the current intensity level and freqency
    ctype= Data{i}.cue.CurrentSource.Level.dB_re_1_Vrms;
    ctype=find(currents == ctype);
    allTrials{ctype}(i,1) = Data{i}.cue.CurrentSource.Level.dB_re_1_Vrms;
    
    %Trial type
    if strcmp(Data{i}.Type,'adapt')
        allTrials{ctype}(i,2) = 1; %Regular trial
    else
        allTrials{ctype}(i,2) = 0; %Catch trial
    end
  %Define the outcome of the trial
    outcome = Data{i}.Result;
    if strcmp(outcome,'Go')
        allTrials{ctype}(i,3) = 1;
    else
        allTrials{ctype}(i,3) = 0;
    end
    
    
end

%Now eliminate all NaNs
for i = 1:length(allTrials)
    allTrials{i}(any(isnan(allTrials{i}), 2), :) = [];
end


%We will plot the results for each set of frequencies

for j = 1:length(currents)
    
    %Get the trials from the current frequency type
    nums = allTrials{j};
    for i = 1:length(nums)
        nums(i,4) = i;
    end
    
    %Decide the y axis limits
    yMin = min(nums(nums(:,1)~=-Inf,1)) - 5;
    yMax = max(nums(nums(:,1)~=-Inf,1)) + 5;

    figure

    subplot(1,2,1)
    %Plot regular trials
    tmp = nums(nums(:,2) == 1,:);
    tmp1 = tmp(tmp(:,3) == 1,:);
    p1 = plot(tmp(:,4),tmp(:,1),'-ko','MarkerSize',6,'LineWidth',1.5);
    hold on
    p2 = plot(tmp1(:,4),tmp1(:,1),'ko','MarkerFaceColor','k','MarkerSize',6);

    %Plot catch trials not -1000
    tmp = nums(nums(:,2) == 0,:);
    tmp = tmp(tmp(:,1) ~= -Inf,:);
    tmp1 = tmp(tmp(:,3) == 1,:);
    p3 = plot(tmp(:,4),tmp(:,1),'bo','MarkerSize',6);
    p4 = plot(tmp1(:,4),tmp1(:,1),'bo','MarkerFaceColor','b','MarkerSize',6);

    %Plot catch trials that are -1000
    tmp = nums(nums(:,2) == 0,:);
    tmp = tmp(tmp(:,1) == -Inf,:);
    tmp1 = tmp(tmp(:,3) == 1,:);
    p5 = plot(tmp(:,4),zeros(1,length(tmp(:,4)))+yMin,'ro','MarkerSize',6);
    p6 = plot(tmp1(:,4),zeros(1,length(tmp1(:,4)))+yMin,'ro','MarkerFaceColor','r','MarkerSize',6);
%     ylim([yMin yMax])
    xlabel('Trial Number','FontSize',16)
    ylabel('Current Level (dBV))', 'FontSize', 16)
    title(strcat(num2str(currents(j)),' kHz Track'), 'FontSize',16)

    legend([p1 p2 p3 p4 p5 p6],{'CS+ miss' 'CS+ hit' 'Range miss' 'Range hit' 'CS- withold' 'CS- FA'},...
             'FontSize',12)
    hold off


%     subplot(1,2,2)
%     tmp = nums(nums(:,1) ~= -1000);
%     [counts, binds] = hist(tmp,20);
%     barh(binds,counts)
% %     ylim([yMin yMax])
%     xlabel('# Trials','FontSize',16)
%     ylabel('Sound Level (dB SPL)', 'FontSize', 16)
%     
end

