function [block] = FreqDisc_Behavior_singleblock(block)
% DOCUMENTATION IN PROGRESS
% 
% What does this function do?
% 
% Argument(s): 
%   block (struct)
% 
% Returns:
%   block (struct)
% 
% Notes:
%
% TODO: Remove magic numbers 
% Search 'TODO'

%% Only run this function if this is a behavior block (i.e. stim_protocol == 7)

% find prep trials
if block.setup.stim_protocol==7
    findPrep = block.trialType;
    r = ismember(0, findPrep);
    if r ==0
        block.prepTrial = 1;
    else block.prepTrial = 0;
    end
end

%% find hitrate; determine if it is above threshold
if block.setup.stim_protocol==7
    Hits = find(block.Outcome==1);
    misses =find(block.Outcome==0);
    FA = find(block.Outcome==4);
    withold = find(block.Outcome==3);
    block.HitRate = length(Hits)./(length(Hits)+length(misses));
    block.FARate = length(FA)./(length(FA)+length(withold));
end
%% visualize the data?
%set sound time to zero
% if block.setup.stim_protocol==7
%     sound_time = block.New_sound_times;
%     Trial_time = block.Tosca_times;
%     for i = 1:length(Trial_time)
%         center2sound{i} = Trial_time{1,i}(:) - sound_time(i);
%         centerLoc (i) =  find(center2sound{i}==0);
%         lickWindow(i,:) = (centerLoc(i)-35):(centerLoc(i)+1000);
%         centeredLicks(i,:) = block.lick_time{i,1}(1,lickWindow(1,:));
%     end
%     
%     %cs+ vs cs-
%     CSplus = find(block.trialType==1);
%     CSminus = find(block.trialType==0);
%     for i = 1:length(CSplus)
%         CSplusLicks(i,:) = centeredLicks(CSplus(i),:);
%     end
%     for i = 1:length(CSminus)
%         CSminusLicks(i,:) = centeredLicks(CSminus(i),:);
%     end
    
%     timescale = center2sound{1,1}(lickWindow(1,:),1);
%     
%     subplot(2,1,1)
%     imagesc(CSplusLicks)
%     title('Licks aligned to Sound, CS+')
% 
%     subplot(2,1,2)
%     imagesc(CSminusLicks)
%     title('Licks aligned to Sound, CS-')
    
% end
end