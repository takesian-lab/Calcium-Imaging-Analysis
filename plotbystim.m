function [data]=plotbystim(data,run_redcell)
%
% DOCUMENTATION IN PROGRESS
%
% What does this function do?
%
% Argument(s):
%   data(struct)
%
% Returns:
%   data(struct)
%
% Notes:
%
%
% TODO: Magic numbers
% Search 'TODO'

setup = data.setup;

for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)};
    stimdata = data.([mouseID]).stim_df_f;
    x_green=1:size(stimdata.F7_df_f,3);
    if run_redcell == 1;
        redIDX = data.([mouseID]).parameters.red_idx;
        greenIDX = data.([mouseID]).parameters.green_idx;
    end
    %plot all cells by stim
    figure;
    for V1=1:length(data.([mouseID]).parameters.Var1List)%loop through stim#1
        for V2=1:length(data.([mouseID]).parameters.Var2List)%loop through stim#2
            stim_list=data.([mouseID]).parameters.stimIDX{V1,V2};
            % graph all cells (responsive and non-responsive by stim
            if run_redcell == 0;
                meanAllbyStim=squeeze(mean(mean(stimdata.F7_df_f(:,stim_list,:),2),1));%avg response of positive responsive cells by stim
                std_resp=squeeze(std(mean(stimdata.F7_df_f(:,stim_list,:),2),1));
                subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2))
                subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot),
                shadedErrorBar(x_green,smooth(meanAllbyStim,10),smooth(std_resp,10),'lineprops','m')
                stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
                stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
                if setup.stim_protocol==2
                    title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
                end
                %             axis([0 length(x_green) 0 70])
            else
                meanAllbyStim_g=squeeze(mean(mean(stimdata.F7_df_f(greenIDX,stim_list,:),2),1));%avg response of positive responsive cells by stim
                meanAllbyStim_r=squeeze(mean(mean(stimdata.F7_df_f(redIDX,stim_list,:),2),1));%avg response of positive responsive cells by stim
                
                std_resp_g=squeeze(std(mean(stimdata.F7_df_f(greenIDX,stim_list,:),2),1));
                std_resp_r=squeeze(std(mean(stimdata.F7_df_f(redIDX,stim_list,:),2),1));
                
                subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2));
                subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot);
                
                shadedErrorBar(x_green,smooth(meanAllbyStim_g,10),smooth(std_resp_g,10),'lineprops','g'), hold on
                shadedErrorBar(x_green,smooth(meanAllbyStim_r,10),smooth(std_resp_r,10),'lineprops','m')
                
                stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
                stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
                if setup.stim_protocol==2
                    title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
                end
                
            end
        end
    end
    
    
    %plot positively responsive cells by stim
    figure;
    
    
    for V1=1:length(data.([mouseID]).parameters.Var1List)%loop through stim#1
        for V2=1:length(data.([mouseID]).parameters.Var2List)%loop through stim#2
            stim_list=data.([mouseID]).parameters.stimIDX{V1,V2};
            stimIDXpos=(data.([mouseID]).response.isRespPosStim(V1,V2,:));%index of responsive cells by stim
            
            
            if run_redcell == 0
                meanPosResp=squeeze(mean(mean(stimdata.F7_df_f(stimIDXpos,stim_list,:),2),1));%avg response of positive responsive cells by stim
                std_resp=squeeze(std(mean(stimdata.F7_df_f(stimIDXpos,stim_list,:),2),1));
                subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2))
                subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot),
                shadedErrorBar(x_green,smooth(meanPosResp,10),smooth(std_resp,10),'lineprops','b')
                stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
                stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
                if setup.stim_protocol==2
                    title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
                end
                %axis([0 length(x_green) 0 70])
                
            else
                a = find(greenIDX==1);
                b = find(redIDX==1);
                c = find(stimIDXpos==1);

                g_pos=ismember(a,c);
                r_pos= ismember(b,c);
%                 g_pos = find(greenIDX==1 && stimIDXpos==1); %green, positively responsive cells
%                 r_pos = find(redIDX==1 && stimIDXpos==1);   % red, positively responsive cells
                
                meanPosResp_g=squeeze(mean(mean(stimdata.F7_df_f(g_pos,stim_list,:),2),1));%avg response of positive responsive cells by stim
                meanPosResp_r=squeeze(mean(mean(stimdata.F7_df_f(r_pos,stim_list,:),2),1));%avg response of positive responsive cells by stim
                
                std_resp_g=squeeze(std(mean(stimdata.F7_df_f(g_pos,stim_list,:),2),1));
                std_resp_r=squeeze(std(mean(stimdata.F7_df_f(r_pos,stim_list,:),2),1));
                
                
                subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2));
                subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot);
                
                
                shadedErrorBar(x_green,smooth(meanPosResp_g,10),smooth(std_resp_g,10),'lineprops','b'), hold on
                shadedErrorBar(x_green,smooth(meanPosResp_r,10),smooth(std_resp_r,10),'lineprops','k')
                
                stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
                stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
                if setup.stim_protocol==2
                    title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
                end
                
            end
        end
    end
    %plot negatively responsive cells by stim
    
    
    
    %     x_green=size(data.([mouseID]).traces_G,3);
    %     %             subplot_loc=m+(length(m)*(length(q)-q));
    %     subplot(length(Var2List),length(Var1List),m*q);
    %     meanPosResp=squeeze(mean(mean(data.([mouseID]).traces_G(stimIDXpos,:,:),2),1));
    %     plot(x_green,meanPosResp,'m');
    
    
    
    
    
end
end
