
function [data]=plotbystim(setup,data)
for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    x_green=1:size(data.([mouseID]).traces_G,3);
    %plot all cells by stim
    figure;
    for V1=1:length(data.([mouseID]).parameters.Var1List)%loop through stim#1
        for V2=1:length(data.([mouseID]).parameters.Var2List)%loop through stim#2
            stim_list=data.([mouseID]).parameters.stimIDX{V1,V2};
            meanAllbyStim=squeeze(mean(mean(data.([mouseID]).traces_G(:,stim_list,:),2),1));%avg response of positive responsive cells by stim
            std_resp=squeeze(std(mean(data.([mouseID]).traces_G(:,stim_list,:),2),1));
            subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2))
            subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot),
            shadedErrorBar(x_green,smooth(meanAllbyStim,10),smooth(std_resp,10),'lineprops','m')
            stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
            stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
            if setup.stim_protocol==2
                title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
            end
            axis([0 length(x_green) 0 70])
        end
    end


%plot positively responsive cells by stim
figure;

for V1=1:length(data.([mouseID]).parameters.Var1List)%loop through stim#1
    for V2=1:length(data.([mouseID]).parameters.Var2List)%loop through stim#2
        stim_list=data.([mouseID]).parameters.stimIDX{V1,V2};
        stimIDXpos=(data.([mouseID]).parameters.isRespPosStim(V1,V2,:));%index of responsive cells by stim
        
        meanPosResp=squeeze(mean(mean(data.([mouseID]).traces_G(stimIDXpos,stim_list,:),2),1));%avg response of positive responsive cells by stim
        std_resp=squeeze(std(mean(data.([mouseID]).traces_G(stimIDXpos,stim_list,:),2),1));
        subplotSpot=V1+(length(data.([mouseID]).parameters.Var1List)*(length(data.([mouseID]).parameters.Var2List)-V2))
        subplot(length(data.([mouseID]).parameters.Var2List),length(data.([mouseID]).parameters.Var1List),subplotSpot),
        shadedErrorBar(x_green,smooth(meanPosResp,10),smooth(std_resp,10),'lineprops','b')
        stim1=num2str(data.([mouseID]).parameters.Var1List(V1));
        stim2=num2str(data.([mouseID]).parameters.Var2List(V2));
        if setup.stim_protocol==2
            title({sprintf('%s dB',stim2);sprintf('%s kHz',stim1)});
        end
         %axis([0 length(x_green) 0 70])
        
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
