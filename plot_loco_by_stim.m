setup = data.setup;

for a=1:size(setup.mousename,1)
    for b=1:size(setup.mousename,2)
        
        if isempty(setup.mousename{a,b})
            continue;
        end
        
        mouseID=setup.mousename{a,b};
        FOV=setup.FOVs{a,b};
        
        stimdata = data.([mouseID]).stim_df_f;
        Loc_trial = data.([mouseID]).cat.isLoco_cat;
        noLoc_trial = ~data.([mouseID]).cat.isLoco_cat;
        
        % active trials
        for i=1:size(stimdata.F7_df_f,1)
            loc_avg(i,:,:) = mean(stimdata.F7_df_f(:,Loc_trial,:),2); %Active trials
        end
        
        loc_trace = mean(loc_avg,1);
        loc_sem = std(loc_avg)./sqrt(size(loc_avg,1));
        
        
        %non-active trials
        for i=1:size(stimdata.F7_df_f,1)
            noloc_avg(i) = mean(stimdata.F7_df_f(:,noLoc_trial,:),2); %Active trials
        end
        
        noloc_trace = mean(noloc_avg,1);
        noloc_sem = std(noloc_avg)./sqrt(size(noloc_avg,1));
        
        
        x=1:length(loc_trace);
        
        figure
        
        subplot(2,1,1); hold on
        title([mouseID ' FOV ' num2str(FOV)])
        shadedErrorBar(x,smooth(loc_trace,5),smooth(loc_sem,5),'lineprops','-b','transparent',1);
        legend({'Active trials'})
        xlabel('Frames')
        ylabel('Delta F/F')
        
        subplot(2,1,2); hold on
        shadedErrorBar(x,smooth(noloc_trace,5),smooth(noloc_trace,5),'lineprops','-k','transparent',1);
        legend({'Inactive trials'})
        xlabel('Frames')
        ylabel('Delta F/F')
    end
end
