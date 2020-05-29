function [data] = plot_loco_by_stim(data);
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
        loc_data(:,:,:) = stimdata.F7_df_f(:,Loc_trial,:);
        loc_trace(:,:) = squeeze(mean(loc_data,2)); %Active trials
        loc_sem = std(loc_trace)./sqrt(size(loc_trace,1));
        loc_avg = mean(loc_trace,1);
        
        
        %non-active trials
        noloc_data(:,:,:) = stimdata.F7_df_f(:,noLoc_trial,:);
        noloc_trace(:,:) = squeeze(mean(noloc_data,2)); %Active trials
        noloc_sem = std(noloc_trace)./sqrt(size(noloc_trace,1));
        noloc_avg = mean(noloc_trace,1);
        
        
        x=1:length(loc_avg);
        
        figure
        
        subplot(2,1,1); hold on
        title([mouseID ' FOV ' num2str(FOV)])
        shadedErrorBar(x,smooth(loc_avg,5),smooth(loc_sem,5),'lineprops','-b','transparent',1);
        legend({'Active trials'})
        xlabel('Frames')
        ylabel('Delta F/F')
        
        subplot(2,1,2); hold on
        shadedErrorBar(x,smooth(noloc_avg,5),smooth(noloc_avg,5),'lineprops','-k','transparent',1);
        legend({'Inactive trials'})
        xlabel('Frames')
        ylabel('Delta F/F')
    end
end
end
