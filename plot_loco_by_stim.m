function [data] = plot_loco_by_stim(data);
setup = data.setup;

for a=1:size(setup.mousename,1)
    for b=1:size(setup.mousename,2)
        
        if isempty(setup.mousename{a,b})
            continue;
        end
        
        mouseID=setup.mousename{a,b};
        FOV=setup.FOVs{a,b};
        
        if data.setup.run_redcell==1
            red_idx=data.([mouseID]).parameters.red_idx;
            green_idx=data.([mouseID]).parameters.green_idx;
        end
        
        %%
        stimdata = data.([mouseID]).stim_df_f;
        Loc_trial = data.([mouseID]).cat.isLoco_cat;
        noLoc_trial = ~data.([mouseID]).cat.isLoco_cat;
        %%
%         if data.setup.run_redcell==0
            % active trials
            loc_data(:,:,:) = stimdata.F7_df_f(:,Loc_trial,:);
            loc_trace(:,:) = squeeze(mean(loc_data,2)); %Active trials
            loc_sem = std(loc_trace)./sqrt(size(loc_trace,1));
            loc_avg = mean(loc_trace,1);
            
            
            %non-active trials
            noloc_data(:,:,:) = stimdata.F7_df_f(:,noLoc_trial,:);
            noloc_trace(:,:) = squeeze(mean(noloc_data,2)); %inActive trials
            noloc_sem = std(noloc_trace)./sqrt(size(noloc_trace,1));
            noloc_avg = mean(noloc_trace,1);
            
            
             x=1:length(loc_avg);
        
        figure
        
        subplot(2,1,1); hold on
        title([mouseID ' FOV ' num2str(FOV)])
        shadedErrorBar(x,smooth(loc_avg,10),smooth(loc_sem,10),'lineprops','-b','transparent',1);
        legend({'Active trials'})
        xlabel('Frames')
        ylabel('Delta F/F')
        
        subplot(2,1,2); hold on
        shadedErrorBar(x,smooth(noloc_avg,10),smooth(noloc_sem,10),'lineprops','-k','transparent',1);
        legend({'Inactive trials'})
        xlabel('Frames')
        ylabel('Delta F/F')
        
        data.([mouseID]).response.loco_data = loc_data;
        data.([mouseID]).response.loc_trace = loc_trace;
        data.([mouseID]).response.noloc_data = noloc_data;
        data.([mouseID]).response.noloc_trace = noloc_trace;
        
        clear loc_data loc_trace noloc_data noloc_trace
        
%         end
        
%         %% red cell
%         if data.setup.run_redcell==1
%             % active trials -red
%             r_loc_data(:,:,:) = stimdata.F7_df_f(red_idx,Loc_trial,:);
%             r_loc_trace(:,:) = squeeze(mean(r_loc_data,2)); %Active trials
%             r_loc_sem = std(r_loc_trace)./sqrt(size(r_loc_trace,1));
%             r_loc_avg = mean(r_loc_trace,1);
%             
%             % active trials -green
%             g_loc_data(:,:,:) = stimdata.F7_df_f(green_idx,Loc_trial,:);
%             g_loc_trace(:,:) = squeeze(mean(g_loc_data,2)); %Active trials
%             g_loc_sem = std(g_loc_trace)./sqrt(size(g_loc_trace,1));
%             g_loc_avg = mean(g_loc_trace,1);
%             
%             
%             
%             %non-active trials - red
%             r_noloc_data(:,:,:) = stimdata.F7_df_f(red_idx,noLoc_trial,:);
%             r_noloc_trace(:,:) = squeeze(mean(r_noloc_data,2)); %inActive trials
%             r_noloc_sem = std(r_noloc_trace)./sqrt(size(r_noloc_trace,1));
%             r_noloc_avg = mean(r_noloc_trace,1);
%             
%             %non-active trials - green
%             g_noloc_data(:,:,:) = stimdata.F7_df_f(green_idx,noLoc_trial,:);
%             g_noloc_trace(:,:) = squeeze(mean(g_noloc_data,2)); %inActive trials
%             g_noloc_sem = std(g_noloc_trace)./sqrt(size(g_noloc_trace,1));
%             g_noloc_avg = mean(g_noloc_trace,1);
%             
%             
%             
%             
%              x=1:length(r_loc_avg);
%         
%         figure
%         
%         subplot(2,1,1); hold on
%         title([mouseID ' FOV ' num2str(FOV)])
%         shadedErrorBar(x,smooth(r_loc_avg,10),smooth(r_loc_sem,10),'lineprops','-m','transparent',1); hold on
%         shadedErrorBar(x,smooth(r_noloc_avg,10),smooth(r_noloc_sem,10),'lineprops','-k','transparent',1); 
%         
%         legend({'Active trials','Inactive trials'})
%         xlabel('Frames')
%         ylabel('Delta F/F')
%         
%         subplot(2,1,2); hold on
%         shadedErrorBar(x,smooth(g_loc_avg,10),smooth(g_loc_sem,10),'lineprops','-g','transparent',1); hold on
%         shadedErrorBar(x,smooth(g_noloc_avg,10),smooth(g_noloc_sem,10),'lineprops','-b','transparent',1);
%         legend({'Active trials','Inactive trials'})
%         xlabel('Frames')
%         ylabel('Delta F/F')
%         
%         data.([mouseID]).response.loco_data.red = r_loc_data;
%         data.([mouseID]).response.loco_data.green = g_loc_data;
%         data.([mouseID]).response.loc_trace.red = r_loc_trace;
%         data.([mouseID]).response.loc_trace.green = g_loc_trace;
%         data.([mouseID]).response.noloc_data.red = r_noloc_data;
%          data.([mouseID]).response.noloc_data.green = g_noloc_data;
%         data.([mouseID]).response.noloc_trace.red = r_noloc_trace;
%         data.([mouseID]).response.noloc_trace.green = g_noloc_trace;
%         
%         clear r_loc_data r_loc_trace r_noloc_data r_noloc_trace...
%             g_loc_data g_loc_trace g_noloc_data g_noloc_trace
%         
%         end
%             
        end
        
        %% move to red cell ==0
       
    end
end

