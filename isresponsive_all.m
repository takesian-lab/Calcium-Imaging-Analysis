function [data]=isresponsive_all(data,std_level,run_redcell)
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

data.setup.std_level = std_level;
setup = data.setup;
display('...testing if cells are responsive to all (averaged) stimuli...')


for a=1;%:length(setup.mousename)
    mouseID=setup.mousename{(a)};
    stimdata = data.([mouseID]).stim_df_f;
    
    %go through all of the cells
    for i=1:size(stimdata.F7_df_f,1);
        all_average(i,:) = squeeze(mean(stimdata.F7_df_f(i,:,:), 2)); %mean of all trials for each cell
        SEM_trace(i,:) = std(stimdata.F7_df_f(i,:,:))./sqrt(size(stimdata.F7_df_f(i,:,:),2));
        
        window_avg(i,:) = mean(stimdata.window_trace(i,1:20), 2); %average means responses (sound to 2s post sound) across trials for each cell
        peak_avg(i,:) = mean(stimdata.peak_val(i,:), 2); %average peak response
        maxpeak_avg(i,:) = mean(stimdata.max_peak(i,:), 2);%average around max peak
        minpeak_avg(i,:) = mean(stimdata.min_peak(i,:), 2); %average around negative peak
        
        %determine whether each cell is responsive (defined as average response
        %more than defined STDS above baseline)
        base_mean(i,:) = mean(stimdata.baseline(i,:,:),3);
        base_std(i) = std(base_mean(i),2);
       
        data.([mouseID]).response.isRespPos(i) = maxpeak_avg(i,:) > std_level*mean(base_std(i)) & window_avg(i,:)>0; %will be 0 or 1
        data.([mouseID]).response.isRespNeg(i) = minpeak_avg(i,:) < -std_level*mean(base_std(i)) & window_avg(i,:)<0;
       
       
     
    end
    
    
    
        first_cell = 1;
        last_cell = size(all_average,1);
        num_cells = last_cell-first_cell+1;
        
        figure;
        
        for i=first_cell:last_cell
            %plot mean traces across cells with response means - responsive cells are
            %green, non-responsive are black
            subplot_num = i-first_cell+1;
            subplot(ceil(sqrt(num_cells)),ceil(sqrt(num_cells)),subplot_num);
            x_green = 1:length(all_average(i,:));
            if data.([mouseID]).response.isRespPos(i) == 1
                shadedErrorBar(x_green,smooth((all_average(i,:)),10),smooth((SEM_trace(i,:)),10),'lineprops','-b','transparent',1); hold on;
            else
                if data.([mouseID]).response.isRespNeg(i)== 1
                    shadedErrorBar(x_green,smooth((all_average(i,:)),10),smooth((SEM_trace(i,:)),10),'lineprops','-c','transparent',1); hold on;
                else
                    shadedErrorBar(x_green,smooth((all_average(i,:)),10),smooth((SEM_trace(i,:)),10),'lineprops','-k','transparent',1); hold on;
                end
            end
%             
%             plot(peak_avg(i),'o'); %plot average response
%             plot(maxpeak_avg(i),'o'); %plot average around peak
%             plot(minpeak_avg(i),'o');
        end
        
     
        all_average = squeeze(mean(mean(stimdata.F7_df_f(:,:,:),2),1));%mean across all cells
        isresponsive_avg= squeeze(mean(mean(stimdata.F7_df_f((data.([mouseID]).response.isRespPos),:,:),2),1));%mean across isRespPos cells
        isresponsiveneg_avg= squeeze(mean(mean(stimdata.F7_df_f((data.([mouseID]).response.isRespNeg),:,:),2),1));%mean across isRespNeg cells
        
        a_green_std = squeeze(std(mean(stimdata.F7_df_f(:,:,:),2),1))./sqrt(size(stimdata.F7_df_f,1));%mean across all cells
        b_green_std= squeeze(std(mean(stimdata.F7_df_f((data.([mouseID]).response.isRespPos),:,:),2),1))./sqrt(size(data.([mouseID]).response.isRespPos,1));%mean across isRespPos cells
        c_green_std= squeeze(std(mean(stimdata.F7_df_f((data.([mouseID]).response.isRespNeg),:,:),2),1))./sqrt(size(data.([mouseID]).response.isRespNeg,1));%mean across isRespNeg cells
        
        try
            figure
            x_green=1:length(all_average);
            subplot(1,3,1); hold on
            title('All cells')
            xlabel('Frames')
            ylabel('Delta F/F')
            shadedErrorBar(x_green,smooth(all_average,10),smooth(a_green_std,10),'lineprops','m');
            subplot(1,3,2); hold on
            title('Positively responding cells')
            xlabel('Frames')
            shadedErrorBar(x_green,smooth(isresponsive_avg,10),smooth(b_green_std,10),'lineprops','b');
            subplot(1,3,3); hold on
            title('Negatively responding cells')
            xlabel('Frames')
            shadedErrorBar(x_green,smooth(isresponsiveneg_avg,10),smooth(c_green_std,10),'lineprops','c');
        catch
            disp(['Skipping mouse ' mouseID ' graph, not enough cells?']);
        end
        
        
    
      clear all_average SEM_trace window_avg peak_avg...
            minpeak_avg maxpeak_avg base_mean  base_std   
    end
end
