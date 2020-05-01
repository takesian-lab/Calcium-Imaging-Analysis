function [data] = isresponsive_all_wisam(data,setup,std_level)
% [data] = isresponsive_all_wisam(data,setup,std_level)
%
% THIS DOCUMENTATION IS A WORK IN PROGRESS:
%   This function operates across all mice and ROIs for a
%   given stimulus.
%   This function tests whether ROI (green cells) time series are responsive to a stimulus.
%
% Details:
%   This is calculated by a threshold (in SDs) of the mean responses above baseline (across trials).
% 
% Arguments: 
%   data (struct)
%   setup (struct)
%   std_level (float)
% 
% Return:
%   data (struct)
% 
% Search 'TODO'

% For a given stimulus 
% Loop through all the mice
for a=1:length(setup.mousename)
    
    % Grab the current mouse ID
    mouseID = setup.mousename{(a)};
    
    % Pre-allocate arrays
    % TODO: memory allocation
    a_green = [];
    SEM_a_green = [];
    b_green = [];
    c_green = [];
    d_green = [];
    e_green = [];
    
    % For a given mouse 
    % Loop through all of the ROIs (green cells)
    for i=1:size(data.([mouseID]).traces_G,1)
        
        % Mean of all trials for each ROI
        a_green(i,:) = squeeze(mean(data.([mouseID]).traces_G(i,:,:), 2));
        % Standard error of the mean (SEM) of all trials for each ROI
        SEM_a_green(i,:) = std(data.([mouseID]).traces_G(i,:,:))./sqrt(size(data.([mouseID]).traces_G(i,:,:),2));
        % TODO: What does this mean?
        % Average means responses (sound to 2s post sound) across trials for each cell
        b_green(i,:) = mean(data.([mouseID]).response(i,:), 2);
        % Average peak response
        c_green(i,:) = mean(data.([mouseID]).peak_G(i,:), 2); 
        % Average around max peak
        d_green(i,:) = mean(data.([mouseID]).peak_maxG(i,:), 2);
        % Average around negative peak
        e_green(i,:) = mean(data.([mouseID]).peak_minG(i,:), 2); 
        
        % Determine whether each cell is responsive (defined as average response
        % more than defined STDs above baseline)
        
        % Average baseline STD across trials for each cell
        f = mean(data.([mouseID]).std_baseline(i,:), 2); 
        data.([mouseID]).parameters.isRespPos(i) = d_green(i,:) > std_level*mean(f) & b_green(i,:)>0; %will be 0 or 1
        data.([mouseID]).parameters.isRespNeg(i) = e_green(i,:) < -std_level*mean(f) & b_green(i,:)<0;
    end
    first_cell = 1;
    last_cell = size(data.([mouseID]).response,1);
    num_cells = last_cell-first_cell+1;
    
    figure;
    
    for i=first_cell:last_cell
        %plot mean traces across cells with response means - responsive cells are
        %green, non-responsive are black
        subplot_num = i-first_cell+1;
        subplot(ceil(sqrt(num_cells)),ceil(sqrt(num_cells)),subplot_num);
        x_green = 1:length(a_green(i,:));
        if data.([mouseID]).parameters.isRespPos(i) == 1
            shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-b','transparent',1); hold on;
        else
            if data.([mouseID]).parameters.isRespNeg(i)== 1
                shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-c','transparent',1); hold on;
            else
                shadedErrorBar(x_green,smooth((a_green(i,:)),10),smooth((SEM_a_green(i,:)),10),'lineprops','-k','transparent',1); hold on;
            end
        end
        
        plot(b_green(i),'o'); %plot average response
        plot(d_green(i),'o'); %plot average around peak
        plot(e_green(i),'o');
    end
    
  clear a_green b_green c_green d_green e_green
    a_green = squeeze(mean(mean(data.([mouseID]).traces_G(:,:,:),2),1));%mean across all cells
    b_green= squeeze(mean(mean(data.([mouseID]).traces_G((data.([mouseID]).parameters.isRespPos),:,:),2),1));%mean across isRespPos cells
    c_green= squeeze(mean(mean(data.([mouseID]).traces_G((data.([mouseID]).parameters.isRespNeg),:,:),2),1));%mean across isRespNeg cells
    
    a_green_std = squeeze(std(mean(data.([mouseID]).traces_G(:,:,:),2),1));%mean across all cells
    b_green_std= squeeze(std(mean(data.([mouseID]).traces_G((data.([mouseID]).parameters.isRespPos),:,:),2),1));%mean across isRespPos cells
    c_green_std= squeeze(std(mean(data.([mouseID]).traces_G((data.([mouseID]).parameters.isRespNeg),:,:),2),1));%mean across isRespNeg cells
    
    try
    figure
    x_green=1:length(a_green);
    subplot(1,3,1); hold on
        title('All cells')
        xlabel('Frames')
        ylabel('Delta F/F')
        shadedErrorBar(x_green,smooth(a_green,10),smooth(a_green_std,10),'lineprops','m');
    subplot(1,3,2); hold on
        title('Positively responding cells')
        xlabel('Frames')
        shadedErrorBar(x_green,smooth(b_green,10),smooth(b_green_std,10),'lineprops','b');
    subplot(1,3,3); hold on
        title('Negatively responding cells')
        xlabel('Frames')
        shadedErrorBar(x_green,smooth(c_green,10),smooth(c_green_std,10),'lineprops','c');
    catch
        disp(['Skipping mouse ' mouseID ' graph, not enough cells?']);
    end
end
end