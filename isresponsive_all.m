PC_name = getenv('computername');

switch PC_name
    case 'RD0366' %Maryse
        info_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study';
        compiled_blocks_path = 'D:/Data/2p/VIPvsNDNF_response_stimuli_study/CompiledBlocks';
        info_filename = 'Info';
    case 'RD0332' %Carolyn
        info_path = 'D:\2P analysis\2P local data\Carolyn';
        compiled_blocks_path = 'D:\2P analysis\2P local data\Carolyn\analyzed\Daily Imaging';
        info_filename = 'Info';
    case 'RD0386' %Wisam
        % INSERT PATHS HERE
        info_filename = 'Info';
    otherwise
        disp('Computer does not match known users')
end
%%
cd(info_path)
Info = importfile(info_filename);
block_path = compiled_blocks_path;
cd(block_path)
allfiles=dir('*Block*');
%%
bl=0;
b=0;
for i = 1:size(allfiles,1)
    b=b+1;
    name = ({allfiles(b).name});
    load(name{1});
        bl=bl+1;
        block_number = sprintf('%03d',bl);
        allData.(['block' block_number]) = block; 
end
numBlocks = bl;

%%
for i = 1:numBlocks
    block_number = sprintf('%03d',bl);
    stim_type(i) = allData.(['block' block_number]).setup.stim_protocol;
    FOV (i) = allData.(['block' block_number]).setup.FOV;
    dates (i) = allData.(['block' block_number]).setup.expt_date;
    Mousename (i) = allData.(['block' block_number]).setup.mousename;
end
a_mice = unique(Mousename);
a_dates = unique(dates);
a_FOV = unique(FOV);
a_stim = unique(stim_type);

%%
for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    
    a_green = [];
    SEM_a_green = [];
    b_green = [];
    c_green = [];
    d_green = [];
    e_green = [];
    
    %go through all of the cells
    for i=1:size(data.([mouseID]).traces_G,1);
        a_green(i,:) = squeeze(mean(data.([mouseID]).traces_G(i,:,:), 2)); %mean of all trials for each cell
        SEM_a_green(i,:) = std(data.([mouseID]).traces_G(i,:,:))./sqrt(size(data.([mouseID]).traces_G(i,:,:),2));
        
        b_green(i,:) = mean(data.([mouseID]).response(i,:), 2); %average means responses (sound to 2s post sound) across trials for each cell
        c_green(i,:) = mean(data.([mouseID]).peak_G(i,:), 2); %average peak response
        d_green(i,:) = mean(data.([mouseID]).peak_maxG(i,:), 2);%average around max peak
        e_green(i,:) = mean(data.([mouseID]).peak_minG(i,:), 2); %average around negative peak
        
        %determine whether each cell is responsive (defined as average response
        %more than defined STDS above baseline)
        f = mean(data.([mouseID]).std_baseline(i,:), 2); %average baseline STD across trials for each cell
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
