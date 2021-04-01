
function [imageData,data,parameters] = crop_window(data)
setup = data.setup;
for a=1:size(setup.mousename,1) %Mice
    for b=1:size(setup.mousename,2) %ROIs
        
        if isempty(setup.mousename{a,b})
            continue;
        end
        
        mouseID=setup.mousename{a,b};
        Imaging_Block=setup.Imaging_sets{a,b};
        date=setup.expt_date{a,:};
    
    for i=1:length(Imaging_Block)
        unique_block_name = setup.unique_block_names{a,b}(i);
        block = data.([mouseID]).([unique_block_name]);
        folder = convertStringsToChars(block.setup.block_path); %convertStringsToChars(setup.BOTpath);
        
        %how many tiffs do we need to load?
        tiffnum = length(block.timestamp);
        display(sprintf('...Loading %d tiffs to make window...',tiffnum))
        
        
%         BOT_number = num2str(setup.BOT_maps(i));
%         BOT_number
%         folder = sprintf([setup.Compiled_blocks_path '/' mouseID '/' date '/' BOT_number]); %direct to BOT data folder
        cd(folder);
%         [folder '*/*.ome.tif']
        d = dir([folder '*/*.ome.tif']);%extract tiffs
        
        num_BOT_count = 1;
        
        
        for k=1:length(d);
%             
            num_BOT = sprintf('%06d',((setup.BOT_start-1)+k));
            fig_name = strcat(block.setup.block_name, '_Cycle00001_', setup.imaging_chan, '_', num_BOT, '.ome.tif');
            fig_name = convertStringsToChars(fig_name);
            image = imread(fig_name);
            filtered_image = imresize(image,0.5);%reduce to 256x256
            Full_Tile_Matrix(:,:,k) = filtered_image;
            
            %Print update every 100 tiffs to let the user know the code is working
            num_BOT_count = num_BOT_count + 1;
            if mod(num_BOT_count,100) == 0
                disp(num_BOT_count)
            end
        
        end
        
        % MEAN IMAGE: Generate and plot the mean image 
        figure;
        imageData.Full_Tile_Mean = mat2gray(mean(Full_Tile_Matrix,3));
        imshow(imageData.Full_Tile_Mean);

        title(sprintf('Mean Window Tile'));
%         data.([setup.mousename]).(['Tile' BOT_number]).Window_Image = Full_Tile_Mean;
        
        
        % CROP IMAGE: Crop the image to the cranial window using the mean image generated above
        figure;
        h_im = imshow(imageData.Full_Tile_Mean); hold on;
        message = sprintf('Draw a circle and then hit the space bar,');
        uiwait(msgbox(message));
        h = imellipse;
        pause;
        BW = createMask(h);
        pos_window = getPosition(h); %returns [x min, y min, width, height]
        %csvwrite('[mouseID] 'window_position'']) to mouse folder
        parameters.x_min = pos_window(1)
        parameters.y_min = pos_window(2)
        parameters.x_max = pos_window(3)+parameters.x_min
        parameters.y_max = pos_window(4)+parameters.y_min
        Tile_Mean_ROI = imageData.Full_Tile_Mean;
        Tile_Mean_ROI(BW==0)=0;
        figure; imshow(Tile_Mean_ROI);title(sprintf('Masked Window'));
        pause;
        
        Tile_Mean_ROI=Tile_Mean_ROI(parameters.y_min:parameters.y_max,parameters.x_min:parameters.x_max);
        figure; imshow(Tile_Mean_ROI); title(sprintf('Cropped Window'));
        pause;
        parameters.Window_Postion = [parameters.x_min parameters.x_max parameters.y_min parameters.y_max];
        imageData.Cropped_Imaging_Data = Full_Tile_Matrix(parameters.y_min:parameters.y_max,parameters.x_min:parameters.x_max,:);
        clear Full_Tile_Matrix
        clear Cropped
        clear Tile_Mean_ROI
   
    end
end
end