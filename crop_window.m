
function [imageData,parameters] = crop_window(setup,parameters)

for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    date=setup.expt_date{(a)};
    Imaging_Block=setup.BOT_maps(a,:)
    
    for i=1:length(setup.BOT_maps)
        BOT_number = num2str(setup.BOT_maps(i));
        BOT_number
        folder = sprintf([setup.path_name setup.username '/' mouseID '/' date '/' BOT_number]); %direct to BOT data folder
        cd(folder);
        [folder '/*.ome.tif']
        d = dir([folder '/*.ome.tif']);%extract tiffs
        
        for k=1:length(d);
%             
            num_BOT = sprintf('%06d',((setup.BOT_start-1)+k))
            num_BOT
            image = imread(['BOT_' mouseID '_noiseburst-00' BOT_number '_Cycle00001_Ch2_' num_BOT '.ome.tif']);
            filtered_image = imresize(image,0.5);%reduce to 256x256
            Full_Tile_Matrix(:,:,k) = filtered_image;
        
        end
        
        % MEAN IMAGE: Generate and plot the mean image 
        figure;
        imageData.Full_Tile_Mean = mat2gray(mean(Full_Tile_Matrix,3));
        imshow(imageData.Full_Tile_Mean);

        title(sprintf('Mean Window Tile %d',BOT_number));
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