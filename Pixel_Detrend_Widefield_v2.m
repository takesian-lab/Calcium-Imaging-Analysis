function  [All_Images_df_over_f] = Pixel_Detrend_Widefield_v2(Cropped,loops);
%
%  Usage: data=A1_map_tile(data, Sound_Time, Sound_Frequency,Sound_Level, frequencies,i)
%  Inputs:
%  Note that units of Fs, movinwin have to be consistent.

%% Detrend all the data within the tile

size_dim1 = size(Cropped,1);
size_dim2 = size(Cropped,2);
size_dim3=size(Cropped,3);
% All_Images_df_over_f = zeros(size(Tile_data));
factor=floor(size_dim1/loops);
factor_overflow=size_dim1-(factor*(loops-1));
%All_Images_detrended = zeros(size(Tile_data));
%All_Images_trend = zeros(size(Tile_data));
q=0;
disp(['Pre-processing to get dF/Fo... loop ' num2str(0) '/' num2str(loops)]);
tic
for ll=1:loops
    disp(['Pre-processing to get dF/Fo... loop ' num2str(ll) '/' num2str(loops)])
    toc
    if ll==1
        q=factor;
        x_start=1;
        x_end=q;
        temp_DF_F0=zeros(factor,size_dim2,size_dim3);
    elseif ll==loops
        g=q+1;
        q=q+factor_overflow;
        x_start=g;
        x_end=q;
        temp_DF_F0=zeros(factor_overflow,size_dim2,size_dim3);
    else g=q+1;
        q=q+factor;
        x_start=g;
        x_end=q;
        temp_DF_F0=zeros(factor,size_dim2,size_dim3);
    end
    
    
    x_start
    x_end
     x=x_start:x_end;
     for i=1:length(x)
         toc
         for y=1:size_dim2
             %             x(i)
             f1 = squeeze(double(Cropped(x(i),y,:)));
             locdetrend_temp = locdetrend(f1,1,[300 10]);
             temp_DF_F0(i,y,:) = locdetrend_temp./(f1-locdetrend_temp);
         end
     end
    
    loop_num=num2str(ll);
    All_Images_df_over_f.(['Tile' loop_num])(:,:,:) = temp_DF_F0;
    clear temp_DF_F0
end


end



