function [loops] = memorycheck(imageData)

newFrame=256;
nImages=size(imageData.Cropped_Imaging_Data,3);
memoryInfo = memory;
    memoryMult=0.05;
    if memoryMult*memoryInfo.MaxPossibleArrayBytes/(newFrame*newFrame*10)<nImages
        nImageWriteMax=ceil(memoryMult*memoryInfo.MaxPossibleArrayBytes/(newFrame*newFrame*10));
        nImagesOverflow = rem( nImages,nImageWriteMax);
        loops=ceil(nImages/nImageWriteMax);
        if  nImagesOverflow <= 200
            loops= loops-1;
            nImagesOverflow=nImageWriteMax+nImagesOverflow;
        end
    else
        nImageWriteMax=nImages;
        nImagesOverflow=nImages;
        loops=1;
    end
    disp(['Your machine can only contain ' num2str(memoryInfo.MaxPossibleArrayBytes/1073741824 ) ...
        ' GB in single array. Thefore dF/Fo analysis will be divided into ' num2str(loops) ' loops'])
end