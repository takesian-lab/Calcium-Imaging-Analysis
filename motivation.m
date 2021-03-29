
function [] = motivation()

ImgFolder = '\\apollo\research\ENT\Takesian Lab\Christine\0_Motivations\';


% find all the names of the image files in motivations folder
PathInfo = dir(ImgFolder);
Mask = ismember({PathInfo.name}, {'.', '..'});
PathInfo(Mask) = [];   %get rid of . and .. directories
allImg = {PathInfo.name};


% organize all the image file names into a cell array
PNGfile = dir([ImgFolder '/*.png' ]); % to find the text file with data
JPEGfile = dir([ImgFolder '/*.jpeg' ]);
JPGfile = dir([ImgFolder '/*.jpg' ]);
allName = {PNGfile.name,JPEGfile.name,JPGfile.name};

% draw a random image
randomImg = allName{randi(length(allName))};

figure
imshow([ImgFolder,randomImg]);