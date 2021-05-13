function Frames = tosca_read_video_frames(TL, TrialNum, StateNum, Pre, Post)

[aviNum, aviFrames] = tosca_find_video_frames(TL, TrialNum, StateNum, Pre, Post);

[~, fn] = fileparts(TL.filename);

Frames = {};

curAviNum = -1;
for k = 1:length(aviNum)
   if aviNum(k) ~= curAviNum
      aviFileName = sprintf('%s.%03d.avi', fn, aviNum(k));
      v = VideoReader(fullfile(TL.aviFolder, aviFileName));
      curAviNum = aviNum(k);
   end
   Frames{k} = read(v, aviFrames(k));
end