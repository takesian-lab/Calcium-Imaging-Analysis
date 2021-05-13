function [aviNum, aviFrames] = tosca_find_video_frames(TL, TrialNum, StateNum, Pre, Post)


if Pre > 0
   aviNum = TL.trials{TrialNum}.states(StateNum-1).aviNum(end-Pre:end);
   aviFrames = TL.trials{TrialNum}.states(StateNum-1).frameInAVI(end-Pre:end);
else
   aviNum = [];
   aviFrames = [];
end

aviNum = [aviNum; TL.trials{TrialNum}.states(StateNum).aviNum(1:Post)];
aviFrames = [aviFrames; TL.trials{TrialNum}.states(StateNum).frameInAVI(1:Post)];