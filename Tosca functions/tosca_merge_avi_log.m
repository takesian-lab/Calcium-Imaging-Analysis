function TL = tosca_merge_avi_log(TL, AVI)

for k = 1:length(TL.trials)
   TL.trials{k}.frameInAVI = AVI(k).frameInAVI;
   TL.trials{k}.aviNum = AVI(k).aviNum;
   TL.trials{k}.frames = AVI(k).frames;
   TL.trials{k}.tframe = AVI(k).tframe;
   
   tframe_re_tosca = AVI(k).tframe - AVI(k).tframe(1) + TL.trials{k}.start;

   
   states = TL.trials{k}.states;
   TL.trials{k} = rmfield(TL.trials{k}, 'states');

   for ks = 1:length(states)
      st = states(ks);

      idx = tframe_re_tosca >= st.start & tframe_re_tosca < st.stop;
      
      st.frameInAVI = AVI(k).frameInAVI(idx);
      st.aviNum = AVI(k).aviNum(idx);
      st.frames = AVI(k).frames(idx);
      st.tframe = AVI(k).tframe(idx);
      
      TL.trials{k}.states(ks) = st;
   end
end

