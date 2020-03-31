function [data] = isresponsive_byStim(data,setup,std_level)
for a=1:length(setup.mousename)
    mouseID=setup.mousename{(a)}
    
    Var1List=unique(data.([mouseID]).parameters.variable1);%what is Variable#1
    Var2List=unique(data.([mouseID]).parameters.variable2);%what is Variable#w
    data.([mouseID]).parameters.Var1List=Var1List;%store in a list for later
    data.([mouseID]).parameters.Var2List=Var2List;%store in a list for later
    n1=data.([mouseID]).parameters.variable1;%pull out variable#1 in order presented in experiment
    n2=data.([mouseID]).parameters.variable2;%pull out variable#2 in order presented in experiment
    for m=1:length(Var1List)%loop through variable1 (frequency for TRF
        p=find(n1==Var1List(m));%pull out a particular stimulus (Var #1) (i.e. 4kHz)
        for q=1:length(Var2List)
            r=find(n2==Var2List(q));%pull out a particular stimulus (Var #2) (i.e. 60dB)
            [s]=intersect(p,r); %find specific stim types (i.e. 4khz, 60dB)
            data.([mouseID]).parameters.stimIDX(m,q)={s};%stim index (Var#1xVar#2, i.e. freq x level)
            
            for i=1:size(data.([mouseID]).traces_G,1);
                a_green(i,:) = mean(data.([mouseID]).traces_G(i,s,:), 2); %mean of all trials for given stim type for each cell
                SEM_a_green(i,:) = std(data.([mouseID]).traces_G(i,s,:))./sqrt(size(data.([mouseID]).traces_G(i,:,:),2));
                
                b_green(i,:) = mean(data.([mouseID]).response(i,s), 2); %average means responses (sound to 2s post sound) across trials for each cell
                c_green(i,:) = mean(data.([mouseID]).peak_G(i,s), 2); %average peak response
                d_green(i,:) = mean(data.([mouseID]).peak_maxG(i,s), 2);%average around max peak
                e_green(i,:) = mean(data.([mouseID]).peak_minG(i,s), 2); %average around negative peak
           
            f = mean(data.([mouseID]).std_baseline(i,s), 2); %average baseline STD across trials for each cell
            data.([mouseID]).parameters.isRespPosStim(m,q,i) =d_green(i,:) > std_level*mean(f) & b_green(i,:)>0; %will be 0 or 1
            data.([mouseID]).parameters.isRespNegStim(m,q,i) = e_green(i,:) < -std_level*mean(f) & b_green(i,:)<0;
            end
            

          end
 
        end
        
    end
end
