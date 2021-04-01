% Cluster Analysis
%this was written to cluster VIP/Ndnf cell activity for Carolyn and
%Maryse's project. Written by CGS on October 16, 2020 and revised by AET on October 20-25,2020. 

%% magic numbers and setup
stimTypes = {'NoiseITI','RF','SAM','SAMfreq','FM'}; 
data_type = 0; %0 for Calcium (df/f) raster or 1 for Deconvolved Spike raster 
datapath =('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020\ExtractedData CGS\Inactive');

%% load the traces of the matched data 
cd(datapath)
load('MatchedData_all.mat');

%% build a matrix with responses to all stimuli (All_Stim) 
 %1st col is stim number, 2nd col is cell number, next col are time traces
 %for each cell's response to a stimulus

A = [];

for stim=1:length(stimTypes) %make a matrix all_stim_mat that excludes empty fields
    s = stimTypes{1,stim}
    name = getfield(Match,s) %get access to structure Match.stimType (ie Match.FM)
    if data_type ==0 %determines whether data used is calcium or spike raster
        type = 'Calcium_Raster';
    else
        type = 'Spikes_Raster';
    end
    
    for cell = 1:size(name.(type),1) %number of cells (rows) for each stimulus type
        B(cell)= size(name.(type){cell},2); %number of frames for each response
        if B(cell)>3 %eliminate empty data (NaNs)
            all_stim_mat = [stim cell name.(type){cell}]'; 
            A = [A all_stim_mat];   
        end
    end
end

All_Stim = A(1:size(A,1)-1,1:end)'; %transpose and eliminate last row of trace that is always NaN

%eliminate rows with NaN or Inf
[rows, columns] = find(isnan(All_Stim)| isinf(All_Stim)); %some rows had random NaN or Inf @ end of trace, figure out why but for now eliminated
rows = unique(rows);
All_Stim(rows,:)=[];

%% use PCA to describe responses to all stimuli (reduce dimensions of response traces from 76 frames to 3-7)

%z-score data across cells such that each cell's mean = 0 and SD = 1
z_score_B = [];
cell_list = unique(All_Stim(:,2)) %list of unique cell numbers
for cells = 1:length(cell_list)
    cell_num = cell_list(cells);
    Idx_cells = find(All_Stim(:,2)==cell_num);
    z_score_cell = zscore(All_Stim(Idx_cells,3:end),0,'all'); %z scored across cells
    z_score_A = [All_Stim(Idx_cells,1) All_Stim(Idx_cells,2) z_score_cell]
    z_score_B = [z_score_B; z_score_A];
end

%plot average df/f and z scores across all stimuli
figure;
plot(mean(All_Stim(:,3:end),1)); %plot mean raw df/f of all responses
title('mean df/f activity averaged across stimuli');
ylabel('frames');
xlabel('df_over_f');

figure;
plot(mean(z_score_B(:,3:end),1),'r'); %plot mean z scores
title('mean z scored data');


%Do the PCA
[pc,score,latent,tsquared,explained,mu] = pca(z_score_B(:,3:end));
    %output of PCA: 
    %pc = principle component vectors, each column is a principle component
    %score = data in principal component space, X*V
    %latent, eigenvalues of the covariance matrix
    %tsquared = Hotelling's t-squared statistic for each observation of matrix
    %explained = % of total variance explained by each principle component
    %mu = estimated mean of each variable in x


%Plot component scores - all the responses plotted in 2-D PCA space
    %(stimuli are color coded)
figure;
sz=15;
colours = {'b','c','g','y', 'm', 'r','k'};      %'[1, 0.5, 0]o', 
  for stimulus = 1:length(unique(z_score_B(:,1))) 
        Idx_points = find(z_score_B(:,1)==stimulus); % index of responses to stimulus n
        hi = scatter(score(Idx_points,1),score(Idx_points,2),sz,colours{stimulus},'filled');
        alpha(0.4)
        hold on;
  end
  
title('Component Scores');
ylabel('PC 1');
xlabel('PC 2');
legend(stimTypes);  

%Plot component scores  - all the responses plotted in 3-D PCA space
    %(stimuli are color coded)
figure;
sz=15;
colours = {'b','c','g','y', 'm', 'r','k'};      %'[1, 0.5, 0]o', 
  for stimulus = 1:length(unique(z_score_B(:,1))) 
        Idx_points = find(z_score_B(:,1)==stimulus); % index of responses to stimulus n
        hi = scatter3(score(Idx_points,1),...
            score(Idx_points,2),...
            score(Idx_points,3),sz,colours{stimulus},'filled');
        alpha(0.4)
        hold on;
  end
  
title('Component Scores');
ylabel('PC 1');
xlabel('PC 2');
zlabel('PC 3');
legend(stimTypes);  

%scree plot - plot of eigenvalues against corresponding PCs
    %(used to determine the number of PCs to retain to explain data)
figure;
plot(1:10,latent(1:10),'o')
title('Scree Plot');
xlabel('Principle Component');
ylabel('Eigenvalues');

%variance plot - fraction of total variance as explained by each PC
    %(obtained by taking each eigenvalue/sum of all eigenvalues
figure;
plot(1:10,explained(1:10),'o')
title('Variance explained');
xlabel('Principle Component');
ylabel('Proportion');

figure;
plot(pc(:,1:5));
title('PCs');
legend('PC1','PC2','PC3','PC4','PC5');
xlabel('frames');
ylabel('df/f');

%% K-means analysis - to find clusters across all responses

PCA_first_5 = score(:,1:5); %take first 5 PCs for K analysis

% Trying out different cluster #'s to see what might be best
tryclusters = 10; %how many clusters do you want to try
for tt = 5:tryclusters;
clusters = tt;
idx = kmeans(PCA_first_5,clusters) %index of clusters

%pull out clustered rows
for i = 1:clusters
    krow{i} = find(idx==i) %which rows in PCA_first_5 == each cluster
end


%order by activity 
  for j= 1:clusters %loop through cluster
        clustrow = krow{j}
        time = size(z_score_B,2)
        a = z_score_B(clustrow,12:time);
        m(j) = mean(mean(a,1),2);
  end
  [value,order] = sort(m,'descend'); % order from biggest to smallest

%restructure data to order by cluster 
count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{order(j)}
        for k = 1:length(clustrow)
            clustmap(count,:) = PCA_first_5 (clustrow(k),:); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines(j) = clusterlength(j);
        else  lines(j) = clusterlength(j)+lines(j-1); %cluster end line location
        end
        
    end
     
   figure;
   clims = [-1 1]
   imagesc(clustmap, clims); hold on
   h= hline(lines,'-r');
   title(['Number of clusters =  ',num2str(clusters)])
   colorbar
   xlabel('Principle Component');

% restructure data to order by cluster using raw traces
count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{order(j)}
        for k = 1:length(clustrow) % loop through each cell in cluster
            clustmap_z(count,:) = smooth(z_score_B(clustrow(k),3:end),3); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_z(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines_z(j) = clusterlength_z(j);
        else  lines_z(j) = clusterlength_z(j)+lines(j-1); %cluster end line location
        end
        
    end
   figure;
   clims = [-1 1]
   imagesc(clustmap_z, clims); hold on
   h= hline(lines_z,'-r');
   title(['Number of clusters =  ',num2str(clusters)])
   colorbar
   xlabel('Principle Component');



end
   
%% Plot component scores  - all the responses plotted in 3-D PCA space
    %(clusters are color coded)
figure;
sz=15;
colours = {'b','c',[0.3010, 0.7450, 0.9330],'g','y',[1, 0.5, 0], 'm', 'r',[0.2, 0, 0],'k'};      %'[1, 0.5, 0]o', 
    for j= 1:clusters %loop through cluster
        j
        clustrow = krow{order(j)};
      %   for k = 1:length(clustrow) % loop through each cell in cluster
            hi = scatter3(score(clustrow,1),...
            score(clustrow,2),...
            score(clustrow,3),sz,colours{j},'filled');
        alpha(0.4);
        hold on;
       %  end
        
  end
  
%title('Component Scores');
ylabel('PC 1');
xlabel('PC 2');
zlabel('PC 3');
legend(arrayfun(@num2str, 1:clusters, 'UniformOutput', 0));
set(gcf,'color','w');
hold on;
 
% rotate response

OptionZ.FrameRate=40;OptionZ.Duration=7;OptionZ.Periodic=true; 
%CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'PCs',OptionZ)
CaptureFigVid([0,0;90,90], 'PCs',OptionZ)
% 
% get(hi, 'children'); 
% rotate(hi, [0, 0, 1], 270)
%  FigHandle = gcf();
%  for Frame = 1:30
%   Frames(Frame) = getframe(FigHandle,...
%    [0, 0, 15 * Frame, 15 * Frame]);
%  end
%  clf
%  movie(FigHandle, Frames, 1, 5);
% % movie2avi(M,'WaveMovie.avi');


%% K-means analysis - to find clusters across all cells

%each cell can be described as its responses to the 5 stimuli 
    %5 PCs per stimuli
    
%lets start with cells in which we have responses to all 5

PCA_first_5 = score(:,1:5); %take first 5 PCs for K analysis
cell_list = unique(z_score_B(:,2)); %list of unique cell numbers
time = size(z_score_B(:,3:end),2)

% score_c = zeros(size(cell_list,1),length(unique(z_score_B(:,1)))*5);
% %cells_matched = zeros(size(cell_list,1),6);
% z_c = zeros(size(cell_list,1),time*5);
% raw_c = zeros(size(cell_list,1),time*5);


Score_matched = [];
z_score_matched = [];
raw_matched = [];

cont_score_matched = [];
cont_z_score_matched = [];
cont_raw_matched = [];
   

for cells = 1:length(cell_list) %loop through cells
    cell_num = cell_list(cells);
    Idx_cells_a = find(z_score_B(:,2)== cell_num) %find cells
    Idx_cells_b = find(B(:,2)== cell_num) %find cells
    
    if length(Idx_cells_a) == length(unique(z_score_B(:,1)))  
        score_b = [];
        z_b = [];
        raw_b = [];
        
        score_b_cont = [];
        z_b_cont = [];
        raw_b_cont = [];
        

        for stimulus = 1:length(unique(z_score_B(:,1))) 
            Idx_stim_a = find(z_score_B(:,1)==stimulus) % index of responses to stimulus n
            Idx_stim_b = find(B(:,1)==stimulus) % index of responses to stimulus n
            
            [idx_a loc_a] = ismember(Idx_cells_a,Idx_stim_a); %index of z score vector
            loc_a = loc_a(any(loc_a,2),:); 
            
            [idx_b loc_b] = ismember(Idx_cells_b,Idx_stim_b); %index of raw data
            loc_b = loc_b(any(loc_b,2),:); 
            
            score_a = score(Idx_stim_a(loc_a),1:5); 
            z_a = z_score_B(Idx_stim_a(loc_a),3:end); 
            raw_a = B(Idx_stim_b(loc_b),3:end); 
                       
            score_b = [score_b; score_a];
            z_b = [z_b; z_a];
            raw_b = [raw_b; raw_a];
            
            score_b_cont = [score_b_cont score_a];
            z_b_cont = [z_b_cont z_a];
            raw_b_cont = [raw_b_cont raw_a];
         
        end
        
         Score_matched = [Score_matched; score_b];
         z_score_matched = [z_score_matched; z_b];
         raw_matched = [raw_matched; raw_b];
         
         
         cont_score_matched = [cont_score_matched; score_b_cont];
         cont_z_score_matched = [cont_z_score_matched; z_b_cont];
         cont_raw_matched = [cont_raw_matched; raw_b_cont];
         
   
    end
   
        
end

%each cell's response is characterie by 25 PCs - 5 from each stimulus 
%rows are cells, columns are PCs
% = score_c(any(score_c,2),:); 

%z-score of each cell in same order as cells in B_matched
%rows are cells, columns are z-score data across 5 stimuli (71 pts each)
%z_score_matched = z_score_matched(any(z_score_matched,2),:); 

%index of each cells response to match to original cell
%
%cells_matched = cells_matched(any(cells_matched,2),:);

%raw_matched = raw_matched(any(raw_matched,2),:);

% Trying out different cluster #'s to see what might be best
tryclusters = 10; %how many clusters do you want to try
for tt = 5:tryclusters;
clusters = tt;
idx = kmeans(cont_score_matched,clusters) %index of clusters

%pull out clustered rows
for i = 1:clusters
    krow{i} = find(idx==i) %which rows in PCA_first_5 == each cluster
end

time = size(z_score_matched,2)
%order by activity 
activity_mat =[]
for stimulus = 1:length(unique(z_score_B(:,1))) 
        activity_vector = 10+((stimulus-1)*time):stimulus*time;
        activity_mat = [activity_mat activity_vector];
end
  
for j = 1:clusters %loop through cluster
    clustrow = krow{j}
    a = cont_z_score_matched(clustrow,activity_mat)';
    Q = cumsum(a);
    Q_norm = Q./repmat(max(abs(Q),[],1), size(a,1), 1);
    for i = 1:size(a,2)
        [M I] = min(abs(Q_norm(:,i)-0.4));
        Q_cross(i) = I;
    end
    Q_cross_cluster(j) = mean(Q_cross);
end
[value,order] = sort(Q_cross_cluster,'ascend'); % order from biggest to smallest


%restructure data to order by cluster 
count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{order(j)}
        for k = 1:length(clustrow)
            clustmap_cell(count,:) = smooth(cont_z_score_matched(clustrow(k),:),5); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_cell(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines_cell(j) = clusterlength_cell(j);
        else  lines_cell(j) = clusterlength_cell(j)+lines_cell(j-1); %cluster end line location
        end
        
    end
     
   figure;
   clims = [-1 1]
   imagesc(clustmap_cell, clims); hold on
   h= hline(lines_cell,'-m');
   title(['Number of clusters =  ',num2str(clusters)])
   colorbar
   xlabel('Principle Component');
clear lines_cell
clear clusterlength_cell

% restructure data to order by cluster using z scored traces - all cells
% combined

count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{order(j)}
        for k = 1:length(clustrow) % loop through each cell in cluster
            clustmap_z_cell(count,:) = smooth(cont_z_score_matched(clustrow(k),:),3); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_z_cell(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines_z_cell(j) = clusterlength_z_cell(j);
        else  lines_z_cell(j) = clusterlength_z_cell(j)+lines_z_cell(j-1); %cluster end line location
        end
        
    end
   figure;
   clims = [-1 1]
   imagesc(clustmap_z_cell, clims); hold on
   h= hline(lines_z_cell,'-m');
   title(['Number of clusters =  ',num2str(clusters)])
   colorbar
   xlabel('Principle Component');
clear lines_z_cell
clear clusterlength_z_cell

% restructure data to order by cluster using calcium traces VIP versus NDNF  
idx_NDNF = 1:196; %need to change this, quick and dirty way to separate based on cell_matched 1:336 cells
idx_VIP = 197:294;
z_score_matched_NDNF = cont_z_score_matched(idx_NDNF,:);
z_score_matched_VIP = cont_z_score_matched(idx_NDNF,:);

%plot NDNF neurons
count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{order(j)}
        clustrow = clustrow(clustrow<197)
        for k = 1:length(clustrow) % loop through each cell in cluster
            clustmap_NDNF(count,:) = smooth(cont_z_score_matched(clustrow(k),:),3); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_NDNF(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines_NDNF(j) = clusterlength_NDNF(j);
        else  lines_NDNF(j) = clusterlength_NDNF(j)+lines_NDNF(j-1); %cluster end line location
        end
        
    end
   figure;
   clims = [-1 1]
   imagesc(clustmap_NDNF, clims); hold on
   h= hline(lines_NDNF,'-m');
   title(['NDNF Number of clusters =  ',num2str(clusters)])
   colorbar
   xlabel('Principle Component');
clear lines_NDNF 
clear clusterlength_NDNF

%plot VIP neurons
count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{order(j)}
        clustrow = clustrow(clustrow>196)
        for k = 1:length(clustrow) % loop through each cell in cluster
            clustmap_VIP(count,:) = smooth(cont_z_score_matched(clustrow(k),:),3); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_VIP(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines_VIP(j) = clusterlength_VIP(j);
        else  lines_VIP(j) = clusterlength_VIP(j)+lines_VIP(j-1); %cluster end line location
        end
    end
    
   figure;
   clims = [-1 1]
   imagesc(clustmap_VIP, clims); hold on
   h= hline(lines_VIP,'-m');
   title(['VIP Number of clusters =  ',num2str(clusters)])
   colorbar
   xlabel('Z scored traces');
clear lines_VIP
clear clusterlength_VIP
end

%VIP and NDNF cells in each cluster
%NDNF
 for j= 1:clusters %loop through cluster       
         clustrow = krow{order(j)}
        clustrow = clustrow(clustrow<197)
        for k = 1:length(clustrow) % loop through each cell in cluster
            clustmap_NDNF(count,:) = smooth(cont_z_score_matched(clustrow(k),:),3); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_NDNF(j) = length(clustrow); %I am making lines where each cluster ends
 end


%VIP
 for j= 1:clusters %loop through cluster       
        clustrow = krow{order(j)}
        clustrow = clustrow(clustrow>196)
        for k = 1:length(clustrow) % loop through each cell in cluster
            clustmap_VIP(count,:) = smooth(cont_z_score_matched(clustrow(k),:),3); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength_VIP(j) = length(clustrow); %I am making lines where each cluster ends
 end  
        
proportions_VIP_clusters = clusterlength_VIP./sum(clusterlength_VIP)*100;
proportions_NDNF_clusters = clusterlength_NDNF./sum(clusterlength_NDNF)*100;

%% Plot component scores  - all the neuron clusters plotted in 3-D PCA space
    %(clusters are color coded)
figure;
sz=15;
colours = {'b','c',[0.3010, 0.7450, 0.9330],'g','y',[1, 0.5, 0], 'm', 'r',[0.2, 0.2, 0.2],'k'};      %'[1, 0.5, 0]o', 
    for j= 1:clusters %loop through cluster
        j
        clustrow = krow{order(j)};
         for k = 1:length(clustrow) % loop through each cell in cluster
            hi = scatter3(score(clustrow(k),1),...
            score(clustrow(k),2),...
            score(clustrow(k),3),sz,colours{j},'filled');
        alpha(0.4);
        hold on;
         end
        
  end
  
title('Component Scores');
ylabel('PC 1');
xlabel('PC 2');
zlabel('PC 3');
legend(stimTypes);  


%% plot examples of NDNF cells in each category

for j= 5; %1:clusters %loop through cluster       
    clustrow = krow{order(j)}
    clustrow = clustrow(clustrow>197)
    for k =1:length(clustrow)
    %traces = cells_matched(clustrow(rows),2:6);
 %   figure;
 figure;
           time = size(z_score_B,2)-2;
    %       traces = cont_raw_matched(clustrow(k),:);
            for n = 1:5
              subplot(1,5,n)
              trace_vector = 1+(n-1)*time:n*time;
               
           %   plot(raw_matched(traces(n),3:end))
            trace = smooth(cont_raw_matched(clustrow(k),trace_vector),5);
              h = plot(trace);
              set(h,'LineWidth',2);
%           % plot(cells_matched(clustrow(k),3:end));
%            title(['Ndnf cell  =  ',num2str(k)])
            xlabel('df/f');
            ylim([-0.5 1.5])
%            end 
%           figure;
            end
      
        %   plot(z_score_matched(clustrow(k),:),'r')
       %    hold on;
        %   plot(raw_matched(clustrow(k),:)*5)
    end
end

    