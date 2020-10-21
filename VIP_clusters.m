% VIP Clusters
%this was written to cluster VIP/Ndnf cell activity for Carolyn and
%Maryse's project. Written by CGS on October 16, 2020. 

%the data loaded is the "meanval" from sort_extracted_data.

%% load the matched data
datapath =('\\apollo\research\ENT\Takesian Lab\Carolyn\2P Imaging data\VIPvsNDNF_response_stimuli_study\APAN 2020');
cd(datapath)
load('MatchedData.mat');

NDNFidx = 1:336;
VIPidx = 337:607;



%% Trying out different cluster #'s to see what might be best
tryclusters = 8; %how many clusters do you want to try
   for tt = 1:tryclusters;

clusters = tt;
idx = kmeans(gcamp_AUC(VIPidx,:),clusters); %index of clusters

%pull out clustered rows
for i = 1:clusters
    krow{i} = find(idx==i) %which rows in meanval == each cluster
end
%%
%restructure data to order by cluster
count=1;

    for j= 1:clusters %loop through cluster
        clustrow = krow{j}
        for k = 1:length(clustrow)
            clustmap(count,:) = gcamp_AUC(clustrow(k),:); %this will be the graph of the now grouped data
            count=count+1;
        end
        clusterlength(j) = length(clustrow); %I am making lines where each cluster ends
        if j == 1
            lines(j) = clusterlength(j);
        else  lines(j) = clusterlength(j)+lines(j-1); %cluster end line location
        end
        
    end
    
   figure;
   clims = [-10 20]
   imagesc(clustmap, clims); hold on
   h= hline(lines,'-m');
   title(['Number of clusters =  ',num2str(clusters)])
   colorbar

   end
