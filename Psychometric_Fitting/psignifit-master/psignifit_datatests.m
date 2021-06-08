% TEST OF PSIGNIFIT 4%

%%
f1 = myplot;
f2 = myplot;
f3 = myplot;
f4 = myplot;

%Set up options structure
options.sigmoidName = 'norm'; %cumulative gaussian
options.expType = 'YesNo'; 
options.confP  = 0.95; %confidence level for confidence intervals
options.threshPC = 0.5; %percent correct corresponding to threshold


%Set up plotOptions structure
plotOptions.CIthresh       = true;            % plot a confidence interval at threshold


%results.fit = [threshold,width,lambda,gamma,eta]
%eta scales the variance of the beta distribution for the beta-binomial model 
%% Scenario 1: Fit only GO data
stims = [-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [2.5;4.5;4.5;10.5;10.5]; %(corrected as in Hautus 1995)
n_total = [11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f1);
subplot(2,3,1)
plotPsych(results,plotOptions);
title('Only fit gos')

figure(f2);
subplot(2,3,1);
plotMarginal(results,1,plotOptions);
title('Only fit gos')
%% Scenario 2: Include NOGO data (adjusted for n trials)
stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [0;2.5;4.5;4.5;10.5;10.5]; %(corrected as in Hautus 1995)
n_total = [11;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f1);
subplot(2,3,2)
plotPsych(results,plotOptions);
title('Include adjusted nogo')

figure(f2);
subplot(2,3,2);
plotMarginal(results,1,plotOptions);
title('Include adjusted nogo')
%% Scenario 3: Include NOGO data (unadjusted for n trials)
stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [0;2.5;4.5;4.5;10.5;10.5]; %(corrected as in Hautus 1995)
n_total = [200;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f1);
subplot(2,3,3)
plotPsych(results,plotOptions);
title('Include unadjusted nogo')

figure(f2);
subplot(2,3,3);
plotMarginal(results,1,plotOptions);
title('Include unadjusted nogo')
%% Scenario 4: High FA rate
stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [30;2.5;4.5;4.5;10.5;10.5]; %(corrected as in Hautus 1995)
n_total = [100;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f1);
subplot(2,3,4)
plotPsych(results,plotOptions);
title('High FA rate')

figure(f2);
subplot(2,3,4);
plotMarginal(results,1,plotOptions);
title('High FA rate')
%% Scenario 5: High lapse rate
stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [0;2.5;4.5;4.5;10.5;5.5]; %(corrected as in Hautus 1995)
n_total = [11;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f1);
subplot(2,3,5)
plotPsych(results,plotOptions);
title('High lapse rate')

figure(f2);
subplot(2,3,5);
plotMarginal(results,1,plotOptions);
title('High lapse rate')
%% Scenario 6: perfect hit rate for a stim value presented only once
stims = [-40;-15;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [0;1;2.5;4.5;4.5;10.5;10.5]; %(corrected as in Hautus 1995)
n_total = [200;1;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f1);
subplot(2,3,6)
plotPsych(results,plotOptions);
title('One trial')

figure(f2);
subplot(2,3,6);
plotMarginal(results,1,plotOptions);
title('One trial')


%% Scenario 7: Overdispersed data 
stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [0;11;1;2;0;4]; %(corrected as in Hautus 1995)
n_total = [200;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f3);
subplot(2,3,1)
plotPsych(results,plotOptions);
title(['eta = ',num2str(results.Fit(5))]);

figure(f4);
subplot(2,3,1);
plotMarginal(results,1,plotOptions);
title(['eta = ',num2str(results.Fit(5))]);
%% Scenario 8: Overdispersed data 2
stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_yes = [4;0;2;1;11;0]; %(corrected as in Hautus 1995)
n_total = [200;11;11;11;11;11];

data = [stims,n_yes,n_total];

results = psignifit(data,options);

figure(f3);
subplot(2,3,2)
plotPsych(results,plotOptions);
title(['eta = ',num2str(results.Fit(5))]);

figure(f4);
subplot(2,3,2);
plotMarginal(results,1,plotOptions);
title(['eta = ',num2str(results.Fit(5))]);



%% Scenario 9: Overdispersed data
CI_mat = [];
eta_mat = [];

stims = [-40;-12;-9;-6;-3;0]; %depth (dB re: 100%)
n_total = [11;11;11;11;11;11];

f5 = myplot;
f6 = myplot;

for i = 1:12
    n_yes  = randi(11,[1,5]); %draw 5 random integers between 1 and 11
    n_yes = [randi(2,1),n_yes]'; %cap at 20% false alarm rate
    data = [stims,n_yes,n_total];
    results = psignifit(data,options);
    
    CIrange = results.conf_Intervals(2)-results.conf_Intervals(1);
    eta = results.Fit(5);
    
    CI_mat = [CI_mat;CIrange];
    eta_mat = [eta_mat;eta];
    
    figure(f5);
    subplot(3,4,i)
    plotPsych(results,plotOptions);
    title(['eta = ',num2str(eta)]);
    
end

figure(f6);
subplot(2,2,1)

plot(eta_mat,CI_mat,'ks')
xlabel('eta')
ylabel('confidence range')