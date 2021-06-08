function PsignfitDemo
%---Script that loads a appetitive Go-Nogo behavior file and runs
%psignifit function. Created 10/19/16 JDY.
clear
close all
clc

%---LOAD FILE---%
%(.e.g, load('/Users/justinyao/Sanes Lab/Experiments/PROJECTS/AMRATE/Behavior/236396/236396_28-Aug-2016.mat'))
load('/Users/Melissa/Desktop/236396_28-Aug-2016.mat')
%---Organize Data---%
DATA		=	GetData(Data,Info);
%---Plotting---%
PlotFunctions(DATA);

%---Locals---%
function DATA = GetData(Data,Info)
D						=	Data;
N						=	length(D);
RespLegend				=	Info.Bits;
Resp					=	nan(N,1);			%---1:Hit; 2:Miss; 3:CR; 4:FA---%
Rate					=	nan(N,1);		
Trial					=	nan(N,1);			%---0:Go; 1:Nogo---%
Lat						=	nan(N,1);
remind					=	nan(N,1);
for a=1:N
	dat					=	D(1,a);
	idx					=	nan(4,1);
	remind(a,1)			=	dat.Reminder;
	for b=1:4
		if( b == 1 )
			if( isfield(RespLegend, 'hit') )
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.hit);
			else
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.Hit);
			end
		end
		if( b == 2 )
			if( isfield(RespLegend, 'miss') )
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.miss);
			else
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.Miss);
			end
		end
		if( b == 3 )
			if( isfield(RespLegend, 'cr') )
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.cr);
			else
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.CR);
			end
		end
		if( b == 4 )
			if( isfield(RespLegend, 'fa') )
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.fa);
			else
				idx(b,1)	=	bitget(dat.ResponseCode,RespLegend.FA);
			end
		end	
	end
	Resp(a,1)			=	find(idx,1);
	Rate(a,1)			=	dat.AMrate;

	Trial(a,1)			=	dat.TrialType;
	Lat(a,1)			=	dat.RespLatency;	
end
DATA					=	[Trial Rate Resp Lat];
%---Remove Reminds---%
rel						=	remind == 0;
DATA					=	DATA(rel,:);

function PlotFunctions(DATA)
gel				=	DATA(:,1) == 0;
Go				=	DATA(gel,:);
Rates			=	unique(Go(:,2));
sel				=	Rates > 0;
Rates			=	Rates(sel);
NRates			=	length(Rates);
dp				=	nan(NRates,1);
pHits			=	nan(NRates,1);
Ngo				=	nan(NRates,1);
N				=	nan(NRates,1);
%---For Nogo stimuli---%
nel				=	DATA(:,1) == 1;
Nogo			=	DATA(nel,:);
fel				=	Nogo(:,3) == 4;
pfa				=	sum(fel)/length(fel);
if( pfa == 0 )
	pfa			=	1/(2*length(fel));
end
pFA				=	repmat(pfa,NRates,1);

Lat				=	cell(NRates,1);
%---For Go stimuli---%
for i=1:NRates
	rate		=	Rates(i);
	nel			=	DATA(:,2) == rate;
	N(i,1)		=	sum(nel);
	sel			=	Go(:,2) == rate;
	Ngo(i,1)	=	sum(sel);
	data		=	Go(sel,:);
	hel			=	data(:,3) == 1;
	phits		=	sum(hel)/length(hel);
	pHits(i,1)	=	phits;
	Lat(i,1)	=	{data(hel,end)};
	if( phits > 0.95  )
		phits	=	0.95;
	end
	if( phits < 0.05  )
		phits	=	0.05;
	end		
	if( phits == 1 )
		phits	=	1 - (1/(2*length(hel)));
	end
	dp(i,1)		=	calculatedprime(phits,pFA(i));
end

%---psignifit---%
[x,fitHits]		=	getpsignfit(Rates,Ngo,pHits,N);		%%%%%%%%
fa				=	repmat(pfa,length(fitHits),1);

%---Calculate d-prime---%
dpFit			=	calculatedprime(fitHits,fa);

%---Plot Hit Rates---%
subplot(1,2,1)
plot(Rates,pHits,'ko','MarkerFaceColor','k')
hold on
%---Plot fit for hits---%
plot(x,fitHits,'k-','LineWidth',2)
set(gca,'FontSize',20)
set(gca,'XTick',[0 4 64 256 1024],'XTickLabel',[0 4 64 256 1024]);xlabel('AM Rate (Hz)')
set(gca,'XScale','log');
xlim([2 1424])
ylabel('Proportion of Hits')
ylim([-0.005 1.005])
set(gca,'YTick',0:0.25:1,'YTickLabel',0:0.25:1);
axis square

%---Plot d-primes---%
subplot(1,2,2)
plot(Rates,dp,'ko','MarkerFaceColor','k')
hold on
%---Plot fit for d-primes---%
plot(x,dpFit,'k-','LineWidth',2)
set(gca,'FontSize',20)
set(gca,'XTick',[0 4 64 256 1024],'XTickLabel',[0 4 64 256 1024]);xlabel('AM Rate (Hz)')
set(gca,'XScale','log');
xlim([2 1424])
ylabel('d-prime')
ylim([-0.005 3.005])
set(gca,'YTick',0:1:3,'YTickLabel',0:1:3);
axis square

function dprime = calculatedprime(pHit,pFA)
dprime		=	nan(length(pHit),1);
for i=1:length(pHit)
	phit	=	pHit(i);
	pfa		=	pFA(i);
	zHit	=	sqrt(2)*erfinv(2*phit-1);
	zFA		=	sqrt(2)*erfinv(2*pfa-1);
% 	zHit	=	norminv(pHit,0,1);
% 	zFA		=	norminv(pFA,0,1);
	%-- Calculate d-prime --%
	dprime(i,1) = zHit - zFA;
end

function [X,fitHits] = getpsignfit(Rates,Ngo,pHits,N)
ISI		=	1000./Rates;	%For AM Rate data set%
D		=	[ISI (Ngo.*pHits) N];
Results	=	psignifit(D);
x       =	linspace(min(Results.data(:,1)),max(Results.data(:,1)),1000);
fitHits =	(1-Results.Fit(3)-Results.Fit(4))*arrayfun(@(x) Results.options.sigmoidHandle(x,Results.Fit(1),Results.Fit(2)),x)+Results.Fit(4);
X		=	1000./x;