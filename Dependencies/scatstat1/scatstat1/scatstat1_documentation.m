%% |scatstat1| documentation 
% The |scatstat1| function returns statistical values of all points within a given 
% radius of each value. This is similar to taking a moving mean, but points do not
% have to be equally spaced, nor do x values need to be monotonically increasing. 
% 
%% Syntax 
% 
%  ybar = scatstat1(x,y,radius)
%  ybar = scatstat1(x,y,radius,fun) 
% 
%% Description 
% 
% |ybar = scatstat1(x,y,radius)| returns the mean of all |y| values 
% within specified |radius| at each point |x|.  
% 
% |ybar = scatstat1(x,y,radius,fun)| applies any function |fun| to |y|
% values, default |fun| is |@mean|, but could be |@median|, |@nanstd|, etc. 
% 
%% Example: 
% Get the local median of all points within 10 units of each x point. Start by
% creating some random non-sorted data with N = 10,000 points: 

N = 10000; 
x = randi(300,N,1) + 20+3*randn(N,1) ; 
y = 3*sind(x) + randn(size(x)) + 3;  
plot(x,y,'k.') 
axis tight
box off

%% 
% Now get the moving median of all points within 15 x units of each
% value: 

yb = scatstat1(x,y,15,@median); 

hold on
plot(x,yb,'bo')

%% 
% Note, this function does not interpolate or enforce any equal spacing.  Final
% values correspond to each x point, which need not be equally spaced or sorted: 

axis([167 173 1 6])

%% Author info: 
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the University 
% of Texas at Austin's Institute for Geophysics (UTIG), June 2016. 