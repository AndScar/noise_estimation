function [absolute_noise,perc_noise,Noise_derivative,Apen,tgrid]=Noise_evaluation_fit(X,dim,ris,range)
%   The function computes the accurate noise level in a time series through
%   non-linear fitting.
%   The function takes as inputs:
%   X - the monodimensional time series;
%   dim - the embedding dimension;
%   ris - resolution of the tolerance levels with respect to the range chosen;
%   range - interval where the tolerance varies through resolution-dependent
%   steps.

%   The outputs are:
%   absolute_noise - the level of found noise with nonlinear fitting;
%   perc_noise - the percentage of the level of noise with nonlinear fitting, found wrt the
%   chosen range;
%   Noise_derivative - the raw estimation given by the derivative
%   comparison.
%   Apen - Data vector of Apen values along the grid. The length is of floor(1/resolution); 
%   tgrid - Grid of tolerance values which span from 0 to range with an increment of 
%   resolution x range. 

%   to plot the ApEn profile in a logarithmic scale with respect to the real grid, type on the command window 
%   semilogx(tgrid',Apen)

%   to plot the ApEn profile with in a logarithmic scale respect to the rescaled grid, type on the command window 
%   semilogx(tgrid'./range,Apen) or semilogx((ris:ris:1),Apen)



%   if inputs are not enough, ris=0.001 and range=peak2peak(X) by default.
if(nargin<2)
   disp("Not enough inputs");
   return;
end
if(nargin==2)
    ris=0.001;
    range=peak2peak(X);
end
if(nargin==3)
    range=peak2peak(X);
end

% Raw estimation of the noise level with derivative comparision
[Abs_noise,Noise_derivative,Apen,tgrid]=Raw_noise_evaluation(X,dim,ris,range);
Apen=smooth(Apen);
scala=(ris:ris:1);
q=floor(Noise_derivative/ris);
scala=scala*range;
y=@(r,sigma) -log(r./(sigma*sqrt(pi)));
absolute_noise=[];
perc_noise=[];
%defining the interval in which we perform nonlinear fitting,
% namely [x_max,x_raw]where: 
%x_max - is the point for which Apen is maximum;
%x_raw - is the point where the Apen and -log r are closer.
o=find(max(Apen)==Apen);
p_min=o;

if(find(max(Apen)==Apen)>=q)
    p_min=o; 
end 
 p_max=q;

% theorethical function to be fitted
fun = @(t) -log(scala(p_min:p_max)'./(t*sqrt(pi)))-Apen(p_min:p_max); 
 x0 =Noise_derivative*range;
%compute noise estimation with starting point equal to the raw noise level
%found with the raw estimation method in the interval [x_max,x_raw]
try
x = lsqnonlin(fun,x0,0,range);
absolute_noise=x;
perc_noise(end+1)=x/range;
catch
    % if the previous initial condition is not suitable 
    % initial condition is changed
    try
    x0 = 0.08*Noise_derivative;
    x = lsqnonlin(fun,x0,0,range);
    absolute_noise=x;
     perc_noise(end+1)=x/range;
    catch
        %if nonlinear estimation is not achieved for badly conditioned
        %series or for general numerical issues, take the raw estimation as
        %the final one.
      absolute_noise=Abs_noise;
     perc_noise(end+1)=Noise_derivative;
     end
end
end
