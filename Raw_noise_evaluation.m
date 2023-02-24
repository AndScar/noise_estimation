function [Absolute_Noise,Percentage_Noise,Apen,tgrid]=Raw_noise_evaluation(X,dim,ris,range)
%The function returns a raw estimation of the noise level in the time
%series.
%The function takes as inputs:
%X - the monodimensional time series;
%dim - the embedding dimension;
%ris - resolution of the tolerance levels with respect to the range chosen;
%range - interval where the tolerance varies through resolution-dependent
%steps.

%The outputs are:
%Absolute_Noise - the level of found noise;
%Percentage_Noise - the percentage of the level of noise found wrt the
%chosen range;
%Apen - Data vector of Apen values along the grid. The length is of floor(1/resolution); 
%tgrid - Grid of tolerance values which span from 0 to range with an increment of 
%resolution x range. 

%to plot the ApEn profile in a logarithmic scale with respect to the real grid, type on the command window 
% semilogx(tgrid',Apen)

%to plot the ApEn profile with in a logarithmic scale respect to the rescaled grid, type on the command window 
% semilogx(tgrid'./range,Apen) or semilogx((ris:ris:1),Apen)



%if inputs are not enough, ris=0.001 and range=peak2peak(X) by default.
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

N=length(X);
delta=ris*range;
nbins=round(1/ris);
rgrid=linspace(delta,range,nbins);
tgrid=linspace(delta,range,nbins);
log_phi=zeros(2,length(rgrid));
z=1;
%ApEn profile computation with Histogram method
for m=dim:dim+1
A = zeros(m,N-m+1);  
for i = 1:m
        A(i,:) = X(i:N-m+i);
end

dA=zeros(N-m+1,N-m+1);
for i=1:N-m+1
    for j=i:N-m+1
        dA(i,j)=norm(A(:,i)-A(:,j),inf);
    end
end
clear A;
dA=dA+dA';

CDFA=zeros(N-m+1,nbins);
for k=1:N-m+1
[counts, bins] = histcounts(dA(k,:),(0:delta:range));
cdf = cumsum(counts);
CDFA(k,:)=cdf/(N-m+1);
end
clear dA;
log_phi(z,:)=mean(log(CDFA));
clear CDFA;
z=z+1;
end
clear X;
Apen=log_phi(1,:)-log_phi(2,:);

%Noise evaluation
%Y is the theoretical function to be compare with the ApEn 
Y=-log(rgrid)-Apen;
%Computation of the derivative
dif=diff(Y)./diff(log(rgrid));
dif=dif(1:round(length(dif)));
rgrid=rgrid(1:length(dif))';
%Smoothing procedure
dif1=smooth(dif);


% I exclude cases where the noise value found is due to dynamics - the ApEn profile 
%does not have peaks due to noise.
%The expected noise level is upposed to be less than the 30% of the peak2peak(X). 
try
if (std(Apen)<0.01)
Absolute_Noise=0;
Percentage_Noise=Absolute_Noise./range;
else
Absolute_Noise=rgrid(max(find(max(Apen)==Apen))-1+min(find(abs(dif1(min(find(max(Apen)==Apen)):round(0.2*length(dif))))==min(abs(dif1(min(find(max(Apen)==Apen)):round(0.2*length(dif))))))));
Percentage_Noise=Absolute_Noise./range;%1
if(isempty(Absolute_Noise))
Percentage_Noise=0.3;
Absolute_Noise=0.3*range;
end
end
catch
Absolute_Noise=rgrid(max(find(max(Apen)==Apen))-1+min(find(max(Apen)==Apen))+min(find(abs(dif1(min(find(max(Apen)==Apen)):round(length(dif))))==min(abs(dif1(min(find(max(Apen)==Apen)):round(length(dif))))))));
Percentage_Noise=Absolute_Noise./range;     
end
end