%% Noise estimation example for noisy series
% Generate noise free series of the Logistic or Pomeau-Manneville map X of N samples,
% with a chaotic parameter lambda, starting from a random initial 
% value x in [0,1]; 

N=5000;
lambda=3.5; %lambda=3.5699456 for the onset of chaos; lambda=2.05 for Pomeau-Manneville map
x=rand(1,1);
X=Logistic_dyn(x,lambda,N,0*randn(1,N));
% X=P_M(x,lambda,N,0*randn(1,N)); % For the Pomeau-Manneville map

%Generate the Gaussian noise whose std is a percentage of the amplitude of
%the series
noise_eps1=0.02*peak2peak(X)*randn(size(X));
noise_eps2=0.05*peak2peak(X)*randn(size(X));
noise_eps3=0.1*peak2peak(X)*randn(size(X));
noise_eps4=0.2*peak2peak(X)*randn(size(X));

%Perturb the series with Dynamical noise
X_eps1=Logistic_dyn(x,lambda,N,noise_eps1); %X_eps1=P_M(x,lambda,N,noise_eps1) %for P_M maps
X_eps2=Logistic_dyn(x,lambda,N,noise_eps2); %X_eps2=P_M(x,lambda,N,noise_eps2) %for P_M maps
X_eps3=Logistic_dyn(x,lambda,N,noise_eps3); %X_eps3=P_M(x,lambda,N,noise_eps3) %for P_M maps
X_eps4=Logistic_dyn(x,lambda,N,noise_eps4); %X_eps4=P_M(x,lambda,N,noise_eps4) %for P_M maps

%Estimate the std series noise
[absolute_noise_eps1,perc_noise_eps1,Noise_derivative_eps1,Apen_eps1,tgrid_eps1]=Noise_evaluation_fit(X_eps1,2,0.001,peak2peak(X_eps1));
[estimated_std_noise1,actual_std_noise1]=deal(absolute_noise_eps1,std(noise_eps1))

[absolute_noise_eps2,perc_noise_eps2,Noise_derivative_eps2,Apen_eps2,tgrid_eps2]=Noise_evaluation_fit(X_eps2,2,0.001,peak2peak(X_eps2));
[estimated_std_noise2,actual_std_noise2]=deal(absolute_noise_eps2,std(noise_eps2))

[absolute_noise_eps3,perc_noise_eps3,Noise_derivative_eps3,Apen_eps3,tgrid_eps3]=Noise_evaluation_fit(X_eps3,2,0.001,peak2peak(X_eps3));
[estimated_std_noise3,actual_std_noise3]=deal(absolute_noise_eps3,std(noise_eps3))

[absolute_noise_eps4,perc_noise_eps4,Noise_derivative_eps4,Apen_eps4,tgrid_eps4]=Noise_evaluation_fit(X_eps4,2,0.001,peak2peak(X_eps4));
[estimated_std_noise4,actual_std_noise4]=deal(absolute_noise_eps4,std(noise_eps4))
