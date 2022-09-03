%% Comparison 4: Fluxes calculated by analytical solution vs. Fluxes calcualated by numerical solution in COMSOL %%

% This MATLAB code conducts the comparison 4 analysis: the fluxes
% calculated by analytical solution when all the assumptions are satisfied 
% vs. fluxes calcualated by numerical solution when all the assumptions are satisfied, but with the
% approximated coefficients and their corresonding frequencies pairs;

% This function answers the question: What is the effect of the Fourier series fitting on the fluxes? 
%% Create a synthetic surface (2 frequencies in two direction, 4 frequencies pairs)
% Define the parameters of the synthetic surface;
A0=0; % Average surface elevation;
lambda1=10;
lambda2=30;
H_Coef1=30;
H_Coef2=0;
H_Coef3=0;
H_Coef4=45;
x=1+1000:150+1000;
y=1+1000:150+1000;
[Z_Synthetic, H_Coef, lambda_Synthetic, kx_Synthetic, ky_Synthetic]=Synthetic_Landscape(A0,lambda1,lambda2,H_Coef1,H_Coef2,H_Coef3,H_Coef4,x,y);
surf(x,y,Z_Synthetic);
%% Calculate the velocity field; real hydraulic gradient;
% fact=1;
% N_Z=25; % N_Z = Number of grid steps with depth.
Cond=.00001; % Cond = Hydraulic conductivity.
porosity=0.3; % Porosity = Effective porosity.
% Calculate the velocity field of the original synthetic surface;

av=mean(mean(Z_Synthetic));
dp=10;
W_top=zeros(150,150);
for x=1001:1150
    for y=1001:1150        
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx_Synthetic,ky_Synthetic,H_Coef',Cond,porosity); % solve the fluxes point by point;
W_top(x-1000,y-1000)=W;
    end
end
%% Solve spectral problem
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=32; 
% N=40;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(Z_Synthetic); % Get the x coordinate and y coordinate size of the cropped area of image;
X = 1:nc;
Y = 1:nr; 
X=X+1000;
Y=Y+1000;
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Z_Synthetic,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
%% Calculate the velocity field; Fourier Fitting hydraulic gradient;
% fact=1;
% N_Z=25; % N_Z = Number of grid steps with depth.
% Cond=.00001; % Cond = Hydraulic conductivity.
% porosity=0.3; % Porosity = Effective porosity.
% Calculate the velocity field of the original synthetic surface;
% av=mean(mean(Z_Synthetic));
dp=10;
W_top_FS=zeros(150,150);
for x=1001:1150
    for y=1001:1150        
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat',Cond,porosity); % solve the fluxes point by point;
W_top_FS(x-1000,y-1000)=W;
    end
end
%% Present the Comparison plot;
figure(10)

subplot(1,2,2)
contourf(W_top_FS)
colorbar
title('Fluxes of Fourier Series Fitting Surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,1)
Top_Velocity_Field_True=W_top;
contourf((Top_Velocity_Field_True));
colorbar
title('Fluxes of Real Synthetic Surface','FontSize',40)
ax = gca;
ax.FontSize = 18;
