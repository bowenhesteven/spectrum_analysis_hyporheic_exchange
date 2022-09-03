%%% Test a synthetic surface with two parts: characteristic spatial part
%%% and noisy part using new proposed idea %%%

%% Create a synthetic surface (For characteristic part of the surface: 1 frequencies in two direction, 1 frequencies pairs,
% for noise part of the surface: red noise);

% Define the synthetic characteristic part;
A0=0; % Average surface elevation;
lambda1=10;
lambda2=30;
H_Coef1=30;
H_Coef2=0;
H_Coef3=0;
H_Coef4=0;
x=1+1000:150+1000;
y=1+1000:150+1000;
[Z_Characteristic, H_Coef, lambda_Synthetic, kx_Synthetic, ky_Synthetic]=Synthetic_Landscape(A0,lambda1,lambda2,H_Coef1,H_Coef2,H_Coef3,H_Coef4,x,y);
% surf(x,y,Z_Synthetic);
surfl(x,y,Z_Characteristic);colormap copper
shading interp

% Define the synthetic red noise part of the surface;
Z_noise = 30*rednoise(150, 150);
surfl(x,y,Z_noise);colormap copper
shading interp

% Plot the two part, respectively
subplot(2,1,1)
surfl(x,y,Z_Characteristic);colormap copper
shading interp
title('Characteristic part of the synthetic surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(2,1,2)
surfl(x,y,Z_noise);colormap copper
title('Noise part of the synthetic surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

% Plot the total synthetic surface;
Z_Synthetic=Z_Characteristic+Z_noise;
surfl(x,y,Z_Synthetic);colormap copper
shading interp
title('The whole synthetic surface','FontSize', 40)
ax = gca;
ax.FontSize = 18; 
%% Numerically calculate the boundary flux of characteristic part, noise part and the whole synthetic surface part
% Import the top flux field into MATLAB;
% (1) Characteristic spatial scale

[X,Y] = meshgrid(x,y);
A=[X(:),Y(:)]; % create coordinate matrix;
dlmwrite('Coordinate.txt', A, 'delimiter', '\t'); % export the 2D coordinate measurements;

% Export the synthetic surface;
% DEM_Synthetic=[X1(:),Y1(:),Topo_whole(:)];
% dlmwrite('DEM_Synthetic.txt', DEM_Synthetic, 'delimiter', '\t'); % export the 2D coordinate measurements
% 
DEM_Characteristic=[X(:),flip(Y(:)),Z_Characteristic(:)];
dlmwrite('DEM_Characteristic.txt', DEM_Characteristic, 'delimiter', '\t'); % export the 2D coordinate measurements

filename = 'Boundary_Flux_Characteristic.txt';
M_Characteristic = dlmread(filename,'',9,0); 

Top_Velocity_Field_Comsol_Characteristic=zeros(150,150);
X_Comsol_Characteristic=zeros(150,150);
Y_Comsol_Characteristic=zeros(150,150);
for i=1:150
    Top_Velocity_Field_Comsol_Characteristic(:,i)=flip(M_Characteristic((i-1)*150+1:150*i,3));   
    X_Comsol_Characteristic(:,i)=M_Characteristic((i-1)*150+1:150*i,1);
    Y_Comsol_Characteristic(:,i)=flip(M_Characteristic((i-1)*150+1:150*i,2));
end
rho=998; % The density of water is 998 kg/m^3;
Boundary_Flux_Numerical_Characteristic=Top_Velocity_Field_Comsol_Characteristic/rho;

% Noise part
% Export the noise part of the synthetic surface;
 DEM_Noise=[X(:),flip(Y(:)),Z_noise(:)];
 dlmwrite('DEM_Noise.txt', DEM_Noise, 'delimiter', '\t'); % export the 2D coordinate measurements
 
filename = 'Boundary_Flux_Noise.txt';
M_Noise = dlmread(filename,'',9,0); 

Top_Velocity_Field_Comsol_Noise=zeros(150,150);
X_Comsol_Noise=zeros(150,150);
Y_Comsol_Noise=zeros(150,150);
for i=1:150
    Top_Velocity_Field_Comsol_Noise(:,i)=flip(M_Noise((i-1)*150+1:150*i,3));   
    X_Comsol_Noise(:,i)=M_Noise((i-1)*150+1:150*i,1);
    Y_Comsol_Noise(:,i)=flip(M_Noise((i-1)*150+1:150*i,2));
end
rho=998; % The density of water is 998 kg/m^3;
Boundary_Flux_Numerical_Noise=Top_Velocity_Field_Comsol_Noise/rho;

% Whole part
% Export the noise part of the synthetic surface;
 DEM_Whole=[X(:),flip(Y(:)),Z_Synthetic(:)];
 dlmwrite('DEM_Whole.txt', DEM_Whole, 'delimiter', '\t'); % export the 2D coordinate measurements
 
filename = 'Boundary_Flux_whole.txt';
M_Whole = dlmread(filename,'',9,0); 

Top_Velocity_Field_Comsol_Whole=zeros(150,150);
X_Comsol_Whole=zeros(150,150);
Y_Comsol_Whole=zeros(150,150);
for i=1:150
    Top_Velocity_Field_Comsol_Whole(:,i)=flip(M_Whole((i-1)*150+1:150*i,3));   
    X_Comsol_Whole(:,i)=M_Whole((i-1)*150+1:150*i,1);
    Y_Comsol_Whole(:,i)=flip(M_Whole((i-1)*150+1:150*i,2));
end
rho=998; % The density of water is 998 kg/m^3;
Boundary_Flux_Numerical_Whole=Top_Velocity_Field_Comsol_Whole/rho;

% Plot the characteristic part flux and the noise part flux, respectively
subplot(1,2,1)
contourf(Boundary_Flux_Numerical_Characteristic);
colorbar
title('Boundary flux induced by characteristic part of the surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2);
contourf(Boundary_Flux_Numerical_Noise);
colorbar
title('Boundary flux induced by noise part of the surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

% Plot the whole part flux and the sum of characteristic part and noise
% part of the flux, respectively;

subplot(1,2,1)
contourf(Boundary_Flux_Numerical_Whole);
colorbar
title('Boundary flux of the whole surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2);
Boundary_Sum=Boundary_Flux_Numerical_Noise+Boundary_Flux_Numerical_Characteristic;
contourf(Boundary_Sum);
colorbar
title('Characteristic boundary flux + Noise boundary flux ','FontSize',40)
ax = gca;
ax.FontSize = 18; 
%% Evaluate the accuracy of two Numerical Fluxes: (1) Whole; (2) Sum of the Noise and the Characteristic part
% (1) Histogram plot evaluation
figure(10)
subplot(1,2,1)
% Top_Velocity_Field_Comsol_rho=Top_Velocity_Field_Comsol/rho;
% contourf(Top_Velocity_Field_Comsol_rho)
% colorbar
histogram(Boundary_Sum,100);
title('Boundary flux of the whole surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2)
histogram(Boundary_Flux_Numerical_Whole,100);
title('Characteristic boundary flux + Noise boundary flux','FontSize',40)
ax = gca;
ax.FontSize = 18; 

% (2) 2D spatial mean error evaluation
Dif_in=abs((Boundary_Sum-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
contour(x,y,Dif_in,'Fill','on'); % plot the spatial distribution of the mean-error
colorbar;
title('Spatial distribution of mean error','Fontsize',40);
xlabel('x (m)','FontSize',20);
ylabel('y (m)', 'FontSize',20);
ax = gca;
ax.FontSize = 18; 

% (3) Scatter plot evaluation
figure
Boundary_Flux_Numerical_Whole_vector=reshape(Boundary_Flux_Numerical_Whole,[150*150,1]);
Boundary_Sum_vector=reshape(Boundary_Sum,[150*150,1]);
scatter(Boundary_Flux_Numerical_Whole_vector,Boundary_Sum_vector);
hold on;
x = [-1.5*10^-03 2*10^-03];
y = [-1.5*10^-03 2*10^-03];
plot(x,y,'Color','red','LineStyle','--','LineWidth',3)
title('Scatter plot of boundary fluxes','Fontsize',40);
xlabel('Boundary fluxes of the whole surface ','FontSize',20);
ylabel('Characteristic boundary fluxes+Noise boundary flux', 'FontSize',20);
%% Calculate the analytical solution for boundary flux of the characteristic part

% fact=1;
% N_Z=25; % N_Z = Number of grid steps with depth.
Cond=.00001; % Cond = Hydraulic conductivity.
porosity=0.3; % Porosity = Effective porosity.
% Calculate the velocity field of the original synthetic surface;

av=mean(mean(Z_Characteristic));
dp=10;
W_top=zeros(150,150);
for x=1001:1150
    for y=1001:1150        
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx_Synthetic,ky_Synthetic,H_Coef',Cond,porosity); % solve the fluxes point by point;
W_top(151-(y-1000),x-1000)=W;
    end
end

Boundary_Flux_Analytical_Characteristic = W_top*porosity; % Results from Comsol is considered as true;

subplot(1,2,1)
contourf(Boundary_Flux_Numerical_Characteristic)
colorbar
title('Characteristic Fluxes of Analytical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2)
contourf(Boundary_Flux_Numerical_Characteristic);
colorbar
title('Characteristic Fluxes of Numerical Model','FontSize',40)
ax = gca;
ax.FontSize = 18;
%% Using Fourier fitting to calculate the analytical solution for boundary flux of the noise part
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=8; 
% N=40;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(Z_noise); % Get the x coordinate and y coordinate size of the cropped area of image;
X = 1:nc;
Y = 1:nr; 
X=X+1000;
Y=Y+1000;
Y=flip(Y);
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Z_noise,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

figure(10)
subplot(2,1,1)
FS_Surface=h;
surf(X,Y,FS_Surface);
% surfl(X,Y,h), colormap copper
shading interp
title('Fourier fitted noise part of the surface (N=128)','FontSize',40)

subplot(2,1,2)
surfl(X,Y,Z_noise);colormap copper
shading interp
title('Real noise part of the surface ','FontSize',40)

dp=10;
W_top_FS=zeros(150,150);
for x=1001:1150
    for y=1001:1150      
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat',Cond,porosity); % solve the fluxes point by point;
W_top_FS(151-(y-1000),x-1000)=W;
    end
end

figure(10)
Boundary_Flux_Analytical_Noise_FFS=W_top_FS*porosity;
subplot(1,2,2)
contourf(Boundary_Flux_Analytical_Noise_FFS )
colorbar
title('Fluxes of Fourier Fitted Noise Surface (N=8)','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,1)
contourf((Boundary_Flux_Numerical_Noise));
colorbar
title('Numerical Fluxes of Noise Surface','FontSize',40)
ax = gca;
ax.FontSize = 18;

mean_error_noise_flux=mean(mean(abs(Boundary_Flux_Analytical_Noise_FFS-Boundary_Flux_Numerical_Noise)))/abs(min(min(Boundary_Flux_Numerical_Noise))-max(max(Boundary_Flux_Numerical_Noise)));
%% Using Fourier fitting to calculate the analytical solution for boundary flux of the whole part
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=8; 
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
Y=flip(Y);
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Z_Synthetic,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

figure(10)
subplot(2,1,1)
FS_Surface=h;
surf(X,Y,FS_Surface);
% surfl(X,Y,h), colormap copper
shading interp
title('Fourier fitted whole part of the surface (N=128)','FontSize',40)

subplot(2,1,2)
surfl(X,Y,Z_Synthetic);colormap copper
shading interp
title('Real whole part of the surface ','FontSize',40)

dp=10;
W_top_FS=zeros(150,150);
for x=1001:1150
    for y=1001:1150      
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat',Cond,porosity); % solve the fluxes point by point;
W_top_FS(151-(y-1000),x-1000)=W;
    end
end

figure(10)
Boundary_Flux_Analytical_Whole_FFS=W_top_FS*porosity;
subplot(1,2,2)
contourf(Boundary_Flux_Analytical_Whole_FFS)
colorbar
title('Fluxes of Fourier Fitted Whole Surface (N=32)','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,1)
contourf((Boundary_Flux_Numerical_Whole));
colorbar
title('Numerical Fluxes of Whole Surface','FontSize',40)
ax = gca;
ax.FontSize = 18;

mean_error_whole_flux=mean(mean(abs(Boundary_Flux_Analytical_Whole_FFS-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
%% Combine the characteristic part of the fluxes and the noise part of the fluxes
Boundary_Flux_Analytical_Characteristic=Boundary_Flux_Numerical_Characteristic;
Boundary_Flux_Analytical_Sum_N8=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N8;
Boundary_Flux_Analytical_Sum_N16=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N16;
Boundary_Flux_Analytical_Sum_N32=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N32;
Boundary_Flux_Analytical_Sum_N64=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N64;
Boundary_Flux_Analytical_Sum_N128=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N128;
%% Calculate the mean error flux for each scenario
mean_error_sum_N8=mean(mean(abs(Boundary_Flux_Analytical_Sum_N8-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_sum_N16=mean(mean(abs(Boundary_Flux_Analytical_Sum_N16-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_sum_N32=mean(mean(abs(Boundary_Flux_Analytical_Sum_N32-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_sum_N64=mean(mean(abs(Boundary_Flux_Analytical_Sum_N64-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_sum_N128=mean(mean(abs(Boundary_Flux_Analytical_Sum_N128-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

mean_error_flux_sum=[mean_error_sum_N8 mean_error_sum_N16 mean_error_sum_N32 mean_error_sum_N64 mean_error_sum_N128];
% mean_error_flux_whole=[0.0797 0.08 0.0808 0.0811 15.23];
% mean_error_flux_noise=[0.0742 0.0745 0.0751 0.0756 17.47];

%% Calculate the mean error topography for each scenario

Z_Sum_N8=Z_Characteristic+Z_noise_N8;
Z_Sum_N16=Z_Characteristic+Z_noise_N16;
Z_Sum_N32=Z_Characteristic+Z_noise_N32;
Z_Sum_N64=Z_Characteristic+Z_noise_N64;
Z_Sum_N128=Z_Characteristic+Z_noise_N128;

% Calculate the mean error topography of each scenario
mean_error_toposum_N8=mean(mean(abs(Z_Sum_N8-Z_Synthetic)))/abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
mean_error_toposum_N16=mean(mean(abs(Z_Sum_N16-Z_Synthetic)))/abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
mean_error_toposum_N32=mean(mean(abs(Z_Sum_N32-Z_Synthetic)))/abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
mean_error_toposum_N64=mean(mean(abs(Z_Sum_N64-Z_Synthetic)))/abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
mean_error_toposum_N128=mean(mean(abs(Z_Sum_N128-Z_Synthetic)))/abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));

mean_error_toposum=[mean_error_toposum_N8 mean_error_toposum_N16 mean_error_toposum_N32 mean_error_toposum_N64 mean_error_toposum_N128];
% mean_error_topo_whole=[0.1164 0.1162 0.1137 0.1100 0.1032];
% mean_error_topo_noise=[0.1156 0.1154 0.1118 0.1070 0.1010];
%% Plot the difference in both 1D and 2D surface plot evaluation
% Evaluate the surface topography: noise part, whole scenario, sum scenario

% 1D evaluation
N_WL=[8 16 32 64 128];
plot(N_WL,mean_error_topo_noise,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','b',...
    'MarkerSize',10,'Color','b')
title('Influence of number of fitting frequencies','FontSize',40)
ax = gca;
ax.FontSize = 18; 
hold on

plot(N_WL,mean_error_topo_whole,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','r',...
    'MarkerSize',10,'Color','r')
title('Influence of number of fitting frequencies','FontSize',40)
ax = gca;
ax.FontSize = 18; 
hold on

plot(N_WL,mean_error_toposum,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','g',...
    'MarkerSize',10,'Color','g')
title('Influence of number of fitting frequencies','FontSize',50)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of fitting frequencies','FontSize',20)
ylabel('Mean error','FontSize',20)
legend('Fourier fitted noise part of the surface','Fourier fitted the whole part of the surface','The sum of the characteristic part and the Fourier fitted noise part of the surface')

% Evaluate the boundary flux: noise part, whole scenario, sum scenario
N_WL=[8 16 32 64];
plot(N_WL,mean_error_flux_noise(1:4),'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','b',...
    'MarkerSize',10,'Color','b')
% title('Influence of number of fitting frequencies on boundary fluxes induced by noise part','FontSize',40)
ax = gca;
ax.FontSize = 18; 
hold on

plot(N_WL,mean_error_flux_whole (1:4),'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','r',...
    'MarkerSize',10,'Color','r')
% title('Influence of number of fitting frequencies on boundary fluxes inducd by whole surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 
hold on

plot(N_WL,mean_error_flux_sum (1:4),'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','g',...
    'MarkerSize',10,'Color','g')
title('Influence of number of fitting frequencies on boundary fluxes','FontSize',50)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of fitting frequencies','FontSize',20)
ylabel('Mean error','FontSize',20)
legend('Boundary fluxes of the Fourier fitted noise part of the surface','Boundary fluxes of the Fourier fitted the whole part of the surface','Boundary fluxes of the the sum of the characteristic part and the Fourier fitted noise part of the surface')

% 2D evaluation
% Calculate the mean error topography of each scenario
TwoD_mean_error_toposum_N8=(abs(Z_Sum_N8-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_toposum_N16=(abs(Z_Sum_N16-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_toposum_N32=(abs(Z_Sum_N32-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_toposum_N64=(abs(Z_Sum_N64-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_toposum_N128=(abs(Z_Sum_N128-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));

TwoD_mean_error_topo_N8=(abs(Z_whole_N8-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_topo_N16=(abs(Z_whole_N16-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_topo_N32=(abs(Z_whole_N32-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_topo_N64=(abs(Z_whole_N64-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));
TwoD_mean_error_topo_N128=(abs(Z_whole_N128-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));


surf(TwoD_mean_error_toposum_N64 , 'FaceColor','g');
hold on;
surf(TwoD_mean_error_topo_N64, 'FaceColor','r');
zlabel('Mean error','FontSize', 20);
title('Comparison of approximated topography by two methods (N=64)', 'FontSize',20);
legend('New proposed method','Current method');
diff_mean_error_topo=TwoD_mean_error_topo_N64-TwoD_mean_error_toposum_N64;
count=0;
for i=1:150
    for j=1:1:150
        if diff_mean_error_topo(i,j)>0
            count=count+1;
        end
    end
end
percent_topo=count/(150*150);
            

% Calculate the spatial distribution of mean error flux of each scenario
TwoD_mean_error_fluxsum_N8=(abs(Boundary_Flux_Analytical_Sum_N8-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N16=(abs(Boundary_Flux_Analytical_Sum_N16-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N32=(abs(Boundary_Flux_Analytical_Sum_N32-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N64=(abs(Boundary_Flux_Analytical_Sum_N64-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N128=(abs(Boundary_Flux_Analytical_Sum_N128-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

TwoD_mean_error_flux_N8=(abs(Boundary_Flux_Analytical_Whole_N8-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N16=(abs(Boundary_Flux_Analytical_Whole_N16-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N32=(abs(Boundary_Flux_Analytical_Whole_N32-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N64=(abs(Boundary_Flux_Analytical_Whole_N64-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N128=(abs(Boundary_Flux_Analytical_Whole_N128-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

surf(TwoD_mean_error_fluxsum_N64 , 'FaceColor','g');
hold on;
surf(TwoD_mean_error_flux_N64, 'FaceColor','r');
zlabel('Mean error','FontSize', 20);
title('Comparison of predircting boundary fluxes by two methods (N=64)', 'FontSize',20);
legend('New proposed method','Current method');

diff_mean_error_flux=TwoD_mean_error_flux_N32-TwoD_mean_error_fluxsum_N32;
% surf(diff_mean_error_flux);
% zlabel('Mean error','FontSize', 20);
% title('Comparison of predircting boundary fluxes by two methods (N=64)', 'FontSize',20);

count=0;
for i=1:150
    for j=1:1:150
        if diff_mean_error_flux(i,j)>0
            count=count+1;
        end
    end
end
percent_flux=count/(150*150);

% % pie chart of the mean error;
% Total_mean_error_fluxsum_N8=sum(sum(TwoD_mean_error_fluxsum_N8));
% Total_mean_error_flux_N8=sum(sum(TwoD_mean_error_flux_N8));
% Relative_mean_error_fluxsum_N8=Total_mean_error_fluxsum_N8/(Total_mean_error_fluxsum_N8+Total_mean_error_flux_N8);
% Relative_mean_error_flux_N8=Total_mean_error_flux_N8/(Total_mean_error_fluxsum_N8+Total_mean_error_flux_N8);
% W=[Relative_mean_error_flux_N8 Relative_mean_error_fluxsum_N8];
% pie(W);
% 
% Total_mean_error_fluxsum_N16=sum(sum(TwoD_mean_error_fluxsum_N16));
% Total_mean_error_flux_N16=sum(sum(TwoD_mean_error_flux_N16));
% Relative_mean_error_fluxsum_N16=Total_mean_error_fluxsum_N16/(Total_mean_error_fluxsum_N16+Total_mean_error_flux_N16);
% Relative_mean_error_flux_N16=Total_mean_error_flux_N16/(Total_mean_error_fluxsum_N16+Total_mean_error_flux_N16);
% W=[Relative_mean_error_flux_N16 Relative_mean_error_fluxsum_N16];
% pie(W);
% 
% Total_mean_error_fluxsum_N32=sum(sum(TwoD_mean_error_fluxsum_N32));
% Total_mean_error_flux_N32=sum(sum(TwoD_mean_error_flux_N32));
% Relative_mean_error_fluxsum_N32=Total_mean_error_fluxsum_N32/(Total_mean_error_fluxsum_N32+Total_mean_error_flux_N32);
% Relative_mean_error_flux_N32=Total_mean_error_flux_N32/(Total_mean_error_fluxsum_N32+Total_mean_error_flux_N32);
% W=[Relative_mean_error_flux_N32 Relative_mean_error_fluxsum_N32];
% pie(W);
% 
% Total_mean_error_fluxsum_N64=sum(sum(TwoD_mean_error_fluxsum_N64));
% Total_mean_error_flux_N64=sum(sum(TwoD_mean_error_flux_N64));
% Relative_mean_error_fluxsum_N64=Total_mean_error_fluxsum_N64/(Total_mean_error_fluxsum_N64+Total_mean_error_flux_N64);
% Relative_mean_error_flux_N64=Total_mean_error_flux_N64/(Total_mean_error_fluxsum_N64+Total_mean_error_flux_N64);
% W=[Relative_mean_error_flux_N64 Relative_mean_error_fluxsum_N64];
% pie(W);
%% Identify the significant frequencies of the synthetic landscape;
% pad = 0; % 1 means pad the data, 0 no padding.
% window = 0; % 1 means window the data, 0 no window
% dx=1;
% dy=1;
% %% 2D FFT of the original surface and 2D power spectrum plot;
% [P, Pm, fm, Pv, fv] = fft2D(Z_Synthetic,dx,dy,pad,window); % Apply the 2D fft;
% figure('Name','Fig. 8: Original 2D Spectrum','NumberTitle','off')
% SpecPlot2D(fm,Pm); % Plot the 2D normalized power spectrum probability;
% xt=get(gca,'XTick');
% set(gca,'FontSize',9);
% xt=get(gca,'YTick');
% set(gca,'FontSize',9)
% title('2D power spectrum of the whole part of the surface','FontSize',20);
% 
% 



