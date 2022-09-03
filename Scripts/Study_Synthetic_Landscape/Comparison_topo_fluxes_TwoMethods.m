%%%%%%      Comparison    %%%%%% 
% Boundary fluxes calculated by direcly applying the Fourier fitting method
% on the whole topography using Worman's analytical solution
% VS.
% Boundary fluxes calculated by summing up the boundary fluces induced by
% two components: (1) characteristic spatial scales, which is exact
% equivalent as the numerical result of the corresponding part; (2) noise
% spatial scales, which is calculated by applying Fourier fitting method
% Standard: the boundary fluxes numerically calcualted in COMSOL ny
% applying the real natural surface as the hydraulic head boundary
% condition
% This test is based on a real natural surface
%% Plot the two part, respectively
x=1+100:150+100;
y=1+100:150+100;
y=flip(y);
subplot(2,1,1)
surfl(x,y,Topo_Sig);colormap copper
shading interp
title('Characteristic part of the synthetic surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(2,1,2)
surfl(x,y,Topo_Noise);colormap copper
shading interp
title('Noise part of the synthetic surface','FontSize',40)
ax = gca;
ax.FontSize = 18; 
%% Compare the two methods on approximating the surface
Topo_Whole=Topo_Noise+Topo_Sig;

% Evaluate FFS on whole part of the surface;
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=64; 
% N=40;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(Topo_Whole); % Get the x coordinate and y coordinate size of the cropped area of image;
X = 1:nc;
Y = 1:nr; 
X=X+100;
Y=Y+100;
Y=flip(Y);
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Topo_Whole,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

Cond=.00001; % Cond = Hydraulic conductivity.
porosity=0.3; %
dp=10;
W_top_FS=zeros(150,150);
for x=101:250
    for y=101:250      
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat',Cond,porosity); % solve the fluxes point by point;
W_top_FS(151-(y-100),x-100)=W;
    end
end
Boundary_Flux_Analytical_Whole_FFS=W_top_FS*porosity;
mean_error_whole_flux=mean(mean(abs(Boundary_Flux_Analytical_Whole_FFS-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

% Evaluate FFS on the noise part of the surface;
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=64; 
% N=40;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(Topo_Noise); % Get the x coordinate and y coordinate size of the cropped area of image;
X = 1:nc;
Y = 1:nr; 
X=X+100;
Y=Y+100;
Y=flip(Y);
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Topo_Noise,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic); 

Cond=.00001; % Cond = Hydraulic conductivity.
porosity=0.3; %
dp=10;
W_top_FS=zeros(150,150);
for x=101:250
    for y=101:250      
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat',Cond,porosity); % solve the fluxes point by point;
W_top_FS(151-(y-100),x-100)=W;
    end
end
Boundary_Flux_Analytical_Noise_FFS=W_top_FS*porosity;
mean_error_noise_flux=mean(mean(abs(Boundary_Flux_Analytical_Noise_FFS-Boundary_Flux_Numerical_Noise)))/abs(min(min(Boundary_Flux_Numerical_Noise))-max(max(Boundary_Flux_Numerical_Noise)));
%% Calculate the mean error topography for each scenario

% The mean error of the topo using the sum of two parts
Topo_Sum_N8=Topo_Sig+Topo_Noise_N8;
Topo_Sum_N16=Topo_Sig+Topo_Noise_N16;
Topo_Sum_N32=Topo_Sig+Topo_Noise_N32;
Topo_Sum_N64=Topo_Sig+Topo_Noise_N64;
Topo_Sum_N128=Topo_Sig+Topo_Noise_N128;

% Calculate the mean error topography of each scenario
mean_error_toposum_N8=mean(mean(abs(Topo_Sum_N8-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_toposum_N16=mean(mean(abs(Topo_Sum_N16-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_toposum_N32=mean(mean(abs(Topo_Sum_N32-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_toposum_N64=mean(mean(abs(Topo_Sum_N64-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_toposum_N128=mean(mean(abs(Topo_Sum_N128-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));

mean_error_toposum=[mean_error_toposum_N8 mean_error_toposum_N16 mean_error_toposum_N32 mean_error_toposum_N64];

% The mean error of the noise topo
mean_error_toponoise_N8=mean(mean(abs(Topo_Noise_N8-Topo_Noise)))/abs(min(min(Topo_Noise))-max(max(Topo_Noise)));
mean_error_toponoise_N16=mean(mean(abs(Topo_Noise_N16-Topo_Noise)))/abs(min(min(Topo_Noise))-max(max(Topo_Noise)));
mean_error_toponoise_N32=mean(mean(abs(Topo_Noise_N32-Topo_Noise)))/abs(min(min(Topo_Noise))-max(max(Topo_Noise)));
mean_error_toponoise_N64=mean(mean(abs(Topo_Noise_N64-Topo_Noise)))/abs(min(min(Topo_Noise))-max(max(Topo_Noise)));
mean_error_toponoise_N128=mean(mean(abs(Topo_Noise_N128-Topo_Noise)))/abs(min(min(Topo_Noise))-max(max(Topo_Noise)));

mean_error_topo_noise=[mean_error_toponoise_N8 mean_error_toponoise_N16 mean_error_toponoise_N32 mean_error_toponoise_N64];

% The mean error of the whole topo
mean_error_topowhole_N8=mean(mean(abs(Topo_Whole_N8-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_topowhole_N16=mean(mean(abs(Topo_Whole_N16-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_topowhole_N32=mean(mean(abs(Topo_Whole_N32-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_topowhole_N64=mean(mean(abs(Topo_Whole_N64-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
mean_error_topowhole_N128=mean(mean(abs(Topo_Whole_N128-Topo_Whole)))/abs(min(min(Topo_Whole))-max(max(Topo_Whole)));

mean_error_topo_whole=[mean_error_topowhole_N8 mean_error_topowhole_N16 mean_error_topowhole_N32 mean_error_topowhole_N64];
%% Plot the difference in both 1D and 2D surface plot evaluation
% Evaluate the surface topography: noise part, whole scenario, sum scenario

% 1D evaluation
N_WL=[8 16 32 64];
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
%% Combine the characteristic part of the fluxes and the noise part of the fluxes
Boundary_Flux_Analytical_Characteristic=Boundary_Flux_Numerical_Characteristic;
Boundary_Flux_Analytical_Sum_N8=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N8;
Boundary_Flux_Analytical_Sum_N16=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N16;
Boundary_Flux_Analytical_Sum_N32=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N32;
Boundary_Flux_Analytical_Sum_N64=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N64;
% Boundary_Flux_Analytical_Sum_N128=Boundary_Flux_Analytical_Characteristic+Boundary_Flux_Analytical_Noise_N128;
%% Calculate the mean error flux for each scenario

% boundary fluxes sum of the two parts;
mean_error_flux_sum_N8=mean(mean(abs(Boundary_Flux_Analytical_Sum_N8-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_flux_sum_N16=mean(mean(abs(Boundary_Flux_Analytical_Sum_N16-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_flux_sum_N32=mean(mean(abs(Boundary_Flux_Analytical_Sum_N32-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_flux_sum_N64=mean(mean(abs(Boundary_Flux_Analytical_Sum_N64-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
% mean_error_sum_N128=mean(mean(abs(Boundary_Flux_Analytical_Sum_N128-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

mean_error_flux_sum=[mean_error_flux_sum_N8 mean_error_flux_sum_N16 mean_error_flux_sum_N32 mean_error_flux_sum_N64];

% boundary fluxes of the noise parts;
mean_error_flux_noise_N8=mean(mean(abs(Boundary_Flux_Analytical_Noise_N8-Boundary_Flux_Numerical_Noise)))/abs(min(min(Boundary_Flux_Numerical_Noise))-max(max(Boundary_Flux_Numerical_Noise)));
mean_error_flux_noise_N16=mean(mean(abs(Boundary_Flux_Analytical_Noise_N16-Boundary_Flux_Numerical_Noise)))/abs(min(min(Boundary_Flux_Numerical_Noise))-max(max(Boundary_Flux_Numerical_Noise)));
mean_error_flux_noise_N32=mean(mean(abs(Boundary_Flux_Analytical_Noise_N32-Boundary_Flux_Numerical_Noise)))/abs(min(min(Boundary_Flux_Numerical_Noise))-max(max(Boundary_Flux_Numerical_Noise)));
mean_error_flux_noise_N64=mean(mean(abs(Boundary_Flux_Analytical_Noise_N64-Boundary_Flux_Numerical_Noise)))/abs(min(min(Boundary_Flux_Numerical_Noise))-max(max(Boundary_Flux_Numerical_Noise)));
% mean_error_sum_N128=mean(mean(abs(Boundary_Flux_Analytical_Sum_N128-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

mean_error_flux_noise=[mean_error_flux_noise_N8 mean_error_flux_noise_N16 mean_error_flux_noise_N32 mean_error_flux_noise_N64];

% boundary fluxes of the whole parts;
mean_error_flux_whole_N8=mean(mean(abs(Boundary_Flux_Analytical_Whole_N8-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_flux_whole_N16=mean(mean(abs(Boundary_Flux_Analytical_Whole_N16-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_flux_whole_N32=mean(mean(abs(Boundary_Flux_Analytical_Whole_N32-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
mean_error_flux_whole_N64=mean(mean(abs(Boundary_Flux_Analytical_Whole_N64-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
% mean_error_sum_N128=mean(mean(abs(Boundary_Flux_Analytical_Sum_N128-Boundary_Flux_Numerical_Whole)))/abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

mean_error_flux_whole=[mean_error_flux_whole_N8 mean_error_flux_whole_N16 mean_error_flux_whole_N32 mean_error_flux_whole_N64];
%% Plot the difference in both 1D and 2D surface plot evaluation
% Evaluate the boundary fluxes: noise part, whole scenario, sum scenario

% 1D evaluation
N_WL=[8 16 32 64];
plot(N_WL,mean_error_flux_noise,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','b',...
    'MarkerSize',10,'Color','b')
title('Influence of number of fitting frequencies','FontSize',40)
ax = gca;
ax.FontSize = 18; 
hold on

plot(N_WL,mean_error_flux_whole,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','r',...
    'MarkerSize',10,'Color','r')
title('Influence of number of fitting frequencies','FontSize',40)
ax = gca;
ax.FontSize = 18; 
hold on

plot(N_WL,mean_error_flux_sum,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','g',...
    'MarkerSize',10,'Color','g')
title('Influence of number of fitting frequencies','FontSize',50)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of fitting frequencies','FontSize',20)
ylabel('Mean error','FontSize',20)
legend('Mean error of the boundary fluxes of Fourier fitted noise part of the surface','Mean error of the boundary fluxes of the Fourier fitted whole part of the surface','Mean error of the boundary fluxes of the sum of the characteristic part and the Fourier fitted noise part of the surface')


% 2D evaluation
% Calculate the spatial distribution of mean error topography of each scenario
TwoD_mean_error_toposum_N8=(abs(Topo_Sum_N8-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
TwoD_mean_error_toposum_N16=(abs(Topo_Sum_N16-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
TwoD_mean_error_toposum_N32=(abs(Topo_Sum_N32-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
TwoD_mean_error_toposum_N64=(abs(Topo_Sum_N64-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
% TwoD_mean_error_toposum_N128=(abs(Z_Sum_N128-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));

TwoD_mean_error_topo_N8=(abs(Topo_Whole_N8-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
TwoD_mean_error_topo_N16=(abs(Topo_Whole_N16-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
TwoD_mean_error_topo_N32=(abs(Topo_Whole_N32-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
TwoD_mean_error_topo_N64=(abs(Topo_Whole_N64-Topo_Whole))./abs(min(min(Topo_Whole))-max(max(Topo_Whole)));
% TwoD_mean_error_topo_N128=(abs(Z_whole_N128-Z_Synthetic))./abs(min(min(Z_Synthetic))-max(max(Z_Synthetic)));


surf(TwoD_mean_error_toposum_N32 , 'FaceColor','g');
hold on;
surf(TwoD_mean_error_topo_N32, 'FaceColor','r');
zlabel('Mean error','FontSize', 20);
title('Comparison of approximated topography by two methods (N=32)', 'FontSize',20);
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

% Calculate the spatial distribution of mean error boundary fluxes of each scenario
TwoD_mean_error_fluxsum_N8=(abs(Boundary_Flux_Analytical_Sum_N8-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N16=(abs(Boundary_Flux_Analytical_Sum_N16-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N32=(abs(Boundary_Flux_Analytical_Sum_N32-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_fluxsum_N64=(abs(Boundary_Flux_Analytical_Sum_N64-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
% TwoD_mean_error_fluxsum_N128=(abs(Boundary_Flux_Analytical_Sum_N128-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

TwoD_mean_error_flux_N8=(abs(Boundary_Flux_Analytical_Whole_N8-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N16=(abs(Boundary_Flux_Analytical_Whole_N16-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N32=(abs(Boundary_Flux_Analytical_Whole_N32-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
TwoD_mean_error_flux_N64=(abs(Boundary_Flux_Analytical_Whole_N64-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));
% TwoD_mean_error_flux_N128=(abs(Boundary_Flux_Analytical_Whole_N128-Boundary_Flux_Numerical_Whole))./abs(min(min(Boundary_Flux_Numerical_Whole))-max(max(Boundary_Flux_Numerical_Whole)));

surf(TwoD_mean_error_fluxsum_N64 , 'FaceColor','g');
hold on;
surf(TwoD_mean_error_flux_N64, 'FaceColor','r');
zlabel('Mean error','FontSize', 20);
title('Comparison of predircting boundary fluxes by two methods (N=64)', 'FontSize',20);
legend('New proposed method','Current method');



% Standard=Boundary_Flux_Numerical_Whole;
% Analytical=Boundary_Flux_Analytical_Whole;
% Sum=Boundary_Flux_Analytical_Noise+Boundary_Flux_Numerical_Sig;
% mean_error1=abs((Analytical-Standard))./abs(min(min(Standard))-max(max(Standard)));
% surf(mean_error1 , 'FaceColor','g');
% hold on;
% mean_error2=abs((Sum-Standard))./abs(min(min(Standard))-max(max(Standard)));
% hold on;
% surf(mean_error2, 'FaceColor','r');
% zlabel('Mean error','FontSize', 20)
% title('Comparison of fluxes by two methods', 'FontSize',20)