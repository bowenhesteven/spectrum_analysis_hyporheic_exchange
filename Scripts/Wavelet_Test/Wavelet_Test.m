%% Anaysis 
%% Create two sections of the synthetic topography and combine them to the whole topography;
A0=0;
H_Coef1=30;
H_Coef2=30;
lamda1_x=8;
lamda1_y=10;
lamda2_x=15;
lamda2_y=7;
x1=101:150;
x2=151:200;
y1=101:200;
y2=101:200;
[section1,section2,Topo_whole]=Synthetic_Surface_2Sections(lamda1_x,lamda1_y,lamda2_x,lamda2_y,A0,H_Coef1,H_Coef2,x1,y1,x2,y2);
%% Plot the synthetic topography;
%--- Plot section1 topography
figure(10)
subplot('position',[0.1 0.75 0.65 0.2])
surfl(x1,y1,section1);colormap copper
shading interp
% set(gca,'XLim',xlim(:))
xlabel('x (m)')
ylabel('y (m)')
zlabel('elevation (m)')

title('Section 1 topography')
hold on

%--- Plot section 2 topography
subplot('position',[0.1 0.4 0.65 0.2])
surfl(x2,y2,section2);colormap copper
shading interp
% set(gca,'XLim',xlim(:))
xlabel('x (m)')
ylabel('y (m)')
zlabel('elevation (m)')

title('Section 2 topography')
hold on


%--- Plot the whole topography
subplot('position',[0.1 0.1 0.65 0.2])
x_tot=[x1,x2];
surfl(x_tot,y2,Topo_whole);colormap copper
shading interp
% set(gca,'XLim',xlim(:))
xlabel('x (m)')
ylabel('y (m)')
zlabel('elevation (m)')

title('Whole topography')
hold on
%% Export the synthetic surface (coordinates) to COMSOL to numerically calculate the boundary fluxes
% Export the coordinates of the synthetic surface;
[X1,Y1] = meshgrid(x_tot,y1);
A=[X1(:),Y1(:)]; % create coordinate matrix;
dlmwrite('Coordinate.txt', A, 'delimiter', '\t'); % export the 2D coordinate measurements;

% Export the synthetic surface;
DEM_Synthetic=[X1(:),Y1(:),Topo_whole(:)];
dlmwrite('DEM_Synthetic.txt', DEM_Synthetic, 'delimiter', '\t'); % export the 2D coordinate measurements

% Read the numerical boundary flux from COMSOL
filename = 'Boundary_flux.txt';
M = dlmread(filename,'',9,0); 
Boundary_flux_Comsol=zeros(100,100);
X_Comsol=zeros(100,100);
Y_Comsol=zeros(100,100);
for i=1:100
    Boundary_flux_Comsol(:,i)=M((i-1)*100+1:100*i,3);   
    X_Comsol(:,i)=M((i-1)*100+1:100*i,1);
    Y_Comsol(:,i)=M((i-1)*100+1:100*i,2);
end
%% Fourier fitting whole topography;
mirror=1; 
periodic=0;
% mirror=1;periodic=1; 
% N=40;
N=32;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(Topo_whole); % Get the x coordinate and y coordinate size of the cropped area of image;
X1 = 1:nc;
Y1 = 1:nr; 
X1=X1+100;
Y1=Y1+100;
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X1,Y1,Topo_whole,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

figure(10)
subplot(1,2,1)
 contourf(Topo_whole)
colorbar
 title('Real Synthetic topography','FontSize',40)
 ax = gca;
 ax.FontSize = 18; 
 
 subplot(1,2,2)
 contourf((h));
 colorbar
 title('Fourier fitting topography','FontSize',40)
 ax = gca;
 ax.FontSize = 18;


% surfl(X1,Y1,h), colormap copper
% shading interp
%axis([40 100 0 10 -0.3 0.3])

%%  Fourier fitting on separate section of topography;
mirror=1; 
periodic=0;
% mirror=1;periodic=1; 
% N=40;
N=32;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(section1); % Get the x coordinate and y coordinate size of the cropped area of image;
X1 = 1:nc;
Y1 = 1:nr; 
X1=X1+100;
Y1=Y1+100;
[paramhat1,lambda1,N1,kx1,ky1,h1,mean_error1,av1]=Spectop(X1,Y1,section1,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

[nr,nc]=size(section2); % Get the x coordinate and y coordinate size of the cropped area of image;
X2 = 1:nc;
Y2 = 1:nr; 
X2=X2+150;
Y2=Y2+100;
[paramhat2,lambda2,N2,kx2,ky2,h2,mean_error2,av2]=Spectop(X2,Y2,section2,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
h_whole=[h1,h2];

figure(10)
subplot(1,2,1)
contourf(Topo_whole)
colorbar
 title('Real Synthetic topography','FontSize',40)
 ax = gca;
 ax.FontSize = 18; 
 
 subplot(1,2,2)
 contourf((h_whole));
 colorbar
 title('Fourier fitting topography','FontSize',40)
 ax = gca;
 ax.FontSize = 18;

%% Evaluate the mean_error of different N/number of fitting frequencies using Fourier fitting on whole topograhy
 mean_error_matrix=zeros(1,10);
 % N_WL=[2,4,8,16,32,64];
 N_f=[2,4,8,12,16,20,24,30,38,45];
 for i=1:10
     N=N_f(i);
     [paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X1,Y1,Topo_whole,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
     mean_error_matrix(i)=mean_error;
 end
 plot(N_WL,mean_error_matrix,'-o',...
     'LineWidth',3,...
     'MarkerEdgeColor','b',...
     'MarkerSize',10)
 title('Influence of Number of Frequencies','FontSize',40)
 ax = gca;
 ax.FontSize = 18; 
 xlabel('Number of frequencies in each direction','FontSize',40)
 ylabel('Mean error','FontSize',40)
%% Evaluate the mean_error of different  N/number of fitting frequencies using Fourier fitting on section1 topograhy
% mean_error_matrix=zeros(1,7);
% mean_error_matrix1=zeros(1,6);
% N_WL=[2,4,8,16,32,64];
% for i=1:6
%     N_f=2^i;
%     [paramhat1,lambda1,N1,kx1,ky1,h1,mean_error1,av1]=Spectop(x1,y1,section1,N_f,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
%     % [paramhat2,lambda2,N2,kx2,ky2,h2,mean_error2,av2]=Spectop(x2,y2,section2,N_f,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
%     % h=[h1,h2];
%     % aaa=h-Topo_whole;
%     % mean_error=mean(mean(abs(aaa)))/abs(min(min(Topo_whole))-max(max(Topo_whole)));
%     mean_error_matrix1(i)=mean_error1;
% end
% plot(N_WL,mean_error_matrix1,'-o',...
%     'LineWidth',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerSize',10)
% title('Influence of Number of Frequencies','FontSize',40)
% ax = gca;
% ax.FontSize = 18; 
% xlabel('Number of frequencies in each direction','FontSize',40)
% ylabel('Mean error','FontSize',40)

%% Evaluate the mean_error of different  N/number of fitting frequencies using Fourier fitting on section2 topograhy
% mean_error_matrix=zeros(1,7);
% mean_error_matrix2=zeros(1,6);
% N_WL=[2,4,8,16,32,64];
% for i=1:6
%     N_f=2^i;
%     [paramhat2,lambda2,N2,kx2,ky2,h2,mean_error2,av2]=Spectop(x2,y2,section2,N_f,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
%     % [paramhat2,lambda2,N2,kx2,ky2,h2,mean_error2,av2]=Spectop(x2,y2,section2,N_f,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
%     % h=[h1,h2];
%     % aaa=h-Topo_whole;
%     % mean_error=mean(mean(abs(aaa)))/abs(min(min(Topo_whole))-max(max(Topo_whole)));
%     mean_error_matrix2(i)=mean_error2;
% end
% plot(N_WL,mean_error_matrix2,'-o',...
%     'LineWidth',3,...
%     'MarkerEdgeColor','b',...
%     'MarkerSize',10)
% title('Influence of Number of Frequencies','FontSize',40)
% ax = gca;
% ax.FontSize = 18; 
% xlabel('Number of frequencies in each direction','FontSize',40)
% ylabel('Mean error','FontSize',40)

%% Evaluate the mean_error of different N/number of fitting frequencies using Fourier fitting on separate topograhy
mean_error_matrix=zeros(1,10);
% N_WL=[2,4,8,16,32,64];
 N_f=[2,4,8,12,16,20,24,30,38,45];
for i=1:10
    N=N_f(i);
    [paramhat1,lambda1,N1,kx1,ky1,h1,mean_error1,av1]=Spectop(x1,y1,section1,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
    [paramhat2,lambda2,N2,kx2,ky2,h2,mean_error2,av2]=Spectop(x2,y2,section2,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
    h=[h1,h2];
    aaa=h-Topo_whole;
    mean_error=mean(mean(abs(aaa)))/abs(min(min(Topo_whole))-max(max(Topo_whole)));
    mean_error_matrix(i)=mean_error;
end
plot(N_f,mean_error_matrix,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','b',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',40)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of frequencies in each direction','FontSize',40)
ylabel('Mean error','FontSize',40)

%% Plot the Influence of the number of Frequencies on Fourier fitting of the whole topography and separate sections
% N_WL=[2,4,8,16,32,64];
N_f=[2,4,8,12,16,20,24,30,38,45];
plot(N_f,mean_error_matrix,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','b',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',60)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of frequencies in each direction','FontSize',20)
ylabel('Mean error','FontSize',30)
 hold on
 plot(N_f,mean_error_matrix_separate,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','r',...
    'MarkerSize',10)
legend('Fourier fitting of whole topography','Fourier fitting of separate sections of topography')
%% Compare the fluxes in two scenarios: (1) Fourier fitting on the whole topography (2) Fourier fitting ont the separate topography;
% Fourier fitting on the whole topography
Cond=.00001; % Cond = Hydraulic conductivity.
porosity=0.3; % Porosity = Effective porosity.
rho=998; % THe density of water is 998 kg/m^3;
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=100; 
% N=40;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(Topo_whole); % Get the x coordinate and y coordinate size of the cropped area of image;
X1 = 1:nc;
Y1 = 1:nr; 
X1=X1+100;
Y1=Y1+100;
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X1,Y1,Topo_whole,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

dp=10;
W_top_FS=zeros(100,100);
for x=101:200
    for y=101:200        
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat',Cond,porosity); % solve the fluxes point by point;
W_top_FS(101-(y-100),x-100)=W;
    end
end
Boundary_flux_FS=W_top_FS*porosity;
Boundary_flux_CFD=flip(Boundary_flux_Comsol)/rho;
mean_error_flux=mean(mean(abs(Boundary_flux_FS-Boundary_flux_CFD)))/abs(min(min(Boundary_flux_FS))-max(max(Boundary_flux_CFD)));
%% (2) Fourier fitting on the separate topography;
Cond=.00001; % Cond = Hydraulic conductivity.
porosity=0.3; % Porosity = Effective porosity.
rho=998; % THe density of water is 998 kg/m^3;
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=100; 
% N=40;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
[nr,nc]=size(section1); % Get the x coordinate and y coordinate size of the cropped area of image;
X1 = 1:nc;
Y1 = 1:nr; 
X1=X1+100;
Y1=Y1+100;
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat1,lambda1,N1,kx1,ky1,h1,mean_error1,av1]=Spectop(X1,Y1,section1,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

dp=10;
W_top_FS1=zeros(100,50);
for x=101:150
    for y=101:200        
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx1,ky1,paramhat1',Cond,porosity); % solve the fluxes point by point;
W_top_FS1(101-(y-100),x-100)=W;
    end
end


[nr,nc]=size(section2); % Get the x coordinate and y coordinate size of the cropped area of image;
X2 = 1:nc;
Y2 = 1:nr; 
X2=X2+150;
Y2=Y2+100;

[paramhat2,lambda2,N2,kx2,ky2,h2,mean_error2,av2]=Spectop(X2,Y2,section2,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);

% dp=10;
W_top_FS2=zeros(100,50);
for x=151:200
    for y=101:200        
z=0;
[U,V,W]=Specvel_New(x,y,z,dp,kx2,ky2,paramhat2',Cond,porosity); % solve the fluxes point by point;
W_top_FS2(101-(y-100),x-150)=W;
    end
end
 W_top_Combine=[W_top_FS1,W_top_FS2];

Boundary_flux_FS=W_top_Combine*porosity;
Boundary_flux_CFD=flip(Boundary_flux_Comsol)/rho;
mean_error_flux=mean(mean(abs(Boundary_flux_FS-Boundary_flux_CFD)))/abs(min(min(Boundary_flux_FS))-max(max(Boundary_flux_CFD)));
%% Plot the Influence of the number of Frequencies on Fourier fitting of the whole topography and separate sections
N_WL=[8,16,32,64,82,100,128];
plot(N_WL,mean_error_flux_FFTopo,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','b',...
    'MarkerSize',10)
title('Influence of Number of Frequencies on approximating boundary fluxes','FontSize',70)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of frequencies in each direction','FontSize',20)
ylabel('Mean error','FontSize',30)
 hold on
 plot(N_WL,mean_error_flux_FFSeparate,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','r',...
    'MarkerSize',10)
legend('Fourier fitting of whole topography','Fourier fitting of separate sections of topography')



figure(10)
subplot(1,2,1)
 contourf(Boundary_flux_FS)
colorbar
 title('Fluxes of Fourier Series Fitting Surface','FontSize',40)
 ax = gca;
 ax.FontSize = 18; 
 
 subplot(1,2,2)
 contourf((Boundary_flux_CFD));
 colorbar
 title('Fluxes of Real Synthetic Surface','FontSize',40)
 ax = gca;
 ax.FontSize = 18;

























