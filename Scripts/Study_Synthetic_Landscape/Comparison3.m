%% Comparison 3: Comparison between original surface and the Fourier Fitting Surface with different number of frequencies %%

% This comparison answers the question: How good is the Fourier series
% fitting in reconstructing the topography? Does it depend on the number of
% the frequencies in each direction?

% This MATLAB code conducts the comparison 3 analysis: the original surface vs. 
% Fourier Fitting surface with different number of frequencies in each
% direction;
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
%% Solve spectral problem
mirror=1; 
periodic=0;
% mirror=1;periodic=1;
N=16; 
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
%% Present the comparison
% plot the approximate surface
Surface_FS=h;
surf(X,Y,Surface_FS);
% surfl(X,Y,h), colormap copper
shading interp

figure(10)
subplot(1,2,1)
contourf(X,Y,Z_Synthetic);
colorbar
title('Original Surface','FontSize',40)
ax = gca;
ax.FontSize = 24; 

subplot(1,2,2)
contourf(X,Y,Surface_FS);
colorbar
title('Fourier Fitting Surface','FontSize',40)
ax = gca;
ax.FontSize = 24; 
%% Evaluate the coefficient paramhat
% for i=1:N
%     paramhat_matrix(i,1:N)=paramhat((i-1)*N+1:N*i);   
% end


Kx = reshape(kx,N,N);
Ky = reshape(ky,N,N);
paramhat_matrix = reshape(paramhat,N,N);


contourf(Kx,Ky,paramhat_matrix)
title('Frequencies of Fourier Series Topography','FontSize',40)
xlabel('Frequencies in x direction','FontSize',30);
ylabel('Frequencies in y direction','FontSize',30);
ax = gca;
ax.FontSize = 24;
colorbar
hold on
plot([kx1_Synthetic,kx2_Synthetic],[ky1_Synthetic,ky2_Synthetic],'or','MarkerSize',15);
hold on

abs_paramhat_matrix=abs(paramhat_matrix);
Max_paramhat=max(max(abs_paramhat_matrix));
Imp_Fre_1=abs_paramhat_matrix==Max_paramhat;
kx_Imp1=Kx(Imp_Fre_1);
ky_Imp1=Ky(Imp_Fre_1);

abs_paramhat_matrix(abs_paramhat_matrix==Max_paramhat)=-inf;
Second_Max_paramhat=max(max(abs_paramhat_matrix));
Imp_Fre_2=abs_paramhat_matrix==Second_Max_paramhat;
kx_Imp2=Kx(Imp_Fre_2);
ky_Imp2=Ky(Imp_Fre_2);

plot([kx_Imp1,kx_Imp2],[ky_Imp1,ky_Imp2],'blue','MarkerSize',15);
