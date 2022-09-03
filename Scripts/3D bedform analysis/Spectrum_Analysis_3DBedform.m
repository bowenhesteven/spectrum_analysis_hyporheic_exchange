%% Spectrum analysis on the simulated synthetic 3D bedforms %%
% Implement the spectrum analysis on the simulated synthetic 3D bedforms
% and identify the significant frequency signals of the bedform
% Compare the real signals with the assigned frequencies of the signal of
% the Fourier fitted surface (N=32 fitting frequencies)
%% Identify the significant frequencies of the synthetic landscape;
pad = 0; % 1 means pad the data, 0 no padding.
window = 0; % 1 means window the data, 0 no window
dx=1;
dy=1;
%% 2D FFT of the original surface and 2D power spectrum plot;
h=Z_H;
[P, Pm, fm, Pv, fv] = fft2D(h,dx,dy,pad,window); % Apply the 2D fft;
figure('Name','Fig. 8: Original 2D Spectrum','NumberTitle','off')
SpecPlot2D(fm,Pm); % Plot the 2D normalized power spectrum probability;
xt=get(gca,'XTick');
set(gca,'FontSize',25);
xt=get(gca,'YTick');
set(gca,'FontSize',25)
title('2D power spectrum of the 3D Simulated Synthetic Surface','FontSize',40);
%% Fourier fitting surface
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
[nr,nc]=size(Z_H); % Get the x coordinate and y coordinate size of the cropped area of image;
X = 1:nc;
Y = 1:nr; 
% X=X+100;
% Y=Y+100;
%[paramhat,lambda,kx,ky,h,mean_error,av]=Spectop_old(X,Y,Z,N,f1,f2,f3,mirror,ratio_low,ratio_high,cross);% Use with spectop_old
[paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(X,Y,Z_H,N,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
figure(10)
subplot(2,1,1)
Surface_FS=h;
surf(X,Y,Surface_FS);
% surfl(X,Y,h), colormap copper
shading interp
title('Fourier Fitted Surface (N=32)','FontSize',40)

subplot(2,1,2)
surfl(x,y,Z_H);colormap copper
shading interp
title('Real simulated 3D Surface ','FontSize',40)
%% Plot the assigned lamda/frequency used in each Fourier fitted surface
frequency=1./lambda;
Kx=frequency;
Ky=frequency;
[Kx,Ky]=meshgrid(Kx,Ky);
KxKy_1=[Kx(:),Ky(:)];
KxKy_2=[-Kx(:),Ky(:)];
KxKy_3=[-Kx(:),-Ky(:)];
KxKy_4=[Kx(:),-Ky(:)];
FrePairs=[KxKy_1;KxKy_2;KxKy_3;KxKy_4];
% Plot the evaluating frequency pairs;
scatter(FrePairs(:,1),FrePairs(:,2));
hold on
sz=200;
scatter(1/40,1/60,sz,'MarkerFaceColor','r');
hold on
scatter(1/40,-1/60,sz,'MarkerFaceColor','r');
hold on
scatter(-1/40,1/60,sz,'MarkerFaceColor','r');
hold on
scatter(-1/40,-1/60,sz,'MarkerFaceColor','r');
xt=get(gca,'XTick');
set(gca,'FontSize',25);
xt=get(gca,'YTick');
set(gca,'FontSize',25)
xlim([-0.5 0.5])
xticks(-0.5:(1/29):0.5);
ylim([-0.5 0.5])
yticks(-0.5:(1/29):0.5);
title('2D assigned frequency pairs on the Fourier fitted surface (N=32)','FontSize',40)
xlabel('x frequency','FontSize',30);
ylabel('y frequency','FontSize',30);

