%% Wavelet analysis of a bed-form profile retrieved from the synthetic landscape %%
y=101;
BFP=Retrieve_Bedform(y,Topo_whole);
dlmwrite('BFP.txt', BFP, 'delimiter', '\t'); % export the 2D coordinate measurements;
% BedformsATM_App1_WaveletAnalysis

BFP=BFP(:,3);

mother='morl';
t=1:150;
w0=6;
dt=1;
ds=0.4875;
s0=2*dt;
nb=fix(log2(length(BFP'))/ds)+1;
scales=s0*2.^((0:nb-1)*ds);

% Calculate the temporal-frequency domain coeffs;
coefs=cwt(BFP,scales, mother);

% Transform the scales to the frequency;
f=scal2frq(scales, mother, dt);
w=1./f;

figure(10)
% x_tot=[x1,x2];
subplot('position',[0.1 0.75 0.65 0.2])
x1=1:150;
y1=1:150;
surfl(x1,y1,topopiece_vector.z);colormap copper
shading interp
% set(gca,'XLim',xlim(:))
xlabel('x (m)','FontSize',20)
ylabel('y (m)','FontSize',20)
zlabel('Elevation (m)','FontSize',20)
hold on

x=1:100;
y=150*ones(1,100);
plot3(x,y,topopiece_vector.z(150,21:120),'Color','blue','LineWidth',5);
% plot3(x1,101*ones(1,100),Topo_whole(100,:),'LineWidth',3,'Color','r');
% text(x(ind),y(ind),'Z = 0','fontweight','bold','color','y')
% plot(Topo_whole(100,:),'magenta','linewidth',4)
title('Whole Topography','FontSize',25)

subplot('position',[0.1 0.4 0.61 0.2])
plot(t,BFP)
% set(gca,'XLim',xlim(:))
xlabel('x (m)','FontSize',20)
ylabel('Elevation (m)','FontSize',20)
title('Real Bedform Profile2','FontSize',25)

subplot('position',[0.1 0.1 0.65 0.2])
% Plot the 2D temporal-frequency power spectrum;
contour(t,w,abs(coefs),'lineStyle','none','LineColor',[0 0 0],'Fill','on')
xlabel('x (m)','FontSize',20)
ylabel('Wavelength','FontSize',20)
title('Wavelet Power Spectrum','FontSize',25)
set(gcf,'Colormap',jet)
set(gca, 'Ylim', [min(w) max(w)],'XGrid','On','YGrid','On')
colorbar
hold on
xval=150;
x=[xval,xval];
y=[2.5,40];
plot(x,y,'LineWidth',1,'Color','r');


