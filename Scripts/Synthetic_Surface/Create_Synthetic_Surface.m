% Create a simple harmonic synthetic surface with known frequency (1
% frequency in x-,y-direction, respectively)
%--------------------------------------------------------------------------
A0=0; % Average surface elevation;
lambda1=10.00029493;
lambda2=30;
H_Coef1=30;
H_Coef2=0;
H_Coef3=0;
H_Coef4=0;
x=1+1000:150+1000;
y=1+1000:150+1000;
[Z_Synthetic, H_Coef, lambda_Synthetic, kx_Synthetic, ky_Synthetic]=Synthetic_Landscape(A0,lambda1,lambda2,H_Coef1,H_Coef2,H_Coef3,H_Coef4,x,y);
% surf(x,y,Z_Synthetic);
surfl(x,y,Z_Synthetic);colormap copper
shading interp
xlabel('x-coordinate (m)','FontSize',20)
ylabel('y-coordinate (m)','FontSize',20)
zlabel('Elevation (m)','FontSize',20)
title('Synthetic Surface','FontSize',30)