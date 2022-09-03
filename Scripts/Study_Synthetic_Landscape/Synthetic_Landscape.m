%% Create a synthetic landscape%%

% Function:
% Create a synthetic landscape, with two frequencies in each direction (x,
% y direction)


% Input:
% A0= average surface elevation, m
% lambda1= the first wavelength/lambda in two direction,
% lambda2= the second wavelength/lambda in two direction,
% H_Coef1= the coefficient of the first frequency pair,
% H_Coef2= the coefficient of the second frequency pair,
% H_Coef3= the coefficient of the third frequency pair,
% H_Coef4= the coefficient of the fourth frequency pair,
% x= coordinate in x-direction,
% y= coordinate in y-direction,

% Output:
% Z_Synthetic= the z elevation of the synthetic surface
% H_Coef= the generation of coefficient vector
% lambda_Synthetic= the generation of lambda/wavelength vector
% kx_Synthetic= the generation of x frequencies vector
% ky_Synthetic= the generation of y frequencies vector

function [Z_Synthetic, H_Coef, lambda_Synthetic, kx_Synthetic, ky_Synthetic]=Synthetic_Landscape(A0,lambda1,lambda2,H_Coef1,H_Coef2,H_Coef3,H_Coef4,x,y)

% A0=0; % Average surface elevation;
lambda_Synthetic=zeros(2,1); %
lambda_Synthetic(1,1)=lambda1;
lambda_Synthetic(2,1)=lambda2;
kx1_Synthetic=2*pi/lambda_Synthetic(1,1);
kx2_Synthetic=2*pi/lambda_Synthetic(2,1);
kx_Synthetic=zeros(1,4);
kx_Synthetic(1,1:2)=kx1_Synthetic;
kx_Synthetic(1,3:4)=kx2_Synthetic;
ky1_Synthetic=2*pi/lambda_Synthetic(1,1);
ky2_Synthetic=2*pi/lambda_Synthetic(2,1);
ky_Synthetic=zeros(1,4);
ky_Synthetic(1,1)=ky1_Synthetic;
ky_Synthetic(1,2)=ky2_Synthetic;
ky_Synthetic(1,3)=ky1_Synthetic;
ky_Synthetic(1,4)=ky2_Synthetic;

H_Coef=zeros(4,1);
H_Coef(1,1)=H_Coef1;
H_Coef(2,1)=H_Coef2;
H_Coef(3,1)=H_Coef3;
H_Coef(4,1)=H_Coef4;
X=x;
Y=y;Y = fliplr(Y);
[X,Y]=meshgrid(X,Y);
% Produce the synthetic landscape using the Worman's analytic solution
Z_Synthetic=A0+H_Coef(1,1)*sin(kx_Synthetic(1)*X).*sin(ky_Synthetic(1)*Y)+H_Coef(2,1)*sin(kx_Synthetic(2)*X).*sin(ky_Synthetic(2)*Y)+H_Coef(3,1)*sin(kx_Synthetic(3)*X).*sin(ky_Synthetic(3)*Y)+H_Coef(4,1)*sin(kx_Synthetic(4)*X).*sin(ky_Synthetic(4)*Y);

end