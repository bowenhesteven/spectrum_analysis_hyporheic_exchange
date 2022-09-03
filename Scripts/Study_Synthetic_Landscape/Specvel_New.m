%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%       SPECVEL 2.0      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%* Matlab function, Specvel, that provides an exact (series) solution to 
% groundwater flow in three dimensions based on spectral information of the 
% boundary condition.
% The routine is intended to be executed with input parameters (spectral
% information) from the routine SPECTOP.
% The input is e.g. the original grid (X,Y) and spectral variables from SPECTOP. 
% The routine returns the velocity components and the hydraulic potential in a 
% three-dimensional grid. 
% This function intends to calculate the fluxes point by point, i.e. given
% each 

%* INPUT VARIABLES:
%
% 1. x = x coordinate value
% 2. y = y coordinate value
% 3. z = z coordinate value
% 4. dp = depth of the study domain
% 5. kx = Wave number, 2*pi/lambda, in x-direction with dimension (1:NxN). 
% 6. ky = Wave number in  y-direction with dimension (1:NxN). 
% 7. paramhat = Amplitude coefficients with dimension (1:NxN).
% 8. Cond = Hydraulic conductivity.
% 9. Porosity = Effective porosity.
%
%* RETURN VARIABLES
%
% 1. U = Velocity vector in x-direction of position (x,y,z)
% 
% 2. V = Velocity vector in y-direction of position (x,y,z)

% 3. W = Velocity vector in z-direction of position (x,y,z)

%  Copyright: Anders Worman
%             Environmental Physics Group, SLU
%             email: anders.worman@bt.slu.se
%  10 Jan 2006, MATLAB 7.1 version 
%  IF YOU PUBLISH WORK BENEFITING FROM THIS M-FILE, PLEASE CITE IT AS:
%    Worman, A. (2006) Specvel: A matlab function that provides a Fourier series 
%    solution to groundwater flow, 
%    ftp://ftp.agu.org/apend/" (Username = "anonymous", Password = "guest") 


function [U,V,W]=Specvel_New(x,y,z,dp,kx,ky,paramhat,Cond,porosity)
Cond_ef=Cond/porosity; % conductivity/porosity;
 N=size(paramhat,2); 
U=0;
V=0;
W=0;
for p=1:N % 1:28, N=the number wavelength in x-direction or y-direction;
	Umx=Cond_ef*paramhat(p)*kx(p); % Umx=Cond_ef*paramhat(1)*kx(1)....Umx=Cond_ef*paramhat(28)*kx(28); paramhat=Amplitude coefficients with dimension (1:NxN), a number;
	Umy=Cond_ef*paramhat(p)*ky(p); % Umy=Cond_ef*paramhat(1)*ky(1)....Umy=Cond_ef*paramhat(28)*ky(28); a number;
	Umz=Cond_ef*paramhat(p)*sqrt(kx(p)^2+ky(p)^2); % Umz=Cond_ef*paramhat(1)*sqrt(kx(1)^2+ky(1)^2)....Cond_ef*paramhat(28)*sqrt(kx(28)^2+ky(28)^2); a number;
     nom=exp(z*sqrt(kx(p)^2+ky(p)^2))+exp(sqrt(kx(p)^2+ky(p)^2)*(-2*dp-z)); % exp(z3D*sqrt(kx(1)^2+ky(1)^2))+exp(sqrt(kx(1)^2+ky(1)^2)*(-2*dp-z3D))...exp(z3D*sqrt(kx(28)^2+ky(28)^2))+exp(sqrt(kx(28)^2+ky(28)^2)*(-2*dp-z3D)); [151,151,26]
     denom=1+exp(-sqrt(kx(p)^2+ky(p)^2)*2*dp); % denom=1+exp(-sqrt(kx(1)^2+ky(1)^2)*2*dp)...1+exp(-sqrt(kx(28)^2+ky(28)^2)*2*dp); a number;
     F1=nom/denom; % [151,151,26];
	% h=h+paramhat(p)*F1.*(sin(x3D*kx(p)).*cos(y3D*ky(p))); % h=h+paramhat(1)*F1.*(sin(x3D*kx(1)).*cos(y3D*ky(1)))...h+paramhat(28)*F1.*(sin(x3D*kx(28)).*cos(y3D*ky(28))); [151,151,26];
    U=U-Umx*F1.*(cos(x*kx(p)).*cos(y*ky(p))); % U=U-Umx*F1.*(cos(x3D*kx(1)).*cos(y3D*ky(1)))...U-Umx*F1.*(cos(x3D*kx(28)).*cos(y3D*ky(28))); [151,151,26];
    V=V+Umy*F1.*(sin(x*kx(p)).*sin(y*ky(p))); % V=V+Umy*F1.*(sin(x3D*kx(1)).*sin(y3D*ky(1)))...V+Umy*F1.*(sin(x3D*kx(28)).*sin(y3D*ky(28))); [151,151,26];
     nom=exp(sqrt(kx(p)^2+ky(p)^2)*z)-exp(sqrt(kx(p)^2+ky(p)^2)*(-2*dp-z)); % nom=exp(sqrt(kx(1)^2+ky(1)^2)*z3D)-exp(sqrt(kx(1)^2+ky(1)^2)*(-2*dp-z3D))....exp(sqrt(kx(28)^2+ky(28)^2)*z3D)-exp(sqrt(kx(28)^2+ky(28)^2)*(-2*dp-z3D)); [151,151,26];
     denom=1+exp(-sqrt(kx(p)^2+ky(p)^2)*2*dp); % denom=1+exp(-sqrt(kx(1)^2+ky(1)^2)*2*dp)...1+exp(-sqrt(kx(28)^2+ky(28)^2)*2*dp), a number;
     F2=nom/denom; % [151,151,26];
    W=W-Umz*F2.*(sin(x*kx(p)).*cos(y*ky(p))); % W=W-Umz*F2.*(sin(x3D*kx(1)).*cos(y3D*ky(1)))...W-Umz*F2.*(sin(x3D*kx(28)).*cos(y3D*ky(28))); [151,151,26];
end