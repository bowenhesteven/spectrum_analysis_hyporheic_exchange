%% Comparison 1: Fluxes calculated by analytical solution vs. Fluxes calcualated by numerical solution in COMSOL %%

% This MATLAB code conducts the comparison 1 analysis: the fluxes
% calculated by analytical solution vs. fluxes calcualated by numerical
% solution in COMSOL
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
%% Calculate the velocity field; two scenarios;
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
%% Comparison of topflux with COMSOL Model/Darcy flux

% Import the top flux field into MATLAB;
filename = 'Top_Flux_Synthetic.txt';
M = dlmread(filename,'',9,0); 

Top_Velocity_Field_Comsol=zeros(150,150);
X_Comsol=zeros(150,150);
Y_Comsol=zeros(150,150);
for i=1:150
    Top_Velocity_Field_Comsol(:,i)=M((i-1)*150+1:150*i,3);   
    X_Comsol(:,i)=M((i-1)*150+1:150*i,1);
    Y_Comsol(:,i)=M((i-1)*150+1:150*i,2);
end
%% Plot the comparison (qualitatively)
Top_Velocity_Field_True = W_top*porosity; 
x=1001:1150;
y=1001:1150;
[X_True,Y_True] = meshgrid(x,y);
rho=998; % The density of water is 998 kg/m^3;
figure(10)

subplot(1,2,1)
Top_Velocity_Field_Comsol_rho=Top_Velocity_Field_Comsol/rho;
contourf(Top_Velocity_Field_Comsol_rho)
colorbar
title('Fluxes of Numerical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2)
contourf((Top_Velocity_Field_True));
colorbar
title('Fluxes of Analytical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 
%% Plot the comparison (quantitatively, use histogram plot)
figure(10)
subplot(1,2,1)
% Top_Velocity_Field_Comsol_rho=Top_Velocity_Field_Comsol/rho;
% contourf(Top_Velocity_Field_Comsol_rho)
% colorbar
histogram(Top_Velocity_Field_Comsol_rho);
title('Fluxes of Numerical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2)
% contourf((Top_Velocity_Field_True));
% colorbar
histogram(Top_Velocity_Field_True);
title('Fluxes of Analytical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 
%% Calculate the total mass flux for different scenarios;
% Calculate the total mass flux out of the top of the surface/analytic solution;
Index_out=Top_Velocity_Field_True>0; % retrieve the index of the grid/position that the flux out of the surface;
Top_Velocity_Field_True_Out=Top_Velocity_Field_True.*Index_out; % create the matrix presenting the grid that flux out of the surface;
% W_Surf_out_vector=W_Surf(Index_out); % create the vector presenting all the grid that flux out of the surface;
Mass_Flux_out_matrix=rho*Top_Velocity_Field_True_Out; % Calculate the total mass surface flux out of the surface (matrix);
Total_Mass_Flux_Out=sum(sum(Mass_Flux_out_matrix)); % sum up the total mass flux on the top of the surface;

% Calculate the total mass flux into the top of the surface/analytic
% solution;
Index_in=Top_Velocity_Field_True<0; % retrieve the index of the grid/position that the flux out of the surface;
Top_Velocity_Field_True_In=Top_Velocity_Field_True.*Index_in; % create the matrix presenting the grid that flux out of the surface;
% W_Surf_out_vector=W_Surf(Index_out); % create the vector presenting all the grid that flux out of the surface;
Mass_Flux_In_matrix=rho*Top_Velocity_Field_True_In; % Calculate the total mass surface flux out of the surface (matrix);
Total_Mass_Flux_In=sum(sum(Mass_Flux_In_matrix)); % sum up of the total mass flux on the top surface;


