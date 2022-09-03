%% Comparison between the analytical solution and the numerical solution: 
%(1) qualitatively %
% -------------------------------------------------------------------------
%% Plot the comparison (qualitatively)
Boundary_Flux_Analytical = W_top*porosity; % Results from Comsol is considered as true;
x=1001:1150;
y=1001:1150;
[X_Ana,Y_Ana] = meshgrid(x,y);
rho=998; % The density of water is 998 kg/m^3;
figure(10)

subplot(1,2,1)
% Boundary_Flux_Numerical=Top_Velocity_Field_Comsol/rho;
contourf(Boundary_Flux_Numerical)
colorbar
title('Numerically Solved Fluxes','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2)
contourf((Boundary_Flux_Analytical));
colorbar
title('Analytically Solved Fluxes','FontSize',40)
ax = gca;
ax.FontSize = 18; 

% Plot the absolute difference
Flux_Difference=Boundary_Flux_Analytical-Boundary_Flux_Numerical;
histogram(Flux_Difference);
Flux_Difference(abs(Flux_Difference)>0.5E-04)=0.5E-04;
contourf(Flux_Difference);
colorbar
title('Analytical Boundary Flux-Numerical Boundary Flux','FontSize',40)
ax = gca;
ax.FontSize = 18; 
