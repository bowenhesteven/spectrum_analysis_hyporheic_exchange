%% Comparison between the analytical solution and the numerical solution: 
%(1) Boxplot %
% -------------------------------------------------------------------------
%% Boxplot of the difference of the two results;
Differ_Fluxes=((Boundary_Flux_Analytical-Boundary_Flux_Numerical))./abs(min(min(Boundary_Flux_Numerical))-max(max(Boundary_Flux_Numerical)));
boxplot(Differ_Fluxes);
xlabel('x (m)','FontSize',30)
ylabel('Mean error (%)', 'FontSize',30)
title('Boxplot of the mean error of the boundary fluxes calculated by two models','FontSize',35)

% Boxplot of the fluxes at all position;
Fluxes_allposition=[Boundary_Flux_Numerical_vector(:),Boundary_Flux_Analytical_vector(:)];
h=boxplot(Fluxes_allposition);
xlabel('All position','FontSize',30)
ylabel('Fluxes (m/s)', 'FontSize',30)
title('Boxplot of the boundary fluxes at all position calculated by two models','FontSize',35)
legend(findobj(gca,'Tag','Box'),'1-Numerical model','2-Analytical model');
% set(findobj(gca,'Type','text'),'FontSize',40);
