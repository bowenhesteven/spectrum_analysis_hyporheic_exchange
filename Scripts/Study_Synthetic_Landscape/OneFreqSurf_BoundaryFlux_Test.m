%%% Example of Calculating the boundary fluxes using Worman analytical
%%% solution, the synthetic landscape has just one frequency pair %%%
% Instead of using his original code, this version code calculates the
% boundary fluxes in a for loop function %%%

h_m=30;
K=0.00001;
z=0;
dp=10;
Kx=0.6283;
Ky=0.6283;

% Plot the synthetic surface;
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z=30*sin(0.6283*X).*sin(0.6283*Y);
surfl(x,y,z);colormap copper
shading interp
%---------------------------



Nom1=exp(sqrt(Kx^2+Ky^2)*z)+exp(sqrt(Kx^2+Ky^2)*(-2*dp-z));
Nom2=exp(sqrt(Kx^2+Ky^2)*z)-exp(sqrt(Kx^2+Ky^2)*(-2*dp-z));
Denom=1+exp(-2*sqrt(Kx^2+Ky^2)*dp);
KyKx_2=sqrt(Kx^2+Ky^2);
W_test_matrix=zeros(150,150);

for x=1001:1150
    for y=1001:1150
W_test=h_m*sin(Kx*x)*cos(Ky*y)*(KyKx_2/Denom)*Nom2*(-K);
W_test_matrix(151-(y-1000),x-1000)=W_test;
    end
end


% filename = 'Boundary_Flux_OneFre.txt';
% M = dlmread(filename,'',9,0); 
% 
% Top_Velocity_Field_Comsol_OneFre=zeros(150,150);
% X_Comsol=zeros(150,150);
% Y_Comsol=zeros(150,150);
% for i=1:150
%     Top_Velocity_Field_Comsol_OneFre(:,i)=flip(M((i-1)*150+1:150*i,3));   
%     X_Comsol(:,i)=M((i-1)*150+1:150*i,1);
%     Y_Comsol(:,i)=flip(M((i-1)*150+1:150*i,2));
% end
% Import the top flux field into MATLAB;
filename = 'Boundary_Flux.txt';
M = dlmread(filename,'',9,0); 
% 
% Top_Velocity_Field_Comsol=zeros(150,150);
% X_Comsol=zeros(150,150);
% Y_Comsol=zeros(150,150);
% for i=1:150
%     Top_Velocity_Field_Comsol(:,i)=flip(M((i-1)*150+1:150*i,3));   
%     X_Comsol(:,i)=M((i-1)*150+1:150*i,1);
%     Y_Comsol(:,i)=flip(M((i-1)*150+1:150*i,2));
% end

Top_Velocity_Field_Comsol=zeros(150,150);
X_Comsol=zeros(150,150);
Y_Comsol=zeros(150,150);
for i=1:150
    Top_Velocity_Field_Comsol(:,i)=flip(M((i-1)*150+1:150*i,3));   
    X_Comsol(:,i)=M((i-1)*150+1:150*i,1);
    Y_Comsol(:,i)=flip(M((i-1)*150+1:150*i,2));
end

subplot(1,2,1)
rho=998;
Boundary_Flux_Numerical=Top_Velocity_Field_Comsol/rho;
contourf(Boundary_Flux_Numerical)
colorbar
title('Fluxes of Numerical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 

subplot(1,2,2)
Boundary_Flux_Analytical=W_test_matrix;
contourf((Boundary_Flux_Analytical));
colorbar
title('Fluxes of Analytical Model','FontSize',40)
ax = gca;
ax.FontSize = 18; 
 
% subplot(1,2,2)
% Boundary_Flux_Numerical_OneFre=Top_Velocity_Field_Comsol_OneFre/rho;
% contourf((Boundary_Flux_Numerical_OneFre));
% colorbar
% title('Fluxes of Analytical Model','FontSize',40)
% ax = gca;
% ax.FontSize = 18; 
%  






