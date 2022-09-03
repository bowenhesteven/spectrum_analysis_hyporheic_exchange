%%----------Create 3D synthetic surface and do the spectrum analysis-------%

% Instruction:
% The same algorithm used in (Chen et al., 2018) to simulate 3D bed form
% dunes, the input variables can be changed.
%% Use Rubin and Carter (2005) algorithm to simulate the 3D synthetic surface (Chen et al., 2018)
lambda_L=40;
lambda_T=60;
phi1_bf=30;
phi2_bf=45;
alpha1=60;
alpha2=60;
beta1=330;
beta2=beta1;
eta=0.053;
H=eta*lambda_L;
A1_f=3;
A1_s=9;
A2_f=4;
A2_s=6;
l=100;

% Calculate the m1
x=-50:50;
y=-50:50;
[X,Y]=meshgrid(x,y);
m1=X-lambda_L*(phi1_bf/360)-A1_f*sin((2*pi/lambda_T)*Y+(2*pi/360)*alpha1)-A1_s*sin((4*pi/lambda_T)*Y+(2*pi/360)*beta1);

Z_Surf1=(eta*lambda_L/l)*(7.5-6*sin(pi*l*m1/(50*lambda_L))-1.5*sin(pi*l*m1/(25*lambda_L)));

% Calculate the m2
x=-50:50;
y=-50:50;
[X,Y]=meshgrid(x,y);
m2=X-lambda_L*(phi2_bf/360)-A2_f*sin((2*pi/lambda_T)*Y+(2*pi/360)*alpha2)-A2_s*sin((4*pi/lambda_T)*Y+(2*pi/360)*beta2);
Z_Surf2=(eta*lambda_L/l)*(7.5-6*sin(pi*l*m2/(50*lambda_L))-1.5*sin(pi*l*m2/(25*lambda_L)));

Z_H=max(Z_Surf1,Z_Surf2);

surfl(X,Y,Z_H);
