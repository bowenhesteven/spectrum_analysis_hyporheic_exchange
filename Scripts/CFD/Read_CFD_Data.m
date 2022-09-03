%%%%%%%%% Read the CFD data %%%%%%%%%%
%% read the mean shear, mean depth, and mean velocity to the matlab, plot them
filename='mass2_10m_modern.xlsx';
CFD = xlsread(filename);
mean_shear=CFD(:,17);
mean_depth=CFD(:,7);
mean_velocity=CFD(:,12);
Easting=CFD(:,3);
Northing=CFD(:,4);
scatter(X,Y,1,mean_velocity,'filled');
xlabel('Easting','FontSize',20);
ylabel('Northing','FontSize',20);
title('CFD model-mean flow velocity (m/s)','FontSize',30);
%-----------------------------------------------------------------------------