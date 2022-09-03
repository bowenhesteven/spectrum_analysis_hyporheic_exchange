%%%% Cut the piece in order to test conformal mapping function %%%%

%%% Read the image  %%%

I = imread('hanlidar_nad88.tif'); % read the index image which is the first image of the file;
imshow(I); % show the image;
hold on

%%% Cut the section that contains meander part %%%
% The coordinate of the rectangle of the image;
xmin=2.792e+04; % The min x-direction coordinate;
height=4999; % The height of the section;
ymin=1; % The min y-direction coordinate;
width=6799; % The width of the section;

% cut the piece from the image;
topopiece_vector=cut_piece(I,xmin,ymin,width,height);

% Plot the piece of the topo, marker the coodinates of the rectangle/the cutted section;
plot(xmin,ymin+height,'r+', 'MarkerSize', 30)
hold on
plot(xmin,ymin,'r+', 'MarkerSize', 30)
hold on
plot(xmin+width,ymin,'r+', 'MarkerSize', 30)
hold on
plot(xmin+width,ymin+height,'r+', 'MarkerSize', 30)
%% Set up the Matlab desktop
% Make new figures docked and white
set(0,'DefaultFigureWindowStyle','docked','DefaultFigureColor','w')
%% Import gridded data
dx = abs(topopiece_vector.x(2) - topopiece_vector.x(1)); % grid spacing in the x-direction
dy = abs(topopiece_vector.y(2) - topopiece_vector.y(1)); % grid spacing in the y-direction
[Ny,Nx] = size(topopiece_vector.z); % grid dimensions
%% View data
% Display a shaded relief map (Fig. 2)
figure('Name','Fig. 2: Shaded relief','NumberTitle','off')
ShadePlot(topopiece_vector.x,topopiece_vector.y,topopiece_vector.z) % ShadePlot of the topo/surface;
hold on;
plot(6000,1,'r+', 'MarkerSize', 30)
hold on
plot(6800,2500,'r+', 'MarkerSize', 30)
hold on
plot(4000,4700,'r+', 'MarkerSize', 30)
hold on
plot(2000,5000,'r+', 'MarkerSize', 30)
hold on;
plot(833,3611,'r+', 'MarkerSize', 30)
hold on
plot(1,2500,'r+', 'MarkerSize', 30)
hold on
plot(1,1,'r+', 'MarkerSize', 30)
hold on
plot(833,1,'r+', 'MarkerSize', 30)
hold on;
plot(2000,1700,'r+', 'MarkerSize', 30)
hold on
plot(4000,2200,'r+', 'MarkerSize', 30)
hold on
plot(4000,600,'r+', 'MarkerSize', 30)
hold on;
x=[6000 6800 4000 2000 833 1 1 833 2000 4000 4000 6000];
y=[1 2500 4700 5000 3611 2500 1 1 1700 2200 600 1];
line(x,y,'LineStyle', '-','Color','r');
%% Retrieve the coordinates
[X,Y]=meshgrid(topopiece_vector.x,topopiece_vector.y);
topopiece_vector.z(topopiece_vector.z==0) = nan; 
Index=~isnan(topopiece_vector.z); % The logic index of the 2D filtered power spectrum;
X_Coordinate=X.*Index;
Y_Coordinate=Y.*Index;
% Create valid X_coordiante vector
X_Coordinate_vector = nonzeros(X_Coordinate);
Y_Coordinate_vector=nonzeros(Y_Coordinate);
Coordinate=[X_Coordinate_vector,Y_Coordinate_vector];
%% Test Conformal Mapping on this Meander Section %%
 scgui;
zp=[6800+1i 6800+2500i 4000+4700i 2000+5000i 833+3611i 1+2500i 1+1i 833+1i 2000+1700i 4000+2200i 4000+600i];
p=polygon(zp);
plot(p);
f=diskmap(p);
rectmap(p);
[h,val1,val2]=plot(ans);
h=plot(f,val1,val2);
















