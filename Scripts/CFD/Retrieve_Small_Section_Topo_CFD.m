%%% Retrieve a small piece/section from both topography data, CFD data %%%%
%% Cut the small section from the topo;
[I,map] = imread('hanlidar_nad88.tif',1); % read the index image which is the first image of the file;
imshow(I); % show the image;
hold on

% The coordinate of the rectangle of the image;
xmin=1.40e+04; % The min x-direction coordinate;
height=149; % The height of the section;
ymin=9300; % The min y-direction coordinate;
width=149; % The width of the section;

% cut the piece from the image;
topopiece_vector=cut_piece(I,xmin,ymin,width,height);

% Plot the piece of the topo, marker the coodinates of the rectangle/the cutted section;
plot(xmin,ymin+height,'r+', 'MarkerSize', 10)
hold on
plot(xmin,ymin,'r+', 'MarkerSize', 10)
hold on
plot(xmin+width,ymin,'r+', 'MarkerSize', 10)
hold on
plot(xmin+width,ymin+height,'r+', 'MarkerSize', 10)
%% Cut the small section from the CFD model;
% Mapping the CFD model to matlab matrix;

Easting_min=min(Easting);
Northing_min=min(Northing);
Easting_max=max(Easting);
Northing_max=max(Northing);

% Find the real Easting and Northing position;
Easting_Section_min=Easting_min+((Easting_max-Easting_min)/51786)*(xmin-1);
Easting_Section_max=Easting_min+((Easting_max-Easting_min)/51786)*(xmin+width-1);
Northing_Section_min=Northing_max-((Northing_max-Northing_min)/51714)*(ymin-1);
Northing_Section_max=Northing_max-((Northing_max-Northing_min)/51714)*(ymin+height-1);

% for x=xmin:xmin+width
%  for y=ymin:ymin+height  
%  x_map(i)=Easting_min+(xmin-1)*((Easting_max-Easting_min)/51786);
%  y_map(i)=Northing_min+(ymin-1)*((Northing_max-Northinging_min)/51714);
%  % shear_matrix(y_map,x_map)=mean_shear(i);
%  end
% end
% Plot the piece of the topo, marker the coodinates of the rectangle/the cutted section;
scatter(Easting,Northing,1,mean_velocity,'filled');
xlabel('Easting','FontSize',20);
ylabel('Northing','FontSize',20);
title('CFD model-mean flow velocity (m/s)','FontSize',30);
hold on;

plot(Easting_Section_min,Northing_Section_min,'r+', 'MarkerSize', 10)
hold on
plot(Easting_Section_min,Northing_Section_max,'r+', 'MarkerSize', 10)
hold on
plot(Easting_Section_max,Northing_Section_min,'r+', 'MarkerSize', 10)
hold on
plot(Easting_Section_max,Northing_Section_max,'r+', 'MarkerSize', 10)


