%%%%%%%%%%%%%%%%%%%%%%%%% Rotation Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function extracts a piece small section (rectanglar) from the topo_vector.mat
% (based on readOriginalImage.m) and rotate the section to the normal
% x-diretion/horizontal and y-direction/vertical based on the river flow direction angle theta;

% Input:
% x1: the x value of the first point/point1 which is in the middle of the rectangular's width;
% y1: the y value of the first point/point1;
% x2: the x value of the second point/point2, is the opposite point on the
% other side of the rectanglar, also the middle point of the rectangular's width;
% y2: the y value of the second point/point2, is the opposite point on the
% other side of the rectanglar;
% half_width: define the half width of the rectanglar section, the length
% of the rectanglar is defined by the Euclidean distance of the two points;

% Output:
% topo_vector_rotation: A struct consists of rotated x coordinate, rotated
% y coordinate and the z-elecation associated with those corresponding x
% and y coordinates;

function[topo_vector_rotation] = Rotation(x1,y1,x2,y2,half_width,topo_vector)

% Example to test;
% imshow(I); % show the image;
% hold on;

% x1=24500;
% y1=(51714+1)-8000;
% 
% x2=22520;
% y2=(51714+1)-10000;
% 
% half_width=100;


% length=sqrt((x1-x2)^2+(y1-y2)^2); % Calculate the length of the rectanglar based on the 2 points' position;
tang_theta=(y1-y2)/(x1-x2); % Calculate the tangent value of the angular theta;
theta=atand(tang_theta); % Calculate the theta based on the tang_theta we just calculated;
thetaInRadians = deg2rad(theta);
theta1=thetaInRadians+pi/2; % Calculate another theta which is perpendicular to the theta;

% verify the product of the two tangent theta equals -1;
product_two_tangent=tan(theta1)*tang_theta; % this value should be -1;

% Start to calculate/generate the four points of the rectangular/polygon;
% The four points of the rectangular labels a, b, c, d in an anticlockwise direction;

% a point;
Xa=x1-sin(thetaInRadians)*half_width;
Ya=y1+cos(thetaInRadians)*half_width;

% b point;
Xb=x2-sin(thetaInRadians)*half_width;
Yb=y2+cos(thetaInRadians)*half_width;

% c point;
Xc=sin(thetaInRadians)*half_width+x2;
Yc=y2-cos(thetaInRadians)*half_width;

% d point;
Xd=sin(thetaInRadians)*half_width+x1;
Yd=y1-cos(thetaInRadians)*half_width;

% plot(Xa,2*y1-Ya,'r+', 'MarkerSize', 20)
%  hold on
%  plot(Xb,2*y2-Yb,'r+', 'MarkerSize', 20)
%  hold on
%  plot(Xc,2*y2-Yc,'r+', 'MarkerSize', 20)
%  hold on
%  plot(Xd,2*y1-Yd,'r+', 'MarkerSize', 20)


% Create polygons/plot the filled polygons;
xv=[Xa Xb Xc Xd];
yv=[Ya Yb Yc Yd];
patch(xv,yv,'red');

% Decide the points that locates inside or on edge of polygonal region;
xq=topo_vector.x; % query points/x coordinate;
yq=topo_vector.y; % query points/y coordinate;
zq=topo_vector.z;
[in,on]=inpolygon(xq,yq,xv,yv);

% Generate the points that locates inside or on edge of polygonal region;
x_in=xq(in);
x_on=xq(on);


y_in=yq(in);
y_on=yq(on);

z_in=zq(in);
z_on=zq(on);


x_true=[x_in;x_on];
y_true=[y_in;y_on];
z_true=[z_in;z_on];


% Do the rotation calculation to the chosen points;
Old_Matrix=[x_true';y_true'];

scatter(x_true',y_true');
hold on;

R=[cos(2*pi-thetaInRadians) -sin(2*pi-thetaInRadians);sin(2*pi-thetaInRadians) cos(2*pi-thetaInRadians)];
New_Matrix=R*Old_Matrix;
x_rotated=New_Matrix(1,:);
y_rotated=New_Matrix(2,:);

scatter(x_rotated,y_rotated);


topo_vector_rotation = struct('x',x_rotated','y',y_rotated','z',z_true);
end

