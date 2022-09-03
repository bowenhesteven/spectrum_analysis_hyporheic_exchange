%% Create several synthetic surfaces with different number of frequencies 
% in each direction to see how different number of Fourier fitting
% fequencies influence the Fouier fitted surface

mirror=1; 
periodic=0;
% mirror=1;periodic=1;
f1=3.5; % f1 = a coefficient for distribution of wavelengths (generally 1<f1<2.5)
f3=1+mirror;
ratio_high=300;
ratio_low=1;
cross=0;
%% The synthetic surface with one frequency; analysis
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z1=30*sin((2*pi/10)*X).*sin((2*pi/10)*Y);
surfl(x,y,z1);colormap copper
shading interp
%% The synthetic surface with two frequency
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z2=30*sin((2*pi/10)*X).*sin((2*pi/10)*Y)+45*sin((2*pi/5)*X).*sin((2*pi/5)*Y);
surfl(x,y,z2);colormap copper
shading interp
%% The synthetic surface with four frequency
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z4=30*sin((2*pi/10)*X).*sin((2*pi/10)*Y)+45*sin((2*pi/5)*X).*sin((2*pi/5)*Y)+20*sin((2*pi/15)*X).*sin((2*pi/15)*Y)+5*sin((2*pi/8)*X).*sin((2*pi/8)*Y);
surfl(x,y,z4);colormap copper
shading interp
%% The synthetic surface with eight frequency
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z8=30*sin((2*pi/10)*X).*sin((2*pi/10)*Y)+45*sin((2*pi/5)*X).*sin((2*pi/5)*Y)+20*sin((2*pi/15)*X).*sin((2*pi/15)*Y)+5*sin((2*pi/8)*X).*sin((2*pi/8)*Y)+8*sin((2*pi/25)*X).*sin((2*pi/25)*Y)+12*sin((2*pi/45)*X).*sin((2*pi/45)*Y)+36*sin((2*pi/6)*X).*sin((2*pi/6)*Y)+1*sin((2*pi/80)*X).*sin((2*pi/80)*Y);
surfl(x,y,z8);colormap copper
shading interp
%% The synthetic surface with sixteen frequency
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z16=30*sin((2*pi/10)*X).*sin((2*pi/10)*Y)+45*sin((2*pi/5)*X).*sin((2*pi/5)*Y)+20*sin((2*pi/15)*X).*sin((2*pi/15)*Y)+5*sin((2*pi/8)*X).*sin((2*pi/8)*Y)+8*sin((2*pi/25)*X).*sin((2*pi/25)*Y)+12*sin((2*pi/45)*X).*sin((2*pi/45)*Y)+36*sin((2*pi/6)*X).*sin((2*pi/6)*Y)+1*sin((2*pi/80)*X).*sin((2*pi/80)*Y)...
    +2*sin((2*pi/13)*X).*sin((2*pi/13)*Y)+54*sin((2*pi/7)*X).*sin((2*pi/7)*Y)+9*sin((2*pi/51)*X).*sin((2*pi/51)*Y)+50*sin((2*pi/18)*X).*sin((2*pi/18)*Y)+61*sin((2*pi/42)*X).*sin((2*pi/42)*Y)+20*sin((2*pi/3)*X).*sin((2*pi/3)*Y)+8*sin((2*pi/60)*X).*sin((2*pi/60)*Y)+11*sin((2*pi/32)*X).*sin((2*pi/32)*Y);
surfl(x,y,z16);colormap copper
shading interp
%% The synthetic surface with 32 frequency
x=1001:1150;
y=1001:1150;
y=fliplr(y);
[X,Y]=meshgrid(x,y);
z32=30*sin((2*pi/10)*X).*sin((2*pi/10)*Y)+45*sin((2*pi/5)*X).*sin((2*pi/5)*Y)+20*sin((2*pi/15)*X).*sin((2*pi/15)*Y)+5*sin((2*pi/8)*X).*sin((2*pi/8)*Y)+8*sin((2*pi/25)*X).*sin((2*pi/25)*Y)+12*sin((2*pi/45)*X).*sin((2*pi/45)*Y)+36*sin((2*pi/6)*X).*sin((2*pi/6)*Y)+1*sin((2*pi/80)*X).*sin((2*pi/80)*Y)...
    +2*sin((2*pi/13)*X).*sin((2*pi/13)*Y)+54*sin((2*pi/7)*X).*sin((2*pi/7)*Y)+9*sin((2*pi/51)*X).*sin((2*pi/51)*Y)+50*sin((2*pi/18)*X).*sin((2*pi/18)*Y)+61*sin((2*pi/42)*X).*sin((2*pi/42)*Y)+20*sin((2*pi/3)*X).*sin((2*pi/3)*Y)+8*sin((2*pi/60)*X).*sin((2*pi/60)*Y)+11*sin((2*pi/32)*X).*sin((2*pi/32)*Y)...
    +3*sin((2*pi/4)*X).*sin((2*pi/4)*Y)+5*sin((2*pi/6)*X).*sin((2*pi/6)*Y)+7*sin((2*pi/48)*X).*sin((2*pi/48)*Y)+54*sin((2*pi/100)*X).*sin((2*pi/100)*Y)+19*sin((2*pi/88)*X).*sin((2*pi/88)*Y)+19*sin((2*pi/56)*X).*sin((2*pi/56)*Y)+64*sin((2*pi/63)*X).*sin((2*pi/63)*Y)+53*sin((2*pi/89)*X).*sin((2*pi/89)*Y)...
    +21*sin((2*pi/78)*X).*sin((2*pi/78)*Y)+52*sin((2*pi/91)*X).*sin((2*pi/91)*Y)+93*sin((2*pi/77)*X).*sin((2*pi/77)*Y)+3*sin((2*pi/9)*X).*sin((2*pi/9)*Y)+65*sin((2*pi/45)*X).*sin((2*pi/45)*Y)+22*sin((2*pi/23)*X).*sin((2*pi/23)*Y)+81*sin((2*pi/66)*X).*sin((2*pi/66)*Y)+16*sin((2*pi/36)*X).*sin((2*pi/36)*Y);
surfl(x,y,z32);colormap copper
shading interp
%% The real natural surface with the size of the synthetic surface
I = imread('hanlidar_nad88.tif'); % read the index image which is the first image of the file;
imshow(I); % show the image;
hold on

% The coordinate of the rectangle of the image;
xmin=1.40e+04; % The min x-direction coordinate;
height=149; % The height of the section;
ymin=9300; % The min y-direction coordinate;
width=149; % The width of the section;

% cut the piece from the image;
topopiece_vector=cut_piece(I,xmin,ymin,width,height);
surf(topopiece_vector.x,topopiece_vector.y,topopiece_vector.z);
colorbar;
title('A small section of the real topograhy surface','FontSize',30);
xlabel('x (m)','FontSize',25);
ylabel('y (m)','FontSize',25);
surfl(topopiece_vector.x,topopiece_vector.y,topopiece_vector.z);colormap copper
shading interp
title('A small section of the real topograhy surface','FontSize',30);
xlabel('x (m)','FontSize',25);
ylabel('y (m)','FontSize',25);
zlabel('Elevation (m)', 'FontSize',25);
% Plot the piece of the topo, marker the coodinates of the rectangle/the cutted section;
plot(xmin,ymin+height,'r+', 'MarkerSize', 30)
hold on
plot(xmin,ymin,'r+', 'MarkerSize', 30)
hold on
plot(xmin+width,ymin,'r+', 'MarkerSize', 30)
hold on
plot(xmin+width,ymin+height,'r+', 'MarkerSize', 30)
%% Plot the result (influence of Fourier fitting frequency on Fourier fitted surface)
mean_error_matrix=zeros(1,7);
N_WL=[2,4,8,16,32,64,128];
for i=1:7
    N_f=2^i;
    [paramhat,lambda,N,kx,ky,h,mean_error,av]=Spectop(x,y,topopiece_vector.z,N_f,f1,f3,mirror,ratio_low,ratio_high,cross,periodic);
    mean_error_matrix(i)=mean_error;
end
plot(N_WL,mean_error_matrix,'-o',...
    'LineWidth',3,...mea 
    'MarkerEdgeColor','b',...
    'MarkerSize',10)
title('Influence of Number of Frequencies','FontSize',40)
ax = gca;
ax.FontSize = 18; 
xlabel('Number of frequencies in each direction','FontSize',20)
ylabel('Mean error','FontSize',20)

