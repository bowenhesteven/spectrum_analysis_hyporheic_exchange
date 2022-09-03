clear all;
%cd '/scratch/heb1'; % cd to the image directory;
[hanlidar_indexed,map] = imread('hanlidar_nad88.tif',1); % read the index image which is the first image of the file;
%imshow(hanlidar_indexed,map);
%imshow('hanlidar_nad88.tif','Reduce',true,'InitialMagnification','fit');
topo=hanlidar_indexed; % save the index into the matrix Index;
topo(topo==0)=nan;
[nr,nc]=size(topo);
x=1:nc;
y=1:nr; y=fliplr(y);
[X,Y]=meshgrid(x,y);
goods = ~isnan(topo);
xp = X(goods);
yp = Y(goods);
zp = topo(goods);
topo_vector = struct('x',xp,'y',yp,'z',zp);
topoVector_coordinates=struct('x',xp,'y',yp);
%filedirectory='/scratch/heb1'; % save the result to this directory;
save('topoVector.mat','topo_vector', '-v7.3');
save('topoVector_coordinates.mat','topoVector_coordinates','-v7.3');
save('topoVector_goods.mat','goods','-v7.3');
plot(xp,yp,'.r'); % Only plot the position with positive elevation;


%imageinfo('hanlidar_nad88.tif')