clear all;
%cd '/scratch/heb1'; % cd to the image directory;
[hanlidar_aspect_indexed,map] = imread('hanlidar_nad88_aspect.tif',1); % read the index image which is the first image of the file;
topo_aspect=hanlidar_aspect_indexed; % save the index into the matrix Index;
%topo_aspect(topo_aspect==0)=nan;
%[nr,nc]=size(topo_aspect);
%x=1:nc;
%y=1:nr; y=fliplr(y);
%[X,Y]=meshgrid(x,y);
%goods = ~isnan(topo);
zp_aspect = topo_aspect(goods);
topo_aspect_vector = struct('x',topoVector_coordinates.x,'y',topoVector_coordinates.y,'z',zp_aspect);
%topoVector_coordinates=struct('x',xp,'y',yp);
%filedirectory='/scratch/heb1'; % save the result to this directory;
save('?topoVector_aspect.mat','topo_aspect_vector', '-v7.3');
%save('topoVector_coordinates.mat','','-v7.3');
plot(topoVector_coordinates.x,topoVector_coordinates.y,'.r');