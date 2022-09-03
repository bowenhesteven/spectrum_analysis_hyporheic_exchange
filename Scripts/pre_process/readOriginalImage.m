% This function reads the images and

function[] = readOriginalImage(pathImage,pathOut)

pathImage = '/Users/gomezvjd/Downloads/dataBedforms/hanlidar_nad88.tif';
% pathOut = '/Users/gomezvjd/Downloads/dataBedforms/';



%% Read image
I = imread(pathImage,1);
I2 = imcrop(I,[3e4 0 1e4 2e4]);

figure(1)
subplot(1,2,1)
image(I)
subplot(1,2,2)
image(I2)

I2new = I2;
bads = I2new == 0;
I2new(bads) = nan;
imagesc(I2new)
axis 'equal'
info = imfinfo(pathImage)

%% Store coordinates with topographic information
goods = I >0;
[nr,nc] = size(I);
x = 1:nc;
y = 1:nr; y = fliplr(y);
j = 1:nc;
i = 1:nr; 


[X,Y] = meshgrid(x,y);
[C,R] = meshgrid(j,i);

xp = X(goods);
yp = Y(goods);
zp = I(goods);
cp = C(goods);
rp = R(goods);

topo_vector = struct('x',xp,'y',yp,'z',zp,'c',cp,'r',rp);

%% Create a subset of points for testing
fracPoints = 0.0005;
npoints = length(xp);

ind_subset = 1:npoints;
s = RandStream('mlfg6331_64'); 
ind_subset = datasample(s,ind_subset,round(npoints*fracPoints),'Replace',false);

xps = xp(ind_subset);
yps = yp(ind_subset);
zps = zp(ind_subset);
cps = cp(ind_subset);
rps = rp(ind_subset);

topo_vector_subset = struct('x',xps,'y',yps,'z',zps,'c',cps,'r',rps);


%plot(xp,yp,'.r')
%plot(xps,yps,'.r')

%% Save final vector
save(strcat(pathOut,'topoVector.mat'),'topo_vector', '-v7.3');
save(strcat(pathOut,'topoVector_subset.mat'),'topo_vector_subset', '-v7.3');


end

% pathData = '/Users/gomezvjd/Downloads/IndexTopography.mat';
%  
% data = load(pathData);
% topo = data.hanlidar_indexed;
% topo(topo == 0) = nan;
%  
% [nr,nc] = size(topo);
% x = 1:nc;
% y = 1:nr; y = fliplr(y)
%  
% [X,Y] = meshgrid(x,y);
%  
% goods = ~isnan(topo);
% xp = X(goods);
% yp = Y(goods);
% zp = topo(goods);
%  
% topo_vector = struct('x',xp,'y',yp,'z',zp);
% save('/Users/gomezvjd/Downloads/topoVector.mat','topo_vector', '-v7.3');
%  
% plot(xp,yp,'.r')
%  