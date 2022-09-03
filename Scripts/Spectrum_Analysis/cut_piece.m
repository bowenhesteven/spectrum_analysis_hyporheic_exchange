function[topopiece_vector]=cut_piece(I,xmin,ymin,width,height)
rect=[xmin ymin width height];
J = imcrop(I,rect);
topo_piece=J; % The elevation of the cropped area;
[nr,nc]=size(topo_piece); % Get the x coordinate and y coordinate size of the cropped area of image;
x = 1:nc;
y = 1:nr; y = fliplr(y);
[X,Y] = meshgrid(x,y);
topopiece_vector = struct('x',X(1,:),'y',Y(:,1),'z',topo_piece); % Reconstruct the index of the cropped area of the image;