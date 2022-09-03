% This function creates fractal landscapes by using diamond-square algorithm
% described by Paul Marts in his web page 
% http://www.gameprogrammer.com/fractal.html
% His code was in C something so I decided to write a much simpler one in
% MATLAB for fun. You can just run this function to create fractal
% landscapes. It doesn't depend on any toolbox so it should work any MATLAB
% version.
%
% You can also export the landscape to Povray by using the exportToPovray
% function at the end of the main function.
%
% If you want to increase the resolution, change maxstep parameter. 8 or 9
% would be fine.
% 
% author: Huseyin Kaya (http://www.youtube.com/user/TemelBilisim)
function fractal_landscape
close all
clc

figure('Position', [0, 0, 800, 800]);
% x = normrnd(0,10,[2,2]);
x = zeros(2,2); % Initial 2 by 2 zero matrix.
maxstep=9;
showit(x,1,maxstep,0,0,'Initial matrix');

h = 0.89; % [0 1], Close to zero results smooth landscapes.
range=1; 
for step=2:maxstep
    x = diamond_square_algorithm(x,range,step,maxstep);
    % In each step, the range for random numbers 
    range = range*h*2^(-h); 
end

% If you open this, then the final output matrix is converted to a POVRAY
% input file. If you have Povray installed, then you can render it.
exportToPovray(x)
end


% This is the diamond-square algorithm described by Paul Marts in his web page
% http://www.gameprogrammer.com/fractal.html
% My implemention focuses on readibility rather than performance so it
% might run slower than expected.
function y = diamond_square_algorithm(x,range,step,maxstep)

x = addborders(x);
y = enlarge(x);
y = centerofsquares(y,range,step,maxstep);
y = wrap(y);
y = centerofdiamonds(y,range,step,maxstep);
y = removeborders(y);

end

% I add a border around x matrix. So the new matrix has two extra rows and
% two extra columns with zeros. These borders are used to wrap the corners
% of diamonds which lies outside of the input matrix.
%
% Example:
%
%   x             y
% -----       ---------
%             0 0 0 0 0
% a b c       0 a b c 0
% d e f  ->   0 d e f 0
% g h i       0 g h i 0
%             0 0 0 0 0
%
function y = addborders(x)
y = zeros(size(x)+2);
y(2:end-1,2:end-1)=x;
end

% In each step of diamond-square algorithm, the size of the matrix
% increases to accommodate newly calculated center of squares. So just
% before the calculation of the center of squares the matrix should be
% enlarged. Remember that input matrix has also borders.
%
% Example:
%
%     x                y
% ---------      -------------
%                0 0 0 0 0 0 0     
% 0 0 0 0 0      0 a 0 b 0 c 0
% 0 a b c 0      0 0 0 0 0 0 0
% 0 d e f 0  ->  0 d 0 e 0 f 0
% 0 g h i 0      0 0 0 0 0 0 0
% 0 0 0 0 0      0 g 0 h 0 i 0
%                0 0 0 0 0 0 0
%
function y = enlarge(x)
n = size(x,1);
y = zeros(2*n-3,2*n-3);
for i=2:n-1
    for j=2:n-1
        y(2*i-2,2*j-2)=x(i,j);
    end
end
end

% After enlarging the matrix, now we can calculate the center of the
% squares by using the surrounding four corners. The height of the centers
% is shifted with a random number between [-range +range] interval.
%
% Example: "q","x","y" and "z" represent the calculated center of squares.
%
%       x                  y
% -------------      -------------
% 0 0 0 0 0 0 0      0 0 0 0 0 0 0
% 0 a 0 b 0 c 0      0 a 0 b 0 c 0
% 0 0 0 0 0 0 0      0 0 q 0 x 0 0
% 0 d 0 e 0 f 0  ->  0 d 0 e 0 f 0
% 0 0 0 0 0 0 0      0 0 y 0 z 0 0
% 0 g 0 h 0 i 0      0 g 0 h 0 i 0
% 0 0 0 0 0 0 0      0 0 0 0 0 0 0
function y = centerofsquares(x,range,step,maxstep)
n=size(x,1);
y = x;
for i=3:2:n-2
    for j=3:2:n-2
        y(i,j) = (y(i-1,j-1)+y(i+1,j+1)+y(i+1,j-1)+y(i-1,j+1))*0.25 + 2*range*rand-range; % random number between [-range +range]
        showit(removeborders(y),step,maxstep,i,n-2,'Square ');
    end
end
end

% Before going into diamonds phase, the corners of the prospected diamonds
% which lies on the borders should be warped so that we are able to
% calculate the center of diamons.
%
% Example: 
%
% Row       x                  y
% ---  -------------      -------------
% 1    0 0 0 0 0 0 0      0 0 y 0 z 0 0   (row 5 [end-2])
% 2    0 a 0 b 0 c 0      0 a 0 b 0 c 0
% 3    0 0 q 0 x 0 0      x 0 q 0 x 0 q
% 4    0 d 0 e 0 f 0      0 d 0 e 0 f 0
% 5    0 0 y 0 z 0 0      z 0 y 0 z 0 y
% 6    0 g 0 h 0 i 0      0 g 0 h 0 i 0
% 7    0 0 0 0 0 0 0      0 0 q 0 x 0 0   (row 3 [3])
%                         |           |
%                         |           |___(col 3 [3])
%                         |___ (col 5 [end-2])
function y=wrap(x)
y = x;
y(:,1)   = x(:,end-2); 
y(:,end) = x(:,3)    ; 
y(1,:)   = x(end-2,:);
y(end,:) = x(:,3)    ;
end

% Now we are ready to calculate the center of diamonds. It is very similar
% to the center of squares except that we now use the point up, down, left
% and right i.e. horizontal and vertical locations.
%
% Example: "." represents the center of diamonds.
%
%      x                  y
% -------------      -------------
% 0 0 y 0 z 0 0      0 0 y 0 z 0 0
% 0 a 0 b 0 c 0      0 a . b . c 0
% x 0 q 0 x 0 q      x . q . x . q
% 0 d 0 e 0 f 0  ->  0 d . e . f 0
% z 0 y 0 z 0 y      z . y . z . y
% 0 g 0 h 0 i 0      0 g . h . i 0
% 0 0 q 0 x 0 0      0 0 q 0 x 0 0
%
function y = centerofdiamonds(x,range,step,maxstep)
n=size(x,1);
y = x;
for i=2:1:n-1
    for j=3-mod(i,2):2:n-1
        y(i,j) = (y(i-1,j)+y(i+1,j)+y(i,j-1)+y(i,j+1))*0.25+2*range*rand-range;
        showit(removeborders(y),step,maxstep,i,n-1,'Diamond');
    end
end
end

% Now we can remove the borders.
%
% Example:
%
%      x                   y 
% -------------        ---------
% 0 0 y 0 z 0 0      
% 0 a . b . c 0        a . b . c
% x . q . x . q        . q . x .
% 0 d . e . f 0   ->   d . e . f
% z . y . z . y        . y . z .
% 0 g . h . i 0        g . h . i
% 0 0 q 0 x 0 0     
%
function y = removeborders(x)
y=x(2:end-1,2:end-1);
end

% This function creates a povray input file by using the landscape
% represented by the matrix x.
function exportToPovray(x)
nm=5/size(x,1);
fd=fopen('fractal_landscape.pov','w+');
fprintf(fd,'#include "colors.inc"\n');
fprintf(fd,'#include "stones.inc"\n');
fprintf(fd,'#include "textures.inc"\n');
fprintf(fd,'#include "shapes.inc"\n');
fprintf(fd,'#include "glass.inc"\n');
fprintf(fd,'#include "metals.inc"\n');
fprintf(fd,'#include "woods.inc"\n');
fprintf(fd,'#declare c1 = texture { T_Stone1 scale 1 }\n');
fprintf(fd,'light_source { <0, 2, 2> color White}\n');
fprintf(fd,'camera { location <3, 1, -3> look_at  <3,0,3>  }\n');
fprintf(fd,'mesh {\n');
for i=1:size(x,1)-1
    for j=1:size(x,1)-1
        fprintf(fd,'triangle{<%8.5f,%8.5f,%8.5f>,<%8.5f,%8.5f,%8.5f>,<%8.5f,%8.5f,%8.5f> texture {c1}}\n',...
            i*nm,x(i,j),j*nm,...
            (i+1)*nm,x(i+1,j),j*nm,...
            (i+1)*nm,x(i+1,j+1),(j+1)*nm);
        fprintf(fd,'triangle{<%8.5f,%8.5f,%8.5f>,<%8.5f,%8.5f,%8.5f>,<%8.5f,%8.5f,%8.5f> texture {c1}}\n',...
            i*nm,x(i,j),j*nm,...
            i*nm,x(i,j+1),(j+1)*nm,...
            (i+1)*nm,x(i+1,j+1),(j+1)*nm);
    end
end
fprintf(fd,'}\n');
end


% This function is used to make a nice animation during the algorithm. It
% is called by square and diamond steps in their inner loops so that we can
% visualize what is going on.
function showit(x,step,maxstep,scan,maxscan,phase)
y = zeros(2^(maxstep-1)+1);
insertlocations = linspace(1,2^(maxstep-1)+1,2^(step-1)+1);
y(insertlocations,insertlocations)=x;
subplot(211)
surf(y);
colormap gray
zlim([-1 1])
title(sprintf('Surface Plot, Step: %d/%d Scan: %3d/%3d, Phase: %9s',step,maxstep,scan,maxscan,phase));
set(gca,'box','off');
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'ztick',[]);
axis off;
view(60,50)

subplot(212)
imagesc(y')
axis equal
set(gca,'box','off');
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'ztick',[]);
axis off;colormap gray
title(sprintf('Projection  , Step: %d/%d Scan: %3d/%3d, Phase: %9s',step,maxstep,scan,maxscan,phase));
pause(0.01)
end