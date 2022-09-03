function ax = SpecPlot2D(f,P)

% ax = SpecPlot2D(f,P)
%
% Plots a two-dimensional power spectrum with frequency matrix f and power
% spectrum P. Returns a handle to the axes.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

ax = axes;
% imagesc(P); axis image
% j = jet; % j(1,:) = [ 1 1 1 ];
% colormap(j);
%  colormap(grayscale);
% colorbar;
% imagesc(C) displays the data in array C as an image that uses the full 
% range of colors in the colormap. Each element of C specifies the color 
% for 1 pixel of the image. The resulting image is an m-by-n grid of pixels 
% where m is the number of rows and n is the number of columns in C. 
 contour(P);
 contourcbar;
xlabel('x frequency','FontSize',20) % write the x frequency;
ylabel('y frequency','FontSize',20) % write the y frequency;

[nfy , nfx] = size(f); % 1024*1024;
nyq = f(nfy/2+1,1); % the Nyquist frequency in the x direction

% Note that these frequency labels are approximate due to the offset of the
% DC (zero frequency) element. The frequencies in f are exact.
numticks=30;
set(gca,'XTick',linspace(1,nfx,numticks)) % linspace(x1,x2,n) generates n points. 
                                          % The spacing between the points is (x2-x1)/(n-1).
xt=get(gca,'XTick');
set(gca,'FontSize',20);
set(gca,'YTick',linspace(1,nfy,numticks))
xt=get(gca,'YTick');
set(gca,'FontSize',20);
set(gca,'XTicklabel',linspace(-nyq,nyq,numticks)) % xticklabels(labels) sets 
                                                  % the x-axis tick labels 
                                                  % for the current axes. 
                                                  % Specify labels as a string 
                                                  % array or a cell array of character vectors; 
                                                  % for example, {'January','February','March'}. 
                                                  % If you specify the labels, then the x-axis tick values and tick labels no longer update automatically based on changes to the axes.
set(gca,'YTicklabel',linspace(-nyq,nyq,numticks))
set(gca,'TickDir','out')
