function [axf axw] = SpecPlot1D(f,P)

% [axf axw] = SpecPlot1D(f,P)
%
% Plots spectral power or periodogram P against frequency f and wavelength
% 1/f using two x-axes. Returns handles to the frequency and wavelength
% axes.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

loglog(f,P,'or','markersize',3); % Log-log scale plot;
axf = gca; % get current axes as axf;

frange = get(axf,'xlim'); % The limit of x-axis;
wrange = 1./fliplr(frange); % The flip of the inverse of x-axis;

axw=axes('Position',get(axf,'Position'),...% Line continuation;
         'XAxisLocation','top',...
         'YAxisLocation','right',...
         'Color','none',...
         'XColor','k','YColor','k',...
         'xlim',wrange,'xdir','reverse',...
         'ytick',[],'xscale','log');
% XAxisLocation=x-axis location, YAxisLocation=y-axis location;
% Color=Background color;
% XColor=x-axis grid color;
% YColor=y-axis grid color;
% xdir=x-axis direction;
% Tick values, specified as a vector of increasing values
set(axf,'box','off','tickdir','out') % box=axes outline; 
set(axw,'box','off','tickdir','out') % Tick mark direction, specified as one of these values:
%'in' ? Direct the tick marks inward from the axis lines. (Default for 2-D views)
%'out' ? Direct the tick marks outward from the axis lines. (Default for 3-D views)
%'both' ? Center the tick marks over the axis lines.

set(get(axf,'xlabel'),'string','Radial frequency (m^{-1})','FontSize',20)
xt=get(axf,'XTick');
set(axf,'FontSize',20);
set(get(axf,'ylabel'),'string','Mean squared amplitude (m^2)','FontSize',20)
set(get(axw,'xlabel'),'string','Wavelength (m)','FontSize',20) % set x-axis, y-axis's name
xt=get(axw,'XTick');
set(axw,'FontSize',20);
