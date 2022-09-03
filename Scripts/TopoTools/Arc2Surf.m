function Arc2Surf(agrid,sgrid)

% converts ESRI ASCII grid agrid.asc to Surfer ASCII grid sgrid.grd.
% agrid and sgrid should be character strings

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

% Read Arc grid
M = ReadArcGrid(agrid);

% Write Surfer grid
WriteSurferGrid(M,sgrid)

% delete the variable
clear M