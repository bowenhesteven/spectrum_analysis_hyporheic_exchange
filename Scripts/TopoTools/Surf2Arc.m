function Surf2Arc(sgrid,agrid)

% converts Surfer ASCII grid sgrid.grd to ESRI ASCII grid agrid.asc.
% sgrid and agrid should be character strings
%
% example:
%   Surf2Arc('gridname.grd','gridname.asc');

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

% Read Surfer grid
M = ReadSurferGrid(sgrid);

% Write Arc grid
WriteArcGrid(M,agrid);

% delete the variable
clear M