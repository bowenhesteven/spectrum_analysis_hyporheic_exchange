function [A, B, D] = D8(M,flood)

% [A, B, D] = D8(M,flood)
% 
% Computes upslope contributing area using the D8 (steepest descent) 
% algorithm of O'Callaghan & Mark (1984)
% 
%      Input arguments:
%        M      a struct variable containing at least the fields
%        
%               M.grid: a K x J matrix of elevations
%               M.dx:   grid spacing in x-direction
%               M.dy:   grid spacing in y-direction
%        flood  if 1, routes flow through local minima in M (potentially
%               time-consuming). if zero, skips this step (faster).
%               (optional: default = 0)
% 
%      Return arguments:
%        A      a matrix of total contributing areas (in units of dx*dy)
%        B      a matrix of cells that receive drainage from boundaries
%        D      a matrix of D8 drainage directions (1-8, CCW from east) 
 
% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.
 

if (nargin < 1 || nargin > 2 || nargout > 3), help(mfilename), return, end

if nargin == 1, flood = 0; end

% convert integer data to double precision if necessary
if isa(M.grid,'integer')
    M.grid = double(M.grid);
end

% Replace NaNs with ArcInfo NODATA value (-9999)
ndv = -9999;
M.grid(isnan(M.grid)) = ndv;

% Calculate drainage area
[A, B, D] = mexD8(M.grid,M.dy/M.dx,flood,ndv);

% Convert to true area units
A = A*M.dy*M.dy;
