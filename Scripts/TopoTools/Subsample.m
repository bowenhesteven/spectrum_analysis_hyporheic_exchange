function Msub = Subsample(M,rows,cols)

% Msub = Subsample(M,rows,cols) 
%
% Subsamples the grid stored in the struct array M by extracting the values
% with row and column coordinates rows and cols. Returns the subsampled
% struct array Msub with updated values for the fields in M, which may
% include:
%
% M.grid      (the matrix of grid values)
% M.ncols     (# of columns in grid)
% M.nrows     (# of rows in grid)
% M.x         (x coordinates of centers of pixels)
% M.y         (y coordinates of centers of pixels)
% M.xllcorner (x coordinate of lower-left element)
% M.yllcorner (y coordinate of lower-left element)
% M.dx        (grid spacing in x direction)
% M.dy        (grid spacing in y direction)
%
% If the subsampled points are unevenly spaced, Msub.dx and Msub.dy will
% contain the mean spacing.
%

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses. 

Msub.grid = M.grid(rows,cols);
[Msub.nrows Msub.ncols] = size(Msub.grid);

if isfield(M,'x')
    Msub.x = M.x(cols);
    Msub.dx = abs(mean(diff(Msub.x)));
end

if isfield(M,'y')
    Msub.y = M.y(rows);
    Msub.dy = abs(mean(diff(Msub.y)));
end
