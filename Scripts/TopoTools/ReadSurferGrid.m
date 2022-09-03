function M = ReadSurferGrid(filename)

% [M dimensions] = ReadSurferGrid('myfile') reads the Surfer ASCII grid
% myfile.grd into a matrix of grid values with associated spatial reference 
% information. The grid is stored as a struct array with elements: 
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
% If the input character array filename ends in '.grd', ReadSurferGrid 
% strips the '.grd' extension.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

if (nargin ~= 1), help(mfilename), return, end

[sizefnr sizefnc] = size(filename);
if strcmpi(filename(sizefnc-3:sizefnc),'.grd')
    filename = filename(1:sizefnc-4);    
end

if ~exist([filename '.grd'],'file')
    error(['File ' filename '.grd does not exist.']);
end

% open the file and check the first line. Should read "DSAA". If it doesn't,
% tell the user it isn't a valid file & bail out.

fid = fopen([filename '.grd'],'r');
if ~strcmpi(fscanf(fid, ' %s',1),'DSAA')
    fclose(fid);
    error(['File ' filename '.grd is not a Surfer ASCII grid file.']);
end

% read in the grid description
dim = fscanf(fid, ' %g', 8);
fclose(fid);

M.ncols = dim(1);
M.nrows = dim(2);
M.dx = abs(dim(4)-dim(3))/(dim(1)-1);
M.dy = abs(dim(6)-dim(5))/(dim(2)-1);
M.xllcorner = dim(3)-0.5*M.dx;
M.yllcorner = dim(5)-0.5*M.dy;
M.x = M.xllcorner + M.dx * (0.5 + (0:M.ncols-1));
M.y = flipud(M.yllcorner + M.dy * (0.5 + (0:M.nrows-1)'));

% read in the grid values: open the file again and load the grid into the 
% variable M (skipping 5 header lines)
M.grid = zeros(M.nrows,M.ncols);
datafile = fopen([filename '.grd'],'r');
fscanf(datafile, '%*s', 9); % throw out the header in the Surfer grid file

% show a progress bar
h = waitbar(0,'Reading Surfer grid...');

for i = M.nrows:-1:1
    temp = fscanf(datafile, '%f', M.ncols);
    M.grid(i,1:M.ncols) = temp';
    f = (M.nrows-i+1)/M.nrows;
    if ~rem(round(f*100),10)
        waitbar(f)
    end
end
fclose(datafile);
clear temp;
close(h);

% Replace Surfer NODATA values (>=1.70141e38) with NaNs
M.grid(M.grid >= 1.70141e38) = NaN;
