function WriteSurferGrid(M,filename)

% WriteSurferGrid(M,limits,'filename') writes the values in matrix M.grid 
% with x and y coordinates M.x and M.y to the Surfer ASCII grid 
% filename.grd. 
% Note that filename must be a character string. 
%
% See Surfer documentation for details on Surfer ASCII grid file format.

% Copyright (C) 2004-2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

if nargin ~= 2, help(mfilename), return, end

% strip the '.grd' from filename if it's there; we'll add it back on later
[sizefnr sizefnc] = size(filename);
if strcmpi(filename(sizefnc-3:sizefnc),'.grd')
    filename = filename(1:sizefnc-4);    
end

% if filename.grd exists, delete it
if exist([filename '.grd'],'file')
    delete([filename '.grd'])
end

% write the header file
xmin = min(M.x);
xmax = max(M.x);
ymin = min(M.y);
ymax = max(M.y);

fid=fopen([filename '.grd'],'w');
fprintf(fid, 'DSAA\n%d %d\n%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n', M.ncols, M.nrows, ...
    xmin, xmax, ymin, ymax, min(M.grid(:)), max(M.grid(:)));

% convert NaNs to Surfer nodata value
M.grid(isnan(M.grid)) = 1.70141e38;

% Write the matrix values to the file

% show a progress bar
h = waitbar(0,'Writing Surfer grid...');

for i = M.nrows:-1:1
    fprintf(fid, ' %g',M.grid(i,:));
    fprintf(fid, '\n');
    f = (M.nrows-i+1)/M.nrows;
    if ~rem(round(f*100),10)
        waitbar(f)
    end
end

close(h);
fclose(fid);
