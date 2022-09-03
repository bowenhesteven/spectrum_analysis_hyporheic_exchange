TopoTools for Matlab v1.0
2010/07/07

These Matlab functions perform various tasks related to I/O and analysis of gridded topographic data. Please acknowledge the use of this software in any publications with a citation along the lines of the following:

Perron, J.T. (2010) TopoTools for Matlab (Version 1.0) [Computer program]. Available at http://eaps.mit.edu/faculty/perron (Accessed DD Month Year).

All code is intended for use with Matlab. See comments in individual m-files for usage instructions. For help with any function, type 'help [function name]' at the Matlab prompt.

All code copyright (C) 2004-2010 Taylor Perron (perron@mit.edu), except where otherwise indicated. These programs are free software: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation. You should have received a copy of the GNU General Public License along with these programs. If not, see http://www.gnu.org/licenses.


Included files:

Documentation/notes:
  aaREADME.txt (this file)

License:
  license.txt

Matlab m-files:
  Arc2Surf.m			Convert ESRI ASCII grid to Surfer ASCII grid
  ColorShade.m (1)		Calculate and plot colored shaded relief
  D8.m					Calculate drainage area using steepest descent (D8)
  Dinf.m				Calculate drainage area using Tarboton's D-infinity method
  ReadArcGrid.m			Read ESRI ASCII grid into Matlab
  ReadSurferGrid.m		Read Surfer ASCII grid into Matlab
  Shade.m (2)			Calculate and plot shaded relief
  Surf2Arc.m			Convert Surfer ASCII grid to ESRI ASCII grid 
  WriteArcGrid.m		Write ESRI ASCII grid from Matlab
  WriteSurferGrid.m		Write Surfer ASCII grid from Matlab

C mex-files for Matlab:
  mexD8.c
  mexDinf.c

Notes
(1) Based on imshade.m by Kelsey Jordahl. See m-file header for license information. Obtained from the Matlab Central File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/25691-imshade
(2) Based on hillshade.m by Felix Hebeler. Obtained from the Matlab Central File Exchange:
http://www.mathworks.com/matlabcentral/fileexchange/14863-hillshade
