-------------------------------------------------------------
file: README.txt
-------------------------------------------------------------
-------------------------------------------------------------
Battelle Memorial Institute
Pacific Northwest Laboratory
-------------------------------------------------------------
-------------------------------------------------------------
Created January 20, 2016 by William A. Perkins
Last Change: 2016-03-08 10:08:46 d3g096
-------------------------------------------------------------

This folder contains the results of cell-by-cell statistical
computation for transient (hourly) 10-m MASS2 simulations of the
Hanford Reach.  For each MASS2 cell, the first four zero moments were
accumulated over two different time periods. From those, the 2nd, 3rd,
and 4th central moments were computed.  These were imported into a GIS
and interpolated onto a 5-m raster grid. Additional GIS raster layers
for coefficient of variation, skew, and kurtosis were computed from
the central moment layers.

mass2_10m_modern.csv

Contains statistics for each cell of the 10m MASS2 simulation from
1976 to 2013.  This period is the "modern era" during which all of the
large storage projects, particularly in Canada, were in place. The
fields are

  NWet - count of hourly time slices in which the cell was wet
  WetPercent - percent of time wet
  depth - water depth, m
  vmag - velocity magnitude, m/s
  shear - bottom shear stress, N/m2
  vstar - shear velocity (just a different representation of shear), m/s 

For each of the four latter fields, the first 4 central moments are
included and the standard deviation (square root of second moment).
It should be clear which is which from the field names.  

modern_gis/mass2_10m_modern_*.img

Raster GIS layers of computed statistics in the "extended" domain. 

modern/mass2_10m_modern_*.png

Plots of the modern_gis/mass2_10m_modern_*.img layers.

mass2_10m_premcnary.csv

Contains statistics for each cell of the 10m MASS2 simulation from
October 1917 through September 1953.  This is the period prior to the
construction of McNary Dam, which inundates the Columbia river in the
vicinity of the 300 Area. Fields are the same as
mass2_10m_modern.csv. 

premcnary_gis/mass2_10m_premcnary_*.img

Raster GIS layers of computed statistics in the "extended" domain. 

premcnary/mass2_10m_premcnary_*.png

Plots of the premcnary_gis/mass2_10m_premcnary_*.img layers.


