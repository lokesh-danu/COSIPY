#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 22:13:42 2020

@author: bharat
"""

import sys
import os
import xarray as xr
import numpy as np
from itertools import product
import richdem as rd

static_folder = '../../data/static/'

tile = True
aggregate = True

### input digital elevation model (DEM)
dem_path_tif = static_folder + 'DEM/n30_e090_3arc_v2.tif'
### input shape of glacier or study area, e.g. from the Randolph glacier inventory
shape_path = static_folder + 'Shapefiles/Zhadang_RGI6.shp'
### path were the static.nc file is saved
output_path = static_folder + 'Zhadang_static.nc'

### to shrink the DEM use the following lat/lon corners
longitude_upper_left = '90.62'
latitude_upper_left = '30.48'
longitude_lower_right = '90.66'
latitude_lower_right = '30.46'

### to aggregate the DEM to a coarser spatial resolution
aggregate_degree = '0.003'

### intermediate files, will be removed afterwards
dem_path_tif_temp = static_folder + 'DEM_temp.tif'
dem_path_tif_temp2 = static_folder + 'DEM_temp2.tif'
dem_path = static_folder + 'dem.nc'
aspect_path = static_folder + 'aspect.nc'
mask_path = static_folder + 'mask.nc'
slope_path = static_folder + 'slope.nc'

### If you do not want to shrink the DEM, comment out the following to three lines
if tile:
    os.system('gdal_translate -projwin ' + longitude_upper_left + ' ' + latitude_upper_left + ' ' +
          longitude_lower_right + ' ' + latitude_lower_right + ' ' + dem_path_tif + ' ' + dem_path_tif_temp)
    dem_path_tif = dem_path_tif_temp
    
if aggregate:
    os.system('gdalwarp -tr ' + aggregate_degree + ' ' + aggregate_degree + ' -r average ' + dem_path_tif + ' ' + dem_path_tif_temp2)
    dem_path_tif = dem_path_tif_temp2    

os.system('gdalwarp -of NETCDF  --config GDALWARP_IGNORE_BAD_CUTLINE YES -cutline ' + shape_path + ' ' + dem_path_tif  + ' ' + mask_path)  
os.system('gdal_translate -of NETCDF ' + dem_path_tif  + ' ' + dem_path)