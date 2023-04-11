# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import netCDF4 as nc
import numpy as np
import time
an = './Sato_static.nc'
fn = '../input/SatoInputFile.nc'
ds = nc.Dataset(an,mode='r')
print(ds.variables)
mb = ds.variables['lat'][:]

