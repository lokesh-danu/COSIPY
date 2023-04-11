# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 09:44:46 2020

@author: Bharat Joshi
"""

from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication,QWidget, QVBoxLayout, QPushButton, QFileDialog , QLabel, QTextEdit, QMessageBox
import sys
from PyQt5.QtGui import QPixmap
from CreateStaticFileUI import *
import sys
import os
import xarray as xr
import numpy as np
from itertools import product
import richdem as rd
a =  0.25
b = []
def div(n):
    if n<0.0003:
        return
    else:
        j = n/2
        b.append(j)
        div(j)
        k = n/5
        b.append(k)
        div(k)
    
div(a)
b.sort()
c=list(set(b))
c.append(a)
c.sort()
for i in range(11):
    c.pop(0)

def signals(self):
    self.BDem.clicked.connect(self.browse1)
    self.BShape.clicked.connect(self.browse2)
    self.BOutput.clicked.connect(self.browse3)
    self.CreateStaticFile.clicked.connect(self.create_static_file)
    self.adjustVals.clicked.connect(self.edits)
    
    

def browse1(self):
    filename = QFileDialog.getOpenFileName()
    path1 = filename[0]
    self.TDem.setText(path1)
def browse2(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TShape.setText(path2)
def browse3(self):
    outputPath = str(QFileDialog.getExistingDirectory())
    self.TOutput.setText(outputPath)
def edits(self):
    try:
        LLon = float(self.TLLon.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Left Longitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        Ulat = float(self.TUlat.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Upper Latitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        RLon = float(self.TRLon.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Right Longitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        LLat = float(self.TLLat.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Lower Latitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        aggTemp = float(self.TAggregateDegree.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Aggregate Degree')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    aggTrue = c[min(range(len(c)), key = lambda i: abs(c[i]-aggTemp))]
    self.TAggregateDegree.setText(str(aggTrue))
    self.TLLon.setText(str(LLon-(aggTrue/2)))
    self.TUlat.setText(str(Ulat+(aggTrue/2)))
    self.TRLon.setText(str(RLon+(aggTrue/2)))
    self.TLLat.setText(str(LLat-(aggTrue/2)))                  

            

        

def create_static_file(self):
    static_folder = str(self.TOutput.text())+"/"
    tile = True
    aggregate = True
    
          
    
    ### input digital elevation model (DEM)
    dem_path_tif = str(self.TDem.text())
    if dem_path_tif=="":
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter DEM Path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    else:
        ### input shape of glacier or study area, e.g. from the Randolph glacier inventory
        shape_path = str(self.TShape.text())
        if shape_path=="":
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText('Enter Shape File Path')
            msg.setWindowTitle("Error")
            msg.exec_()
            return
        else:
            ### path were the static.nc file is saved
            output_path = static_folder + str(self.TOutputName.text()) + ".nc"
            if str(self.TOutput.text())=="":
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Error")
                msg.setInformativeText('Enter Output File Path')
                msg.setWindowTitle("Error")
                msg.exec_()
                return
            else:
                ### intermediate files, will be removed afterwards
                dem_path_tif_temp = static_folder + 'DEM_temp.tif'
                dem_path_tif_temp2 = static_folder + 'DEM_temp2.tif'
                dem_path = static_folder + 'dem.nc'
                aspect_path = static_folder + 'aspect.nc'
                mask_path = static_folder + 'mask.nc'
                slope_path = static_folder + 'slope.nc'
                
                ### If you do not want to shrink the DEM, comment out the following to three lines
                if tile:
                    ### to shrink the DEM use the following lat/lon corners
                    try:
                        longitude_upper_left = str(float(self.TLLon.text()))
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('Enter a Valid Left Longitude')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                    try:
                        latitude_upper_left = str(float(self.TUlat.text()))
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('Enter a Valid Upper Latitude')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                    try:
                        longitude_lower_right = str(float(self.TRLon.text()))
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('Enter a Valid Right Longitude')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                    try:
                        latitude_lower_right = str(float(self.TLLat.text()))
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('Enter a Valid Lower Latitude')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                    try:
                        os.system('gdal_translate -projwin ' + longitude_upper_left + ' ' + latitude_upper_left + ' ' +
                                  longitude_lower_right + ' ' + latitude_lower_right + ' ' +'"'+ dem_path_tif +'"'+  ' ' + '"'+ dem_path_tif_temp+'"')
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('GDal Error - Could not shrink the DEM. Check Input DEM file, Lat/Lon values and Output File Path')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                        
                    dem_path_tif = dem_path_tif_temp
                
                ### If you do not want to aggregate DEM, comment out the following to two lines
                if aggregate:
                    try:
                        ### to aggregate the DEM to a coarser spatial resolution
                        aggregate_degree = str(float(self.TAggregateDegree.text()))
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('Enter a Valid Aggregation Degree')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                    try:
                        os.system('gdalwarp -tr ' + aggregate_degree + ' ' + aggregate_degree + ' -r average ' + '"'+ dem_path_tif +'"'+  ' ' + '"'+ dem_path_tif_temp2+'"')
                    except:
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Critical)
                        msg.setText("Error")
                        msg.setInformativeText('GDal Error - Could not aggregate the DEM. Check Input DEM file, Output File path and aggregation degree')
                        msg.setWindowTitle("Error")
                        msg.exec_()
                        return
                    dem_path_tif = dem_path_tif_temp2
                
                ### convert DEM from tif to NetCDF
                try:
                    os.system('gdal_translate -of NETCDF ' + '"'+ dem_path_tif  +'"'+  ' ' +'"'+  dem_path+'"' )
                except:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Error")
                    msg.setInformativeText('GDal Error - Could not convert DEM from TIF to NETCDF - Check DEM File and its path')
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
                    
                
                ### calculate slope as NetCDF from DEM
                try:
                    os.system('gdaldem slope -of NETCDF ' + '"'+ dem_path +'"'+  ' ' + '"'+ slope_path +'"'+  ' -s 111120')
                except:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Error")
                    msg.setInformativeText('GDal Error - Could not calculate slope as NETCDF from DEM - Check DEM File, its path and output file path')
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
                    
                
                ### calculate aspect from DEM
                try:
                    demDS = rd.LoadGDAL(dem_path_tif)
                    aspect = np.flipud(rd.TerrainAttribute(demDS, attrib = 'aspect'))
                except:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Error")
                    msg.setInformativeText('RichDEM Error - Could not calculate Aspect from DEM - Check DEM File and its path')
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
                
                ### calculate mask as NetCDF with DEM and shapefile
                try:
                    os.system('gdalwarp -of NETCDF  --config GDALWARP_IGNORE_BAD_CUTLINE YES -cutline ' +'"'+  shape_path +'"'+  ' -cblend 0.5 ' + '"'+ dem_path_tif  +'"'+  ' ' + '"'+ mask_path+'"')
                except:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Error")
                    msg.setInformativeText('GDal Error - Could not calculate Mask as NETCDF with DEM and Shapefile - Check DEM File and its path, Shape File and its path,  and Output file path')
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
                
                ### open intermediate netcdf files
                try:
                    dem = xr.open_dataset(dem_path)
                    mask1 = xr.open_dataset(mask_path)
                    slope = xr.open_dataset(slope_path)
                except:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Error")
                    msg.setInformativeText('Xarray Error - Install Xarray Package')
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
                
                ### set NaNs in mask to -9999 and elevation within the shape to 1
                mask=mask1.Band1.values
                mask[np.isnan(mask)]=-9999
                mask[mask>0]=1
                print(mask)
                
                ## create output dataset
                ds = xr.Dataset()
                ds.coords['lon'] = dem.lon.values
                ds.lon.attrs['standard_name'] = 'lon'
                ds.lon.attrs['long_name'] = 'longitude'
                ds.lon.attrs['units'] = 'degrees_east'
                
                ds.coords['lat'] = dem.lat.values
                ds.lat.attrs['standard_name'] = 'lat'
                ds.lat.attrs['long_name'] = 'latitude'
                ds.lat.attrs['units'] = 'degrees_north'
                
                ### function to insert variables to dataset
                def insert_var(ds, var, name, units, long_name):
                    ds[name] = (('lat','lon'), var)
                    ds[name].attrs['units'] = units
                    ds[name].attrs['long_name'] = long_name
                    ds[name].attrs['_FillValue'] = -9999
                
                ### insert needed static variables
                insert_var(ds, dem.Band1.values,'HGT','meters','meter above sea level')
                insert_var(ds, aspect,'ASPECT','degrees','Aspect of slope')
                insert_var(ds, slope.Band1.values,'SLOPE','degrees','Terrain slope')
                insert_var(ds, mask,'MASK','boolean','Glacier mask')
                '''
                os.system('rm '+ '"'+ dem_path +'"'+  ' ' + '"'+ mask_path +'"'+  ' ' + '"'+ slope_path + '"'+ ' ' + '"'+ dem_path_tif_temp +'"'+  ' '+ '"'+ dem_path_tif_temp2+'"')
                '''
                dem_path = os.path.abspath(dem_path)
                mask_path = os.path.abspath(mask_path)
                slope_path = os.path.abspath(slope_path)
                dem_path_tif_temp = os.path.abspath(dem_path_tif_temp)
                dem_path_tif_temp2 = os.path.abspath(dem_path_tif_temp2)
                
                dem.close()
                mask1.close()
                slope.close()
                os.remove(dem_path)
                os.remove(mask_path)
                os.remove(slope_path)
                os.remove(dem_path_tif_temp)
                '''
                ### save combined static file, delete intermediate files and print number of glacier grid points
                def check_for_nan(ds,var=None):
                    for y,x in product(range(ds.dims['lat']),range(ds.dims['lon'])):
                        mask = ds.MASK.isel(lat=y, lon=x)
                        if mask==1:
                            if var is None:
                                if np.isnan(ds.isel(lat=y, lon=x).to_array()).any():
                                    msg = QMessageBox()
                                    msg.setIcon(QMessageBox.Critical)
                                    msg.setText("Error")
                                    msg.setInformativeText('ERROR!!!!!!!!!!! There are NaNs in the static fields')
                                    msg.setWindowTitle("Error")
                                    msg.exec_()
                                    return
                                    sys.exit()
                            else:
                                if np.isnan(ds[var].isel(lat=y, lon=x)).any():
                                    msg = QMessageBox()
                                    msg.setIcon(QMessageBox.Critical)
                                    msg.setText("Error")
                                    msg.setInformativeText('ERROR!!!!!!!!!!! There are NaNs in the static fields')
                                    msg.setWindowTitle("Error")
                                    msg.exec_()
                                    return
                                    sys.exit()
                check_for_nan(ds)
                '''
                try:
                    ds.to_netcdf(output_path)
                except:
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Error")
                    msg.setInformativeText('Could not save output as NETCDF file - Check Output Path')
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
                    
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Information)
                msg.setText("Done")
                msg.setInformativeText("Study area consists of "+str(np.nansum(mask[mask==1]))+" glacier points. Done")
                msg.setWindowTitle("Done")
                msg.exec_()
                return


        

     

        

Ui_MainWindow.signals = signals
Ui_MainWindow.browse1 = browse1
Ui_MainWindow.browse2 = browse2
Ui_MainWindow.browse3 = browse3
Ui_MainWindow.edits = edits
Ui_MainWindow.create_static_file = create_static_file

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    ui.signals()
    MainWindow.show()
    sys.exit(app.exec_())
    