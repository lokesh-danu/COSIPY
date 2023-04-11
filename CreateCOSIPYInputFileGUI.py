# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 11:03:20 2021

@author: Bharat Joshi
"""

from PyQt5 import QtGui
from PyQt5.QtWidgets import *
import sys
from CreateCOSIPYInputFileUI import *
import sys
import os
import xarray as xr
import numpy as np
from itertools import product, islice
import csv
import netCDF4 as nc
import math
import datetime as dt
import pandas as pd
from scipy.spatial.distance import cdist
from cosipy.modules.radCor import correctRadiation


def div(n,b):
    if n<0.0003:
        return
    else:
        j = n/2
        b.append(j)
        div(j,b)
        k = n/5
        b.append(k)
        div(k,b)
    

    

def signals(self):
    self.BStatic.clicked.connect(self.browse1)
    self.BTemp.clicked.connect(self.browse2)
    self.BPres.clicked.connect(self.browse3)
    self.BPrec.clicked.connect(self.browse4)
    self.BSnow.clicked.connect(self.browse5)
    self.BTcc.clicked.connect(self.browse6)
    self.BSrd.clicked.connect(self.browse7)
    self.B10u.clicked.connect(self.browse8)
    self.B10v.clicked.connect(self.browse9)
    self.B2d.clicked.connect(self.browse10)
    self.BOutput.clicked.connect(self.browse11)
    self.pushButton.clicked.connect(self.create_cosipy_input_file)
    
    
    

def browse1(self):
    filename = QFileDialog.getOpenFileName()
    path1 = filename[0]
    self.TStatic.setText(path1)
def browse2(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TTemp.setText(path2)
def browse3(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TPres.setText(path2)
def browse4(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TPrec.setText(path2)
def browse5(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TSnow.setText(path2)
def browse6(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TTcc.setText(path2)
def browse7(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.TSrd.setText(path2)
def browse8(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.T10u.setText(path2)
def browse9(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.T10v.setText(path2)
def browse10(self):
    filename = QFileDialog.getOpenFileName()
    path2 = filename[0]
    self.T2d.setText(path2)
def browse11(self):
    outputPath = QFileDialog.getSaveFileName()
    path = outputPath[0]
    self.TOutput.setText(path)
    
    




def create_cosipy_input_file(self):
    try:
        LonL = float(self.TLLon.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Left Longitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        LatU = float(self.TUlat.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Upper Latitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        LonR = float(self.TRLon.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Right Longitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        LatD = float(self.TLLat.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Lower Latitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        lapse_T = float(self.TLRTemp_2.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Temperature Lapse Rate')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        lapse_RH = float(self.TLRRh_2.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Relative Humidity Lapse Rate')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        lapse_RRR = float(self.TLRPrep_2.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Precipitation Lapse Rate')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        lapse_SNOWFALL = float(self.TLRSnow_2.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Snowfall Lapse Rate')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        z0 = float(self.TRoughness_2.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Roughness Length for Momentum')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        timezone_lon = float(self.TTimezone.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Timezone Longitude')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        zeni_thld = float(self.TZenith.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Zenith Threshold. Zenith Threshold is the maximum potential solar zenith angle during the whole year, specific for each location. If you do not know the exact value for your location, set value to 89.0')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    try:
        agg = float(self.TAggregateDegree.text())
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Aggregate Degree')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    if self.TOutput.text()=="":
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Output Path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
        
    a =  0.25
    b = []
    div(a,b)
    b.sort()
    c=list(set(b))
    c.append(a)
    c.sort()
    for i in range(11):
        c.pop(0)
    for i in c: 
        if(i != agg):
            agg = c[min(range(len(c)), key = lambda i: abs(c[i]-agg))]
            self.TAggregateDegree.setText(str(agg))

    
    lat = np.arange(start=LatD, stop=LatU+0.00001, step=agg).T
    lon = np.arange(start=LonL, stop=LonR+0.00001, step=agg) 
    date0 = self.dateFirst.dateTime().toPyDateTime()
    date1 = self.startDate.dateTime().toPyDateTime()
    date2 = self.endDate.dateTime().toPyDateTime()
    try:
        diff0 = date1 - date0
        diff1 = date2 - date0
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Dates & Time')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
        
    days, seconds = diff0.days, diff0.seconds
    hours0 = days * 24 + seconds // 3600
    days, seconds = diff1.days, diff1.seconds
    hours1 = days * 24 + seconds // 3600
    print("Total hours between "+str(date1)+" and "+str(date2)+" are "+str(hours1-hours0))
           
            
    
    
    
    
    
    try:
        fn = self.TStatic.text()
        ds = nc.Dataset(fn,mode='r')
        latSD = ds.variables['lat'][:]
        latS = len(latSD)
        lonSD = ds.variables['lon'][:]
        lonS = len(lonSD)
        Ele = ds.variables['HGT'][:]
        print("Loaded Elevation Data from Static File")
        
        
        Asp = ds.variables['ASPECT'][:]-180
        print("Loaded Aspect Data from Static File")
        
        
        Slo = ds.variables['SLOPE'][:]
        print("Loaded Slope Data from Static File")
        
        
        Mas = ds.variables['MASK'][:]
        print("Loaded Glacier Mask Data from Static File")
        MasVec = Mas.flatten() 
        MasSorted = np.sort(MasVec)
        SortedIndex = np.argsort(MasVec)
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Static File and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
        
    
    latN = len(np.arange(start=LatD, stop=LatU+0.00001, step=0.25))
    lonN = len(np.arange(start=LonL, stop=LonR+0.00001, step=0.25))
    # Load Temperature
    try:
        valsTemp = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.TTemp.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        valsTemp[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Temperature Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Temperature CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Surface Pressure
    try:
        valsPres = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.TPres.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        valsPres[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Surface Pressure Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Surface Pressure CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Precipitation 
    try:
        valsPrec = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.TPrec.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        valsPrec[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Precipitation Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Precipitation CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Snowfall
    try:
        valsSnow = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.TSnow.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        valsSnow[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Snowfall Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Snowfall CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Total Cloud Cover
    try:
        valsTCC = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.TTcc.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        valsTCC[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Cloud Cover Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Cloud Cover CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Incoming Solar Radiation Downward
    try:
        valsSRD = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.TSrd.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        valsSRD[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Incoming Solar Radiation Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Incoming Solar Radiation CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Wind 10mU
    try:
        vals10u = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.T10u.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        vals10u[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Wind Data U10")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Wind U10 CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Wind 10mV
    try:
        vals10v = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.T10v.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        vals10v[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Wind Data V10")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Wind V10 CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Load Dew Point
    try:
        vals2d = np.zeros((hours1-hours0,latN,lonN),dtype=float)
        with open(self.T2d.text()) as fd:
            for i in range(hours0,hours1):
                d=[[]]
                for row in islice(csv.reader(fd), 1, 46):
                    d.append(row)
                d.remove([])
                for f in range(latN-1,-1,-1):
                    for g in range(lonN):
                        vals2d[i-hours0,f,g]=d[(8-f)*5+g][2]
                del d
        
        print("Loaded Dew Point Data")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter a Valid Dew Point CSV and check its path')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
    # Calculate Relative Humidity
    valsRh = np.zeros((hours1-hours0,latN,lonN),dtype=float)
    for a in range(latN):
        for b in range(lonN):
            for c in range(hours0,hours1):
                valsRh[c-hours0,a,b] = 100*(np.power(math.e,((17.625*vals2d[c-hours0,a,b])/(243.04+vals2d[c-hours0,a,b])))/np.power(math.e,((17.625*valsTemp[c-hours0,a,b])/(243.04+valsTemp[c-hours0,a,b]))))
                if valsRh[c-hours0,a,b]>100:
                    valsRh[c-hours0,a,b]=100
                elif valsRh[c-hours0,a,b]<0.1:
                    valsRh[c-hours0,a,b]=0.1
                
    print("Calculated Relative Humidity Data")
    



    
    # Interpolate 
    try:
        ITemp = np.zeros((hours1-hours0,len(lat),len(lon)))
        IPres = np.zeros((hours1-hours0,len(lat),len(lon)))
        IPrec = np.zeros((hours1-hours0,len(lat),len(lon)))
        ISnow = np.zeros((hours1-hours0,len(lat),len(lon)))
        ITcc = np.zeros((hours1-hours0,len(lat),len(lon)))
        ISrd = np.zeros((hours1-hours0,len(lat),len(lon)))
        Iu2 = np.zeros((hours1-hours0,len(lat),len(lon)))
        IRh = np.zeros((hours1-hours0,len(lat),len(lon)))
        trV=[[]]
        tcV=[[]]
        i=0
        while(MasSorted[i]>0):
            f = SortedIndex[i]
            tr,tc = divmod(f,lon.size)
            trV.append(tr)
            tcV.append(tc)
            i+=1
        trV.pop(0)
        tcV.pop(0)
        
        for q in range(min(trV),-1,-1):
            if round(lat[q],10)%0.25==0:
                LatMin = q
                break
        for q in range(max(trV),lat.size):
            if round(lat[q],10)%0.25==0:
                LatMax = q
                break
        for q in range(min(tcV),-1,-1):
            if round(lon[q],10)%0.25==0:
                LonMin = q
                break
        for q in range(max(tcV),lon.size):
            if round(lon[q],10)%0.25==0:
                LonMax = q
                break
        
        
        LatERA=[[]]
        LonERA=[[]]
        for q in range(LatMin,LatMax+1):
            if round(lat[q],10)%0.25==0:
                LatERA.append(q)
        for q in range(LonMin,LonMax+1):
            if round(lon[q],10)%0.25==0:
                LonERA.append(q)
        LatERA.pop(0)
        LonERA.pop(0)
        
        
        ERAPoints = [[]]
        for i in LatERA:
            for j in LonERA:
                ERAPoints.append([i,j])
        ERAPoints.pop(0)
        ERAPoints = np.array(ERAPoints)
        
        if len(np.arange(start=LonL, stop=LonL+0.25+0.00001, step=agg))%2==0:
            N = len(np.arange(start=LonL, stop=LonL+0.25+0.00001, step=agg))/2
            M = N*2-1
        else:
            N = math.floor(len(np.arange(start=LonL, stop=LonL+0.25+0.00001, step=agg))/2)+1
            M = N*2-2
        print("Interpolating Data")
        for t in range(hours0,hours1):
            for i in range(len(tcV)):
                SumTemp = 0
                SumPres = 0
                SumPrec = 0
                SumSnow = 0
                SumTcc = 0
                SumSrd = 0
                Sumu2 = 0
                SumRh = 0
                count = 0
                
                for j in LatERA:
                    for k in LonERA:
                        
                        SumTemp+= valsTemp[t-hours0,int(j/M),int(k/M)]+(Ele[trV[i],tcV[i]]-Ele[j,k])*lapse_T
                        SumPres+= (valsPres[t-hours0,int(j/M),int(k/M)]/np.power((1-(0.0065*Ele[j,k])/(288.15)),5.255))*np.power((1-(0.0065*Ele[trV[i],tcV[i]])/(288.15)),5.255)/100
                        SumPrec+= np.maximum(valsPrec[t-hours0,int(j/M),int(k/M)]+(Ele[trV[i],tcV[i]]-Ele[j,k])*lapse_RRR,0)*1000
                        SumSnow+= np.maximum(valsSnow[t-hours0,int(j/M),int(k/M)]+(Ele[trV[i],tcV[i]]-Ele[j,k])*lapse_SNOWFALL,0)
                        Rh = valsRh[t-hours0,int(j/M),int(k/M)]+(Ele[trV[i],tcV[i]]-Ele[j,k])*lapse_RH
                        if Rh>100:
                            SumRh+= 100
                        elif Rh<0.1:
                            SumRh+= 0.1
                        else:
                            SumRh+= Rh
                        
                        SumTcc+= valsTCC[t-hours0,int(j/M),int(k/M)]
                        
                        Sumu2+= (np.sqrt(vals10v[t-hours0,int(j/M),int(k/M)]**2+vals10u[t-hours0,int(j/M),int(k/M)]**2))*(np.log(2/z0)/np.log(10/z0))
                        
                        '''
                        SumSrdnp.maximum(0.0, correctRadiation(lat[j],lon[k], timezone_lon, doy, hour, Slo[j,k], Asp[j,k], valsSRD[t,j,k], zeni_thld))/3600
                        SumSrd+= valsSRD[t,int(j/M),int(k/M)]/3600
                        '''
                        updatedDate = date1 + dt.timedelta(hours=t)
                        day = updatedDate.day
                        hour = updatedDate.hour
                        SumSrd+= np.maximum(0.0, correctRadiation(lat[trV[i]],lon[tcV[i]], timezone_lon, day, hour, Slo[trV[i],tcV[i]], Asp[trV[i],tcV[i]], valsSRD[t-hours0,int(j/M),int(k/M)], zeni_thld))/3600
                        
                        count+=1
                
                        
                
                
                '''
                #Taking nearest point value instead of Averaging
                Dist = np.argmin(cdist(ERAPoints, np.array([(trV[i],tcV[i])]), 'euclidean'))
                ITcc[t,trV[i],tcV[i]] = valsTCC[t,int(ERAPoints[Dist,0]/M),int(ERAPoints[Dist,1]/M)]
                Iu2[t,trV[i],tcV[i]] = (np.sqrt(vals10v[t,int(ERAPoints[Dist,0]/M),int(ERAPoints[Dist,1]/M)]**2+vals10u[t,int(ERAPoints[Dist,0]/M),int(ERAPoints[Dist,1]/M)]**2))*(np.log(2/z0)/np.log(10/z0))
                ISrd[t,trV[i],tcV[i]] = np.maximum(0.0, correctRadiation(lat[trV[i]],lon[tcV[i]], timezone_lon, day, hour, Slo[j,k], Asp[j,k], valsSRD[t,int(ERAPoints[Dist,0]/M),int(ERAPoints[Dist,1]/M)], zeni_thld))/3600
                '''
                
                ITemp[t-hours0,trV[i],tcV[i]] = SumTemp/count
                IPres[t-hours0,trV[i],tcV[i]] = SumPres/count
                IPrec[t-hours0,trV[i],tcV[i]] = SumPrec/count
                ISnow[t-hours0,trV[i],tcV[i]] = SumSnow/count
                IRh[t-hours0,trV[i],tcV[i]] = SumRh/count
                ISrd[t-hours0,trV[i],tcV[i]] = SumSrd/count
                Iu2[t-hours0,trV[i],tcV[i]] = Sumu2/count
                ITcc[t-hours0,trV[i],tcV[i]] = SumTcc/count
            if (t-hours0)%24==0:
                print("Data interpolated till "+str(date0+dt.timedelta(hours=t)))
                '''
                         
        #Old Code
        
        if len(np.arange(start=LonL, stop=LonL+0.25+0.00001, step=agg))%2==0:
            N = len(np.arange(start=LonL, stop=LonL+0.25+0.00001, step=agg))/2
            M = N*2-1
        else:
            N = math.floor(len(np.arange(start=LonL, stop=LonL+0.25+0.00001, step=agg))/2)+1
            M = N*2-2
        
        print("Interpolating Data")
        
        for t in range(hours):
            if t%1000==0:
                print("Data interpolated till "+str(date1+dt.timedelta(hours=t+1)))
            for i in range(latN):
                for j in range(lonN):
                    if i==0 and j==0:
                        for x in range(N):
                            for y in range(N):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif i==0 and j==lonN-1:
                        for x in range(N):
                            for y in range(0,-N,-1):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif i==8 and j==0:
                        for x in range(0,-N,-1):
                            for y in range(N):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif i==8 and j==lonN-1:    
                        for x in range(0,-N,-1):
                            for y in range(0,-N,-1):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif i==0:
                        for x in range(N):
                            for y in range(N-1,-N,-1):      
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif i==8:
                        for x in range(0,-N,-1):
                            for y in range(N-1,-N,-1):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif j==0:
                        for x in range(N-1,-N,-1):
                            for y in range(N):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif j==lonN-1:
                        for x in range(N-1,-N,-1):
                            for y in range(0,-N,-1):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
                    elif i>0 and i<8 and j>0 and j<lonN-1:
                        for x in range(N-1,-N,-1):
                            for y in range(N-1,-N,-1):
                                ITemp[t,M*i+x,M*j+y] = valsTemp[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_T
                                IPres[t,M*i+x,M*j+y] = (valsPres[t,i,j]/np.power((1-(0.0065*Ele[M*i,M*j])/(288.15)),5.255))*np.power((1-(0.0065*Ele[M*i+x,M*j+y])/(288.15)),5.255)/100
                                IPrec[t,M*i+x,M*j+y] = np.maximum(valsPrec[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RRR,0)*1000
                                ISnow[t,M*i+x,M*j+y] = np.maximum(valsSnow[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_SNOWFALL,0)
                                ITcc[t,M*i+x,M*j+y] = valsTCC[t,i,j]
                                ISrd[t,M*i+x,M*j+y] = valsSRD[t,i,j]/3600
                                Iu2[t,M*i+x,M*j+y] = (np.sqrt(vals10v[t,i,j]**2+vals10u[t,i,j]**2))*(np.log(2/z0)/np.log(10/z0))
                                IRh[t,M*i+x,M*j+y] = valsRh[t,i,j]+(Ele[M*i+x,M*j+y] - Ele[M*i,M*j])*lapse_RH
                                if IRh[t,M*i+x,M*j+y]>100:
                                    IRh[t,M*i+x,M*j+y]=100
                                elif IRh[t,M*i+x,M*j+y]<0.1:
                                    IRh[t,M*i+x,M*j+y]=0.1
        '''           
                        
        print("Data interpolated till "+str(date0+dt.timedelta(hours=t)))
        
        print("Interpolation Completed") 
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Error in Interpolation. Check Data Formats and Paths')
        msg.setWindowTitle("Error")
        msg.exec_()
        return                    
    
    
    
    
    try:
        # Writing Data to netCDF4 File
        print("Writing Data to netCDF4 File")
        
        
        # Create Dataframe
        time = pd.date_range(start=date1,freq='1H',periods=hours1-hours0)
        latF = np.arange(start=LatD, stop=LatU+0.00001, step=agg)
        lonF = np.arange(start=LonL, stop=LonR+0.00001, step=agg) 
        # Create Dataset
        ds = xr.Dataset({
            'HGT': xr.DataArray(
                        data   = Ele,   # enter data here
                        dims   = ['lat', 'lon'],
                        coords = {'lat':latF, 'lon':lonF},
                        attrs  = {
                            '_FillValue' : -9999.0,
                            'units' : 'm',
                            'long_name' : 'Elevation'
                            }
                        ),
            'ASPECT': xr.DataArray(
                        data   = Asp,   # enter data here
                        dims   = ['lat', 'lon'],
                        coords = {'lat':latF, 'lon':lonF},
                        attrs  = {
                            '_FillValue' : -9999.0,
                            'units' : 'degrees',
                            'long_name' : 'Aspect of slope' 
                            }
                        ),
            'SLOPE': xr.DataArray(
                        data   = Slo,   # enter data here
                        dims   = ['lat', 'lon'],
                        coords = {'lat':latF, 'lon':lonF},
                        attrs  = {
                            '_FillValue' : -9999.0,
                            'units' : 'degrees',
                            'long_name' : 'Terrain Slope' 
                            }
                        ),
            'MASK': xr.DataArray(
                        data   = Mas,   # enter data here
                        dims   = ['lat', 'lon'],
                        coords = {'lat':latF, 'lon':lonF},
                        attrs  = {
                            '_FillValue' : -9999.0,
                            'units' : 'boolean',
                            'long_name' : 'Glacier mask' 
                            }
                        ),
            'T2': xr.DataArray(
                        data   = ITemp,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : 'K',
                            'long_name' : 'Temperature at 2 m' 
                            }
                        ),
            'RH2': xr.DataArray(
                        data   = IRh,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : '%',
                            'long_name' : 'Relative humidity at 2 m' 
                            }
                        ),
            'U2': xr.DataArray(
                        data   = Iu2,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : 'm s⁻¹',
                            'long_name' : 'Wind velocity at 2 m' 
                            }
                        ),
            'G': xr.DataArray(
                        data   = ISrd,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : 'W m⁻²',
                            'long_name' : 'Incoming shortwave radiation' 
                            }
                        ),
            'PRES': xr.DataArray(
                        data   = IPres,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : 'hPa',
                            'long_name' : 'Atmospheric Pressure' 
                            }
                        ),
            'RRR': xr.DataArray(
                        data   = IPrec,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : 'mm',
                            'long_name' : 'Total precipitation' 
                            }
                        ),
            'SNOWFALL': xr.DataArray(
                        data   = ISnow,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : 'm',
                            'long_name' : 'Snowfall' 
                            }
                        ),
            'N': xr.DataArray(
                        data   = ITcc,   # enter data here
                        dims   = ['time','lat', 'lon'],
                        coords = {'time':time,'lat':latF, 'lon':lonF},
                        attrs  = {
                            'units' : '%',
                            'long_name' : 'Cloud cover fraction' 
                            }
                        )},
                attrs = {'Description': 'COSIPY Input File'}
            )
        ds.lon.attrs['standard_name'] = 'lon'
        ds.lon.attrs['long_name'] = 'longitude'
        ds.lon.attrs['units'] = 'degrees_east'
        ds.lat.attrs['standard_name'] = 'lat'
        ds.lat.attrs['long_name'] = 'latitude'
        ds.lat.attrs['units'] = 'degrees_north'
        try:
            TOutputTemp=os.path.abspath(self.TOutput.text())
            os.remove(TOutputTemp)
            print("Overwriting file.")
            
        except:
            print("No File to Delete")
            
        ds.to_netcdf(self.TOutput.text())
        print("Wrote Data Successfully")
        
    except:
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setInformativeText('Enter in saving output file')
        msg.setWindowTitle("Error")
        msg.exec_()
        return
    
   
    
    
    
    
    
    
    
    
    
    
Ui_MainWindow.signals = signals
Ui_MainWindow.browse1 = browse1
Ui_MainWindow.browse2 = browse2
Ui_MainWindow.browse3 = browse3
Ui_MainWindow.browse4 = browse4
Ui_MainWindow.browse5 = browse5
Ui_MainWindow.browse6 = browse6
Ui_MainWindow.browse7 = browse7
Ui_MainWindow.browse8 = browse8
Ui_MainWindow.browse9 = browse9
Ui_MainWindow.browse10 = browse10
Ui_MainWindow.browse11 = browse11
Ui_MainWindow.create_cosipy_input_file = create_cosipy_input_file

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    ui.signals()
    MainWindow.show()
    sys.exit(app.exec_())
