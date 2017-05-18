# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:22:21 2017

@author: nanda
"""

import netCDF4

import numpy as np
import datetime
from scipy.interpolate import interp1d

#filename = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/03/S5P_OPER_AUX_MET_2D_20100803T000000_20100803T090000_20100803T000000.nc'
#inp      = (59.99,40.00, 1.851883349E7)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx,array[idx]

def readMeteoS5Pformat(filename,inp):
    """
    examples:
    
    inp      = (59.99,40.00, 3)
    filename = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/03/S5P_OPER_AUX_MET_2D_20100803T000000_20100803T090000_20100803T000000.nc'

    """
#    inp      = (59.99,40.00, 1.851883349E7)
    
    data = netCDF4.Dataset(filename)
        
    lat     = data.variables['latitude'][:]
    lon     = data.variables['longitude'][:]
    a       = data.variables['hyam'][:]
    b       = data.variables['hybm'][:]
    
    sp              = data.variables['sp'][:]
#    sp_addoffset    = data.variables['sp'].add_offset
#    sp_mulScaling   = data.variables['sp'].scale_factor
    
    t              = data.variables['t'][:]
#    t_addoffset    = data.variables['t'].add_offset
#    t_mulScaling   = data.variables['t'].scale_factor
    
    time    = data.variables['time'][2:]*3600
    
    
    inpTime = datetime.datetime.utcfromtimestamp(inp[2])
    timeToInterpolate = inpTime.hour * 3600. + inpTime.minute * 60. + inpTime.second
    
    print('meteo data is interpolated at time {0}:{1}:{2}'.format(inpTime.hour,inpTime.minute,inpTime.second))
    
    nlat = find_nearest(lat,inp[0])
    nlon = find_nearest(lon,inp[1])
    
    
    """
    Half pressure is calculated as P = a + b * SurfacePressure, i.e., the boundaries of the layers.
    source: https://rda.ucar.edu/datasets/ds115.4/docs/levels.hybrid.html
    
    """
    
    #nSP = ( sp_addoffset + sp[2:,nlat[0],nlon[0]]   * sp_mulScaling )/ 100.
    nSP = sp[2:,nlat[0],nlon[0]]/ 100.
    
    nSP_TIME = np.interp(timeToInterpolate,time,nSP)
    
    tempProf = np.vstack(t[2:,:,nlat[0],nlon[0]])
    interpTempProf =  interp1d(time,tempProf,axis=0)(timeToInterpolate)
    
    #nT  =   t_addoffset  +  interpTempProf *  t_mulScaling
    nT = interpTempProf
    
    nP  = a + b * nSP_TIME
    
    return nlat[1], nlon[1], nP[::-1], nT[::-1], nSP_TIME    # array is flipped to conform to disamar configFile
