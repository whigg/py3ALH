# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 11:48:41 2017

@author: nanda
"""

import tables
from pyhdf.SD import SD, SDC
from pyproj import Proj
import scipy.spatial as sp
import numpy as np
from osgeo import gdal

#modisDATA = '/nobackup/users/nanda/Projects/O2A_O2B_O2O2_ALH/MCD43A1_A2002017_h27v06_005_2007121074726.hdf'
#modisData_RussianFires = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/MODIS_MCD4A1/MCD43A1.A2010224.h20v03.006.2016079125358.hdf'

def find_index(data,coord):
    
    dist=sp.distance.cdist(data,coord)
    
    return dist.argmin()

def readMODIS_brdf(filename,band):
    
    """
    example:
    filename = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/MODIS_MCD4A1/MCD43A1.A2010224.h20v03.006.2016079125358.hdf'
    band = 2
    
    """
    
    if ( band < 1 ) or ( band > 7) :
        print('ERROR: band unavailable or does not exist')
        return
    
    hdf_file = gdal.Open(filename)
    subDatasets = hdf_file.GetSubDatasets()
    dataset = gdal.Open(subDatasets[0][0]) # this one is NDVI
    geotransform = dataset.GetGeoTransform()
    
    p_modis_grid = Proj('+proj=sinu +R=6371007.181 +nadgrids=@null +wktext')
    #x, y = p_modis_grid(lon, lat)
    # or the inverse, from x, y to lon, lat
#    x = 2223901.039333 + (2400*geotransform[1])
#    y = 6671703.118000 + (2400*geotransform[-1])
#    lon, lat = p_modis_grid(x,y, inverse=True)
#    print 'upper left lat/lon = {}'.format(p_modis_grid(geotransform[1],geotransform[5], inverse=True))
#    print 'bottom right lat/lon = {}'.format(p_modis_grid(lon,lat, inverse=True))
   
    f =   SD(filename, SDC.READ)
    brdf_b2  = f.select('BRDF_Albedo_Parameters_Band{0}'.format(band)).get()
    brdf_P = []
    
    for iy,y in enumerate(brdf_b2):
        for ix,x in enumerate(y):
            
            X = geotransform[0] + (ix+1)*geotransform[1]
            Y = geotransform[3] + (iy+1)*geotransform[5]
    
            lon, lat = p_modis_grid(X,Y, inverse=True)
    
            
            brdf_P.append([lat,lon,x[0]/1000.,x[1]/1000.,x[2]/1000.])
    
    return np.vstack(brdf_P)

def getBRDFparam(filename,band,nlat,nlon):
    
    """
    example:
    
    slat = 60.
    slon = 40.
    filename = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/MODIS_MCD4A1/MCD43A1.A2010224.h20v03.006.2016079125358.hdf'
    band = 1
    
    """

    data = readMODIS_brdf(filename,band)
        
    return data[find_index(data[:,0:2],np.array([[nlat,nlon]]))]
    
    
def getBRDFparams_allBands(filename,nlat,nlon):
    
    """
    example:
    
    slat = 60.
    slon = 40.
    filename = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/MODIS_MCD4A1/MCD43A1.A2010224.h20v03.006.2016079125358.hdf'
    band = 1
    
    """

    bands = [1,2,3,4]
    brdf_param_mband = []
    for band in bands:
        
        brdf_param_mband.append(getBRDFparam(filename,band,nlat,nlon))
    
    return brdf_param_mband

