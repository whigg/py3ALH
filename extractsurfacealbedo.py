# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 13:35:12 2015

@author: Swadhin Nanda, KNMI
"""

""" 
Paths 

Please set path to the location of the folder containing the GOME-2 or SCIAMACHY file in the variable 'filepath'
Please also set the name of the file in the variable 'filename'
"""

filepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data'
filename = 'GOME2_MetOp-A_MSC_surface_LER_product.hdf5'


""" imports """
import tables
import numpy as np

""" Basic function/definition to find a value in an array """
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]


def extractsurfalbedo(typedata,month,wavelength,latitude,longitude,acceptableflags):
    """
    Example

    Supply 'typedata' as either: 'Mode_LER' or 'Minimum_LER'
    Supply 'wavelength' as wavelength for the surface albedo you are looking for. The program will find the closest wavelength to it in the data
    Supply 'latitude' and 'longitude' as the coordinates you are looking for. The program will find the closest coordinates to it in the data.
    Supply 'month' as the month you want the data for. eg. 12
    Supply 'acceptableflags' as the list of flags acceptable to you, eg. [0.,1.,2.]
    
    """    
    
    
    """ assign hdf5 file to variable """
    f = tables.open_file( '{0}/{1}'.format(filepath,filename)  )
    
    """ find index of the wavelength, latitude, longitude that is the closest to the user-specified input"""
    waveindex, val   = find_nearest(f.root.Wavelength.read(),wavelength)
    latindex, valat  = find_nearest(f.root.Latitude.read(),latitude)
    lonindex, valon  = find_nearest(f.root.Longitude.read(),longitude)
    
    data = f.root._f_get_child(typedata).read()
    
    if f.root.Flag.read()[latindex,lonindex,month-1] in acceptableflags:
        surfalbedo = data[latindex,lonindex,waveindex,month-1]
    
    return valat, valon, surfalbedo
    
def main():
    """
    Example script to run the surface albedo extraction function
    
    """        
    typedata = 'Mode_LER'
    month = 8
    wavelength = 758.
    latitude = 55.
    longitude = 45.
    acceptableflags = [0.,1.]
    
    print(extractsurfalbedo(typedata,month,wavelength,latitude,longitude,acceptableflags))


if __name__ == "__main__":
    main()




