# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:25:23 2017

@author: nanda
"""

import tables
import numpy as np
import datetime


def readGOME2radIrrData(filenameRad,filenameIrr,scanline,pixel,band):
    
    """
    examples:
    
    filenameIrr = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_IR_UVN_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'
    filenameRad = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_RA_BD6_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'

    scanline = 15
    pixel    = 24
    band     = 6
    
    """
    
    with tables.open_file(filenameIrr) as f:

        wavelength = f.root._f_get_child('BAND{0}_IRRADIANCE'.format(int(band))).STANDARD_MODE.INSTRUMENT.calibrated_wavelength.read()[0,int(pixel),:]
        irradiance = f.root._f_get_child('BAND{0}_IRRADIANCE'.format(int(band))).STANDARD_MODE.OBSERVATIONS.irradiance.read()[0,0,int(pixel),:]
        signalToNoiseIrr= f.root._f_get_child('BAND{0}_IRRADIANCE'.format(int(band))).STANDARD_MODE.OBSERVATIONS.irradiance_noise.read()[0,0,int(pixel),:]
        
        wvl_irr = wavelength[wavelength[:]<1e35]
        irr = irradiance[irradiance[:]<1e35] * 6.02214179*10E23
        snr_irr = 10**(signalToNoiseIrr[irradiance[:]<1e35]/10.)
    
    with tables.open_file(filenameRad) as g:
        wavelength = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.INSTRUMENT.nominal_wavelength.read()[0,int(pixel),:]
        radiance = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.OBSERVATIONS.radiance.read()[0,int(scanline),int(pixel),:]
        signalToNoiseRad= g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.OBSERVATIONS.radiance_noise.read()[0,int(scanline),int(pixel),:]

        latB = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.latitude_bounds.read()[0,int(scanline),int(pixel),:]
        lonB = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.longitude_bounds.read()[0,int(scanline),int(pixel),:]
        
        sza = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.solar_zenith_angle.read()[0,int(scanline),int(pixel)]
        saa = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.solar_azimuth_angle.read()[0,int(scanline),int(pixel)]
        vza = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.viewing_zenith_angle.read()[0,int(scanline),int(pixel)]
        vaa = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.viewing_azimuth_angle.read()[0,int(scanline),int(pixel)]

        lat = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.latitude.read()[0,int(scanline),int(pixel)]
        lon = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.GEODATA.longitude.read()[0,int(scanline),int(pixel)]
        
        dTime = g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.OBSERVATIONS.time.read()[:] 
        dTime = dTime + g.root._f_get_child('BAND{0}_RADIANCE'.format(int(band))).STANDARD_MODE.OBSERVATIONS.delta_time.read()[0,int(scanline)]/1000.

        geodata = [sza,vza,saa,vaa,lat,lon]
        
        wvl_rad = wavelength[wavelength[:]<1e35]
        rad = radiance[radiance[:]<1e35] * 6.02214179*10E23
        snr_rad = 10**(signalToNoiseRad[radiance[:]<1e35]/10.)
        
        
    return wvl_irr, irr, snr_irr, wvl_rad, rad, snr_rad, geodata, dTime, latB, lonB

def writeSimFile_GOME2(simfilename,filenameRad,filenameIrr,scanline,pixel,band):

    """
    examples:
    
    filenameIrr = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_IR_UVN_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'
    filenameRad = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_RA_BD6_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'

    scanline = 15
    pixel    = 24
    band     = 6
    
    simfilename = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/GeneratedSimFile_Gome2.sim'
    
    """
    
#    print 'GOME-2 orbit irradiance/radiance data for ALH-summary-data scanline {0} and pixel {1}'.format(scanline,pixel)    

    wvl_irr, irr, snr_irr, wvl_rad, rad, snr_rad, geo, dTime, latB, lonB = readGOME2radIrrData(filenameRad,filenameIrr,scanline,pixel,band)
    dtime_timestamp = dTime
    dTime = datetime.datetime.utcfromtimestamp(dTime)
    
#    print 'GOME-2 orbit irradiance/radiance data is provided for lat {0} and lon {1}'.format(geo[-2],geo[-1])    
#    print 'GOME-2 orbit irradiance/radiance data is provided at time {0}:{1}:{2}'.format(dTime.hour,dTime.minute,dTime.second)    

#    print 'GOME-2 orbit data is derived for angles [deg]:\n'
#    print 'sza = {0}, vza = {1}, saa={2}, vaa = {3} '.format(geo[0],geo[1],geo[2],geo[3])    
    geometry = [geo[0],geo[1],geo[2],geo[3]]
    
    with open(simfilename,'w') as simFile:
        simFile.write('start_channel_irr \n')
        for ix, x in enumerate(irr):
            simFile.write('irr \t {0} \t {1} \t {2} \n'.format(wvl_irr[ix],snr_irr[ix],irr[ix]))
        simFile.write('end_channel_irr \nstart_channel_rad\n')
        for ix, x in enumerate(rad):
            simFile.write('rad \t {0} \t {1} \t {2} \n'.format(wvl_rad[ix],snr_rad[ix],rad[ix]))    
        simFile.write('end_channel_rad \nend_file')
    
    return geo[-2],geo[-1], dtime_timestamp, geometry, latB, lonB

def writeSimFile_GOME2__o2ao2b(simfilename,filenameRad,filenameRad2,filenameIrr,scanline,pixel,band):

    """
    examples:
    
    filenameIrr = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_IR_UVN_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'
    filenameRad = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_RA_BD6_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'
    filenameRad2= '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/GOME2_orbits/S5P_TEST_L1B_RA_BD6_20100803T072413_20100803T090531_19657_02_010000_20160511T022334.nc'

    scanline = 15
    pixel    = 24
    band     = 6 (or 5)
    
    simfilename = '/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/GeneratedSimFile_Gome2.sim'
    
    """
    
    print('GOME-2 orbit irradiance/radiance data for ALH-summary-data scanline {0} and pixel {1}'.format(scanline,pixel)    )

    wvl_irr1, irr1, snr_irr1, wvl_rad1, rad1, snr_rad1, geo, dTime, latB, lonB = readGOME2radIrrData(filenameRad2,filenameIrr,scanline,pixel,5)
    wvl_irr2, irr2, snr_irr2, wvl_rad2, rad2, snr_rad2, geo, dTime, latB, lonB = readGOME2radIrrData(filenameRad,filenameIrr,scanline,pixel,6)
    
    dtime_timestamp = dTime
    dTime = datetime.datetime.utcfromtimestamp(dTime)
    
    print('GOME-2 orbit irradiance/radiance data is provided for lat {0} and lon {1}'.format(geo[-2],geo[-1])    )
    print('GOME-2 orbit irradiance/radiance data is provided at time {0}:{1}:{2}'.format(dTime.hour,dTime.minute,dTime.second)    )

    print ('GOME-2 orbit data is derived for angles [deg]:\n')
    print ('sza = {0}, vza = {1}, saa={2}, vaa = {3} '.format(geo[0],geo[1],geo[2],geo[3])    )
    geometry = [geo[0],geo[1],geo[2],geo[3]]
    
    with open(simfilename,'w') as simFile:
        
        simFile.write('start_channel_irr \n')
        
        for ix, x in enumerate(irr1):
            
            simFile.write('irr \t {0} \t {1} \t {2} \n'.format(wvl_irr1[ix],snr_irr1[ix],irr1[ix]))
            
        simFile.write('end_channel_irr \n')
        
        simFile.write('start_channel_rad\n')
        
        for ix, x in enumerate(rad1):
            
            simFile.write('rad \t {0} \t {1} \t {2} \n'.format(wvl_rad1[ix],snr_rad1[ix],rad1[ix]))  
            
        simFile.write('end_channel_rad \n')
    
        simFile.write('start_channel_irr \n')
        
        for ix, x in enumerate(irr2):
            
            simFile.write('irr \t {0} \t {1} \t {2} \n'.format(wvl_irr2[ix],snr_irr2[ix],irr2[ix]))
            
        simFile.write('end_channel_irr \n')
        simFile.write('start_channel_rad\n')

        for ix, x in enumerate(rad2):
            
            simFile.write('rad \t {0} \t {1} \t {2} \n'.format(wvl_rad2[ix],snr_rad2[ix],rad2[ix]))  
            
        simFile.write('end_channel_rad \n')
        simFile.write('end_file')
        
    return geo[-2],geo[-1], dtime_timestamp, geometry, latB, lonB