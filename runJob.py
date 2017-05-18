# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 14:20:27 2017

@author: nanda
"""

import traceback, argparse
disamar = '/usr/people/nanda/disamar/Disamar.exe'


#Imports:
import sys, os, tables
import rt_cfg, rt_run
from obtainTempPresProfile_fromMeteoData import readMeteoS5Pformat
import readGome2radIrr
from extractsurfacealbedo import extractsurfalbedo
from read_AAI_FRESCO import readFrescoAAI
import numpy as np
import copy

def readPF(filename):
    
    with tables.open_file(filename) as f:
        if f.root._f_getattr('solution_has_converged') == 'true':
            ret_aot = f.root.parameters.aerosol_tau.read()[2]
            ret_aot_precision = f.root.parameters.precision_bias_aerosol_tau.read()[1]
        else:
            ret_aot = -1
            ret_aot_precision = -1
    
    return ret_aot,ret_aot_precision

def runPreFit(args):
    
    layerThickness = 5.0
    
    inarg          = { 'outpath_loc'    : args.outloc, 
                       'scanline'       : args.scanline,
                       'pixel'          : args.pixel,
                       'band'           : args.band,
                       'year'           : args.year, 'month' : args.month, 'day': args.day, 'time': args.time, 
                       'frescoFile'     : args.frescoFile,
                       'aaiFile'        : args.aaiFile,
                       'simfile'        : args.simfile,
                       'rad'            : args.rad, 'irr': args.irr, 'meteofile': args.meteofile,
                       'cfg'            : args.cfg, 'tmp':args.tmp }
                   
    outputfilepathLER = inarg['outpath_loc']

    maincfg             = rt_cfg.RT_configuration(inarg['cfg'])

    maincfg['INSTRUMENT', 'wavelength_range', 'wavelength_start'].setvalue(755.0)
    maincfg['INSTRUMENT', 'wavelength_range', 'wavelength_end'].setvalue(756.0)    
    maincfg['INSTRUMENT', 'wavelength_range', 'wavelength_step'].setvalue(0.2)
    
    maincfg['INSTRUMENT', 'calibrationErrorRefl', 'wavelength'].setvalue([755.0,756.0])
    maincfg['INSTRUMENT', 'SNR_irradiance', 'wavelSNR'].setvalue(755.5)
    maincfg['INSTRUMENT', 'SNR_radiance', 'wavelSNR'].setvalue(755.5)
    maincfg['MUL_OFFSET', 'simulation', 'wavelengths'].setvalue(755.5)
    maincfg['MUL_OFFSET', 'retrieval', 'wavelengths'].setvalue(755.5)
    maincfg['STRAY_LIGHT', 'simulation', 'wavelengths'].setvalue(755.5)
    maincfg['STRAY_LIGHT', 'retrieval', 'wavelengths'].setvalue(755.5)    
    maincfg['GENERAL', 'specifyFitting', 'fitIntervalDP'].setvalue(0)    

     # create simfile and read lat/lon/time data
    lat,lon,dT,geometry,latB,lonB = readGome2radIrr.writeSimFile_GOME2(inarg['simfile'],inarg['rad'],inarg['irr'],inarg['scanline'],inarg['pixel'],inarg['band'])
            
     # pressure temperature profile - SOME ISSUES HERE FOR SURE
    PTdata = readMeteoS5Pformat( inarg['meteofile'],(lat,lon,dT) )
    pressureTemperatureRetr = []
    for iPT,PT in enumerate(PTdata[2]):
        if iPT == 0:
            if PT < PTdata[-1]:
                pressureTemperatureRetr.append([PTdata[-1], PTdata[3][iPT], 1.0])
        else:
            
            pressureTemperatureRetr.append([ PTdata[2][iPT], PTdata[3][iPT], 1.0 ])       
            
    salb758 = extractsurfalbedo('Mode_LER',8,758.,lat,lon,[0.,1.])
    salb772 = extractsurfalbedo('Mode_LER',8,772.,lat,lon,[0.,1.])

    maincfg['SURFACE','brdf','useBRDFSim'].setvalue(0)
    maincfg['SURFACE','brdf','useBRDFRetr'].setvalue(0)
    
    #SURFACE HERE - POSSIBLE TO CHANGE?
    maincfg['SURFACE','wavelDependentRetr','wavelSurfAlbedo'].setvalue([755.,756.])
    maincfg['SURFACE','wavelDependentRetr','surfAlbedo'].setvalue([salb758[-1],salb758[-1]])       

    maincfg['GEOMETRY', 'geometry', 'solar_zenith_angle_retr'].setvalue(geometry[0])
    maincfg['GEOMETRY', 'geometry', 'solar_azimuth_angle_retr'].setvalue(geometry[2])
    maincfg['GEOMETRY', 'geometry', 'instrument_nadir_angle_retr'].setvalue(geometry[1])
    maincfg['GEOMETRY', 'geometry', 'instrument_azimuth_angle_retr'].setvalue(geometry[3])
    maincfg['GEOMETRY', 'geometry', 'solar_zenith_angle_sim'].setvalue(geometry[0])
    maincfg['GEOMETRY', 'geometry', 'solar_azimuth_angle_sim'].setvalue(geometry[2])
    maincfg['GEOMETRY', 'geometry', 'instrument_nadir_angle_sim'].setvalue(geometry[1])
    maincfg['GEOMETRY', 'geometry', 'instrument_azimuth_angle_sim'].setvalue(geometry[3])

    maincfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','topPressureRetr'].setvalue([pressureTemperatureRetr[0][0]-200.+layerThickness, pressureTemperatureRetr[0][0]-200.-layerThickness , 0.3])
 
    # here are the keys were most of the issues may take place:        
    maincfg['PRESSURE_TEMPERATURE', 'PT_retr', 'PT'].setvalue(pressureTemperatureRetr)

    O2VMR = maincfg['O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].values()
    O2vmr_new = []
    O2vmr_new.append([PTdata[-1],209460.0, 10.])
    for ivmr in O2VMR:    
        if ivmr[0] < PTdata[-1]:
            O2vmr_new.append(ivmr)
            
    maincfg['O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].setvalue(O2vmr_new)
    
    O2O2VMR = maincfg['O2-O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].values()
    O2O2vmr_new = []
    O2O2vmr_new.append([PTdata[-1],209460.0, 10.])
    for ivmr in O2O2VMR:            
        if ivmr[0] < PTdata[-1]:
            O2O2vmr_new.append(ivmr)
            
    maincfg['O2-O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].setvalue(O2O2vmr_new)
    
    maincfg['SURFACE','pressure','surfPressureRetr'].setvalue(PTdata[-1])

    maincfg['GENERAL','external_data', 'year'].setvalue(inarg['year'])
    maincfg['GENERAL','external_data', 'month'].setvalue(inarg['month'])
    maincfg['GENERAL','external_data', 'day'].setvalue(inarg['day'])    
    maincfg['GENERAL','external_data', 'latitude'].setvalue(lat)
    maincfg['GENERAL','external_data', 'longitude'].setvalue(lon)            
    maincfg['GENERAL','external_data', 'cornerLatitudes'].setvalue('_'.join([str(idx_latB) for idx_latB in latB]))
    maincfg['GENERAL','external_data', 'cornerLongitudes'].setvalue('_'.join([str(idx_lonB) for idx_lonB in lonB]))
    maincfg['GENERAL','external_data', 'AtrackNumber'].setvalue(inarg['scanline'])
    maincfg['GENERAL','external_data', 'XtrackNumber'].setvalue(inarg['pixel'])
    maincfg['GENERAL','external_data', 'originalIrradianceFileName'].setvalue(inarg['irr'])
    maincfg['GENERAL','external_data', 'originalRadianceVISFileName'].setvalue(inarg['rad'])
    
    print(maincfg)

    
    cfg_pf__1 = copy.deepcopy(maincfg)
    cfg_pf__2 = copy.deepcopy(maincfg)
    del maincfg
    pf__1 = 1.0
    
    """ Double Prefit filter """
                               
    try:
        
        cfg_pf__1['AEROSOL','HGscatteringRetr','opticalThickness'].setvalue([2,pf__1,1.0])
        cfg_pf__1['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltRetr'].setvalue([6,int(np.ceil(1.5*pf__1)),16])
                               
        tmpbase_pf__1 = os.path.join(inarg['tmp'],'PF/lo')
        baseoutfilename_pf__1 = '{0}/{1}_{2}_{3}__{4}__{5}-{6}__PF__1.h5'.format(outputfilepathLER,inarg['year'],inarg['month'],inarg['day'],inarg['time'],inarg['scanline'],inarg['pixel'])
        runner_pf__1 = rt_run.rt_run(cfg=cfg_pf__1, disamar=disamar, output=baseoutfilename_pf__1,
                               spectrum =inarg['simfile'], tempbase = tmpbase_pf__1 )

        # run 1st prefit
        runner_pf__1()
        
        aot_pf__1,aot_precision_pf__1 = readPF(baseoutfilename_pf__1)

        if aot_pf__1 < 1.0:
            
            pf__2 = aot_pf__1/2.0
            
        else:
            
            if aot_pf__1 < 10.0:
                pf__2 = aot_pf__1 + 2.0
            else:
                pf__2 = -1
        

        if pf__2 > 0.0:
            
            cfg_pf__2['AEROSOL','HGscatteringRetr','opticalThickness'].setvalue([2,pf__2,1.0])
            cfg_pf__2['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltRetr'].setvalue([6,int(np.ceil(1.5*pf__2)),16])        
            tmpbase_runner_pf__2 = os.path.join(inarg['tmp'],'PF/hi')
            baseoutfilename_pf__2 = '{0}/{1}_{2}_{3}__{4}__{5}-{6}__PF__2.h5'.format(outputfilepathLER,inarg['year'],inarg['month'],inarg['day'],inarg['time'],inarg['scanline'],inarg['pixel'])
            runner_pf__2 = rt_run.rt_run(cfg=cfg_pf__2, disamar=disamar, output=baseoutfilename_pf__2,
                                   spectrum =inarg['simfile'], tempbase = tmpbase_runner_pf__2)
    
            # run 2nd prefit                          
            runner_pf__2()
            
            aot_pf__2,aot_precision_pf__2 = readPF(baseoutfilename_pf__2)
                
            if (aot_pf__2 > 0.0) and (aot_pf__2 > 0.0) and (np.abs(aot_pf__1 - aot_pf__2) < 0.1*np.max([aot_pf__1,aot_pf__2])):
                
                fitALH = True
                PF = np.mean([aot_pf__1,aot_pf__2])
                
            else:
                
                fitALH = False
                print('prefit Failed')
                PF = -1.0            
        
        else:
            
            fitALH = False
            print('prefit Failed')
            PF = -1.0 

    except Exception as exception:
        
        print('Error while prefitting')
        print(exception)
        fitALH = False
        PF = -1
        
    return fitALH,PF


def runALH(args,PF,donePF):
    
    
    layerThickness = 5.0 #hPa
    
    inarg          = { 'outpath_loc'    : args.outloc, 
                       'scanline'       : args.scanline,
                       'pixel'          : args.pixel,
                       'band'           : args.band,
                       'year'           : args.year, 'month' : args.month, 'day': args.day, 'time': args.time, 
                       'frescoFile'     : args.frescoFile,
                       'aaiFile'        : args.aaiFile,
                       'simfile'        : args.simfile,
                       'rad'            : args.rad, 'irr': args.irr, 'meteofile': args.meteofile,
                       'cfg'            : args.cfg, 'tmp':args.tmp }
    
    cfg             = rt_cfg.RT_configuration(inarg['cfg'])

    fresco_cld_fraction, AAI = readFrescoAAI(inarg['frescoFile'],inarg['aaiFile'],inarg['scanline'],inarg['pixel'])
  
    outputfilepathLER = inarg['outpath_loc']
    baseoutfilename = '{0}/{1}_{2}_{3}__{4}__{5}-{6}__ALH.h5'.format(outputfilepathLER,inarg['year'],inarg['month'],inarg['day'],inarg['time'],inarg['scanline'],inarg['pixel'])
    
    # create simfile and read lat/lon/time data
    lat,lon,dT,geometry,latB,lonB = readGome2radIrr.writeSimFile_GOME2(inarg['simfile'],inarg['rad'],inarg['irr'],inarg['scanline'],inarg['pixel'],inarg['band'])
        
    # pressure temperature profile - SOME ISSUES HERE FOR SURE
    PTdata = readMeteoS5Pformat( inarg['meteofile'],(lat,lon,dT) )
    pressureTemperatureRetr = []
    for iPT,PT in enumerate(PTdata[2]):
        if iPT == 0:
            if PT < PTdata[-1]:
                pressureTemperatureRetr.append([PTdata[-1], PTdata[3][iPT], 1.0])
        else:
            
            pressureTemperatureRetr.append([ PTdata[2][iPT], PTdata[3][iPT], 1.0 ])       
            
    salb758 = extractsurfalbedo('Mode_LER',8,758.,lat,lon,[0.,1.])
    salb772 = extractsurfalbedo('Mode_LER',8,772.,lat,lon,[0.,1.])

    if donePF:
        cfg['AEROSOL','HGscatteringRetr','opticalThickness'].setvalue([2,PF,0.4])   
#        cfg['AEROSOL','HGscatteringRetr','opticalThickness'].setvalue([2,PF,0.000001])      
        cfg['RADIATIVE_TRANSFER','numDivPointsAlt','numDivPointsAltRetr'].setvalue([6,int(np.ceil(1.5*PF)),16])

    # values for russian fires 2010 derived from: http://www.atmos-meas-tech.net/5/557/2012/amt-5-557-2012.pdf : Table 2
    cfg['AEROSOL','HGscatteringRetr','singleScatteringAlbedo'].setvalue([2,0.95])
    cfg['AEROSOL','HGscatteringRetr','HGparameter_g'].setvalue([2,0.7])
    
    cfg['SURFACE','brdf','useBRDFSim'].setvalue(0)
    cfg['SURFACE','brdf','useBRDFRetr'].setvalue(0)
    cfg['SURFACE','wavelDependentRetr','wavelSurfAlbedo'].setvalue([758.,772.])
    cfg['SURFACE','wavelDependentRetr','surfAlbedo'].setvalue([salb758[-1],salb772[-1]])        

    cfg['GEOMETRY', 'geometry', 'solar_zenith_angle_retr'].setvalue(geometry[0])
    cfg['GEOMETRY', 'geometry', 'solar_azimuth_angle_retr'].setvalue(geometry[2])
    cfg['GEOMETRY', 'geometry', 'instrument_nadir_angle_retr'].setvalue(geometry[1])
    cfg['GEOMETRY', 'geometry', 'instrument_azimuth_angle_retr'].setvalue(geometry[3])
    cfg['GEOMETRY', 'geometry', 'solar_zenith_angle_sim'].setvalue(geometry[0])
    cfg['GEOMETRY', 'geometry', 'solar_azimuth_angle_sim'].setvalue(geometry[1])
    cfg['GEOMETRY', 'geometry', 'instrument_nadir_angle_sim'].setvalue(geometry[2])
    cfg['GEOMETRY', 'geometry', 'instrument_azimuth_angle_sim'].setvalue(geometry[3])

 
    # here are the keys were most of the issues may take place:        
    cfg['PRESSURE_TEMPERATURE', 'PT_retr', 'PT'].setvalue(pressureTemperatureRetr)

    O2VMR = cfg['O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].values()
    O2vmr_new = []
    O2vmr_new.append([PTdata[-1],209460.0, 10.])
    for ivmr in O2VMR:    
        if ivmr[0] < PTdata[-1]:
            O2vmr_new.append(ivmr)
            
    cfg['O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].setvalue(O2vmr_new)
    
    O2O2VMR = cfg['O2-O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].values()
    O2O2vmr_new = []
    O2O2vmr_new.append([PTdata[-1],209460.0, 10.])
    for ivmr in O2O2VMR:            
        if ivmr[0] < PTdata[-1]:
            O2O2vmr_new.append(ivmr)
            
    cfg['O2-O2', 'profile', 'P_vmr_ppmv_error_percent_retr'].setvalue(O2O2vmr_new)
    
    cfg['SURFACE','pressure','surfPressureRetr'].setvalue(PTdata[-1])

#     if mid pressure is fitted
    cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','topPressureRetr'].setvalue([pressureTemperatureRetr[0][0]-200.+layerThickness, pressureTemperatureRetr[0][0]-200.-layerThickness , 0.3])

    # if top pressure is fitted and bottom is fixed on the surface
#    cfg['GENERAL','specifyFitting','fitIntervalDP'].setvalue(0)
#    cfg['GENERAL','specifyFitting','fitIntervalTop'].setvalue(1)
#    cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','topPressureRetr'].setvalue([PTdata[-1]-5., pressureTemperatureRetr[0][0]-200.-25. , 0.3])

#     if top and bottom pressures are fitted separately
#    cfg['GENERAL','specifyFitting','fitIntervalDP'].setvalue(0)
#    cfg['GENERAL','specifyFitting','fitIntervalTop'].setvalue(1)
#    cfg['GENERAL','specifyFitting','fitIntervalBot'].setvalue(1)
#    cfg['ATMOSPHERIC_INTERVALS','interval_top_pressures','topPressureRetr'].setvalue([PTdata[-1]-200.+layerThickness, pressureTemperatureRetr[0][0]-200.-layerThickness , 0.3])


    cfg['GENERAL','external_data', 'year'].setvalue(inarg['year'])
    cfg['GENERAL','external_data', 'month'].setvalue(inarg['month'])
    cfg['GENERAL','external_data', 'day'].setvalue(inarg['day'])
    cfg['GENERAL','external_data', 'latitude'].setvalue(lat)
    cfg['GENERAL','external_data', 'longitude'].setvalue(lon)            
    cfg['GENERAL','external_data', 'cornerLatitudes'].setvalue('_'.join([str(idx_latB) for idx_latB in latB]))
    cfg['GENERAL','external_data', 'cornerLongitudes'].setvalue('_'.join([str(idx_lonB) for idx_lonB in lonB]))
    cfg['GENERAL','external_data', 'AtrackNumber'].setvalue(inarg['scanline'])
    cfg['GENERAL','external_data', 'XtrackNumber'].setvalue(inarg['pixel'])
    cfg['GENERAL','external_data', 'originalIrradianceFileName'].setvalue(inarg['irr'])
    cfg['GENERAL','external_data', 'originalRadianceVISFileName'].setvalue(inarg['rad'])
        
    # Retrieve
    runner = rt_run.rt_run(cfg=cfg, disamar=disamar, output=baseoutfilename,
                               spectrum =inarg['simfile'], tempbase = inarg['tmp'])
    runner()

    
def main():
    
    parser = argparse.ArgumentParser(description='Parse input for ALH retrieval from GOME-2 spectra')
   
    parser.add_argument('--outdir', '-o', metavar='DIR', dest='outloc', type=str, default='/net/pc150230/nobackup/users/nanda/Projects/tempOutput',
                        help='Put output files in DIR.')
                        
    parser.add_argument('--tempdir', '-T', metavar='DIR', dest='tmp', type=str, default = '/net/pc150230/nobackup/users/nanda/Projects/tempDir',
                        help='Put temporary files in DIR.')
                        
    parser.add_argument('--scanline', '-scl', metavar='N', dest='scanline', type=int, default = 0,
                        help='Enter scanline to process')
    
    parser.add_argument('--pixel', '-pxl', metavar='N', dest='pixel', type=int, default = 0,
                        help='Enter pixel to process')

    parser.add_argument('--band', '-b', metavar='N', dest='band', type=int, default = 6,
                        help='Enter TROPOMI band to process')

    parser.add_argument('--year', '-y', metavar='N', dest='year', type=str, default = '2010',
                        help='Enter year of data')

    parser.add_argument('--month', '-m', metavar='N', dest='month', type=str, default = '08',
                        help='Enter month of data')

    parser.add_argument('--day', '-d', metavar='N', dest='day', type=str, default = '10',
                        help='Enter day of data')

    parser.add_argument('--time', '-t', metavar='N', dest='time', type=str, default = '0820-1015',
                        help='Enter time of data')

    parser.add_argument('--fresco', '-f', metavar='DIR', dest='frescoFile', type=str, default=  '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/FRESCO_AAI_data/G2A_OFFL_L2__FRESCO_20100810T063917_20100810T082041_19756_02_001002_20161115T020814.nc',
                        help='FRESCO file directory.')

    parser.add_argument('--aai', '-a', metavar='DIR', dest='aaiFile', type=str, default=  '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/FRESCO_AAI_data/G2A_OFFL_L2__AER_AI_20100810T063917_20100810T082041_19756_02_001002_20161115T020754.nc',
                        help='AAI file directory.')
                        
    parser.add_argument('--sim', '-s', metavar='DIR', dest='simfile', type=str, default=  '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/GeneratedSimFile_Gome2.sim',
                        help='sim file directory.')    

    parser.add_argument('--rad', '-R', metavar='DIR', dest='rad', type=str, default=  '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/10/S5P_TEST_L1B_RA_BD6_20100810T063917_20100810T082041_19756_02_010000_20160511T123520.nc',
                        help='rad file directory.')
                        
    parser.add_argument('--irr', '-I', metavar='DIR', dest='irr', type=str, default=  '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/10/S5P_TEST_L1B_IR_UVN_20100810T063917_20100810T082041_19756_02_010000_20160511T123520.nc',
                        help='irr file directory.')  

    parser.add_argument('--meteo', '-M', metavar='DIR', dest='meteofile', type=str, default=  '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/10/S5P_OPER_AUX_MET_2D_20100810T000000_20100810T090000_20100810T000000.nc',
                        help='meteo file directory.') 

    parser.add_argument('--cfg', '-C', metavar='DIR', dest='cfg', type=str, default=  '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/S5P_OPER_CFG_AERLHF_00000000T000000_99999999T999999_20160715T101100.cfg',
                        help='cfg file directory.')

    args = parser.parse_args()
    
    try:
        fitALH,PF = runPreFit(args)
        if fitALH:
            runALH(args,PF,True)
        else:
            print("Pixel did not qualify Prefit filtering")
    except:
        print("failed to run the retrieval.")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()    
    
