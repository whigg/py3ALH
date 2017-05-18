# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:01:12 2017

@author: nanda
"""

#Gome2datascripts = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/scripts'

import os
from read_AAI_FRESCO import readFrescoAAI
#sys.path.insert(0,Gome2datascripts)
import readGome2radIrr
import numpy as np


def createjobfiles(args):
      
    with open(os.path.join(args['outputfilepath'],args['filename']), 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('mkdir -p {tempfilepath}/{tempfilename} \n\n'.format(**args))

        f.write("python3 {inputcodedir}/{script} {arguments}\n\n"\
        .format(**args))

        f.write('rm -r {tempfilepath}/{tempfilename} \n'.format(**args))



#boundingbox = [60.0, 25.0, 52.5, 60.0]
boundingbox = [60.0, 29.0, 52.5, 45.0]

""" 10th August """
#year = '2010'
#month = '08'
#day = '10'
#
#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/10/S5P_OPER_AUX_MET_2D_20100810T000000_20100810T090000_20100810T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/10/S5P_TEST_L1B_IR_UVN_20100810T063917_20100810T082041_19756_02_010000_20160511T123520.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/10/S5P_TEST_L1B_RA_BD6_20100810T063917_20100810T082041_19756_02_010000_20160511T123520.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/10/G2A_OFFL_L2__FRESCO_20100810T063917_20100810T082041_19756_02_001002_20161115T020814.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/10/G2A_OFFL_L2__AER_AI_20100810T063917_20100810T082041_19756_02_001002_20161115T020754.nc'    
#time = '0639-0820'
#scanlines = (120,153)

#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/10/S5P_OPER_AUX_MET_2D_20100810T000000_20100810T090000_20100810T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/10/S5P_TEST_L1B_IR_UVN_20100810T082041_20100810T100159_19757_02_010000_20160511T124125.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/10/S5P_TEST_L1B_RA_BD6_20100810T082041_20100810T100159_19757_02_010000_20160511T124125.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/10/G2A_OFFL_L2__FRESCO_20100810T082041_20100810T100159_19757_02_001002_20161115T024856.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/10/G2A_OFFL_L2__AER_AI_20100810T082041_20100810T100159_19757_02_001002_20161115T024836.nc'    
#time = '0820-1001'
#scanlines = (120,153)


""" 9th August """
year = '2010'
month = '08'
day = '09'

#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/09/S5P_OPER_AUX_MET_2D_20100809T000000_20100809T090000_20100809T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/09/S5P_TEST_L1B_IR_UVN_20100809T070011_20100809T084135_19742_02_010000_20160511T110740.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/09/S5P_TEST_L1B_RA_BD6_20100809T070011_20100809T084135_19742_02_010000_20160511T110740.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/09/G2A_OFFL_L2__FRESCO_20100809T070011_20100809T084135_19742_02_001002_20161115T005324.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/09/G2A_OFFL_L2__AER_AI_20100809T070011_20100809T084135_19742_02_001002_20161115T005303.nc'    
#time = '0700-0841'
#scanlines = (120,153)

#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/09/S5P_OPER_AUX_MET_2D_20100809T000000_20100809T090000_20100809T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/09/S5P_TEST_L1B_IR_UVN_20100809T084135_20100809T102253_19743_02_010000_20160511T111358.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/09/S5P_TEST_L1B_RA_BD6_20100809T084135_20100809T102253_19743_02_010000_20160511T111358.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/09/G2A_OFFL_L2__FRESCO_20100809T084135_20100809T102253_19743_02_001002_20161115T010415.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/09/G2A_OFFL_L2__AER_AI_20100809T084135_20100809T102253_19743_02_001002_20161115T010354.nc'    
#time = '0841-1022'
#scanlines = (120,153)

""" 8th August """
year = '2010'
month = '08'
day = '08'

meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/08/S5P_OPER_AUX_MET_2D_20100808T000000_20100808T090000_20100808T000000.nc'
filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/08/S5P_TEST_L1B_IR_UVN_20100808T072105_20100808T090229_19728_02_010000_20160511T094025.nc'
filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/08/S5P_TEST_L1B_RA_BD6_20100808T072105_20100808T090229_19728_02_010000_20160511T094025.nc'
fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/08/G2A_OFFL_L2__FRESCO_20100808T072105_20100808T090229_19728_02_001002_20161114T221440.nc'    
aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/08/G2A_OFFL_L2__AER_AI_20100808T072105_20100808T090229_19728_02_001002_20161114T221419.nc'    
time = '0721-0902'
scanlines = (120,153)
filenameRad_band5 = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/08/S5P_TEST_L1B_RA_BD5_20100808T072105_20100808T090229_19728_02_010000_20160511T094025.nc'

""" 7th August """
#year = '2010'
#month = '08'
#day = '07'
#
#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/07/S5P_OPER_AUX_MET_2D_20100807T000000_20100807T090000_20100807T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/07/S5P_TEST_L1B_IR_UVN_20100807T074200_20100807T092318_19714_02_010000_20160511T081329.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/07/S5P_TEST_L1B_RA_BD6_20100807T074200_20100807T092318_19714_02_010000_20160511T081329.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/07/G2A_OFFL_L2__FRESCO_20100807T074200_20100807T092318_19714_02_001002_20161114T201451.nc'
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/07/G2A_OFFL_L2__AER_AI_20100807T074200_20100807T092318_19714_02_001002_20161114T201430.nc'    
#time            = '0742-0923'
#scanlines = (130,155)

""" 6th August """
#year = '2010'
#month = '08'
#day = '06'
#
#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/06/S5P_OPER_AUX_MET_2D_20100806T000000_20100806T090000_20100806T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/06/S5P_TEST_L1B_IR_UVN_20100806T080254_20100806T094412_19700_02_010000_20160511T064734.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/06/S5P_TEST_L1B_RA_BD6_20100806T080254_20100806T094412_19700_02_010000_20160511T064734.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/06/G2A_OFFL_L2__FRESCO_20100806T080254_20100806T094412_19700_02_001002_20161114T184905.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/06/G2A_OFFL_L2__AER_AI_20100806T080254_20100806T094412_19700_02_001002_20161114T184844.nc'    
#time            = '0802-0944'
#scanlines = (120,153)

""" 5th August """
#year = '2010'
#month = '08'
#day = '05'
#
#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/05/S5P_OPER_AUX_MET_2D_20100805T000000_20100805T090000_20100805T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/05/S5P_TEST_L1B_IR_UVN_20100805T064224_20100805T082348_19685_02_010000_20160511T051529.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/05/S5P_TEST_L1B_RA_BD6_20100805T064224_20100805T082348_19685_02_010000_20160511T051529.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/05/G2A_OFFL_L2__FRESCO_20100805T064224_20100805T082348_19685_02_001002_20161114T164459.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/05/G2A_OFFL_L2__AER_AI_20100805T064224_20100805T082348_19685_02_001002_20161114T164438.nc'    
#time            = '0642-0823'
#scanlines = (120,153)

""" 2nd August """
#year = '2010'
#month = '08'
#day = '02'
#
#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/02/S5P_OPER_AUX_MET_2D_20100802T000000_20100802T090000_20100802T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/02/S5P_TEST_L1B_IR_UVN_20100802T074507_20100802T092625_19643_02_010000_20160511T005757.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/02/S5P_TEST_L1B_RA_BD6_20100802T074507_20100802T092625_19643_02_010000_20160511T005757.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/02/G2A_OFFL_L2__FRESCO_20100802T074507_20100802T092625_19643_02_001002_20161114T104815.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/02/G2A_OFFL_L2__AER_AI_20100802T074507_20100802T092625_19643_02_001002_20161114T104755.nc'    
#time            = '0745-0926'
#scanlines = (120,153)

""" 1st August """
#year = '2010'
#month = '08'
#day = '01'
#
#meteofile       = '/net/pc160061/nobackup/users/tuinder/sat/projects/GRIB_for_AerosolCCI/2010/08/01/S5P_OPER_AUX_MET_2D_20100801T000000_20100801T090000_20100801T000000.nc'
#filenameIrr     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/01/S5P_TEST_L1B_IR_UVN_20100801T080602_20100801T094720_19629_02_010000_20160510T233152.nc'
#filenameRad     = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/L1b_Tropomi_format/2010/08/01/S5P_TEST_L1B_RA_BD6_20100801T080602_20100801T094720_19629_02_010000_20160510T233152.nc'
#fresco          = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__FRESCO/2010_20161116/08/01/G2A_OFFL_L2__FRESCO_20100801T080602_20100801T094720_19629_02_001002_20161114T090153.nc'    
#aai             = '/net/mos/data/obsrs/projects/o3msaf/experiments/AerosolCCI/metopa/gome2/G2A_OFFL_L2__AER_AI/2010_20161116/08/01/G2A_OFFL_L2__AER_AI_20100801T080602_20100801T094720_19629_02_001002_20161114T090132.nc'    
#time            = '0806-0947'
#scanlines = (120,153)

#configFile = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/S5P_OPER_CFG_AERLHF_00000000T000000_99999999T999999_20160715T101100.cfg'
configFile = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/S5P_OPER_CFG_AERLHF_00000000T000000_99999999T999999_20160715T101100.cfg'
tmp = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/tempdir'
band = 6
outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/tightenedAOT'
simfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/simfiles'
cfgo2b = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/S5P_OPER_CFG_AERLHF_00000000T000000_99999999T999999_20160715T101100__O2A-O2B.cfg'

FilterData = False


#pixels_for_processing = [(i,j) for i in range(scanlines[0],scanlines[1]+1) for j in range(0,24)]
#pixels_for_processing = [(150,18),(150,17),(146,13),(145,12),(147,12),(151,17),(139,19),(148,10),(143,18),(153,18)]
#pixels_for_processing = [(138,13)]
#pixels_for_processing = [(142,15),(140,16),(147,11),(142,16),(143,14),(143,17),(150,17),(138,18),(142,18),(139,16),(140,17),
#                         (147,13),(140,18),(141,18),(152,17),(141,17),(141,15),(141,16),(142,14),(149,6),(146,12),(142,13),
#                         (142,17),(145,13),(144,17),(139,18),(148,13),(146,13),(152,18),(138,17),(138,19),(147,10),(138,13),
#                         (138,16),(146,11),(143,15),(139,17),(145,11),(153,17),(143,13),(144,13),(151,18)] # prefit does not converge at all
                         
pixels_for_processing = [(137,15),(145,12),(153,16),(146,18),(147,12),(143,18),(152,16),(144,14),(145,18),(144,18),(153,18),
                         (150,18),(137,16),(151,17),(139,15),(143,16),(148,10),(145,14),(145,16),(145,7),(139,19)] #prefit runs but alh unqualified

                         
pf=1

for selPxl in pixels_for_processing:
    
    print(selPxl)
    
    i = selPxl[0]
    j = selPxl[1]

    if (os.path.exists(meteofile)) and (os.path.exists(fresco)) and (os.path.exists(filenameRad)) and (os.path.exists(filenameIrr)) and (os.path.exists(aai)):


        pixelSelectionSimFile = os.path.join(simfilepath,'PIXEL_SELECTION_SIMFILE.sim')
        lat,lon = readGome2radIrr.writeSimFile_GOME2(pixelSelectionSimFile,filenameRad,filenameIrr,i,j,band)[0:2]
        runALHprocess = False

        if FilterData:

            fresco_cld_fraction, AAI = readFrescoAAI(fresco,aai,i,j)
            
            if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                    and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ) \
                    and ( (fresco_cld_fraction<0.2) and (AAI>2.0) ):
                    
                runALHprocess = True
                    
        else:

            if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                    and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):
           
                runALHprocess=True

        if runALHprocess:
            
            simfilename    = 'GeneratedSimFile_GOME2-A__{0}_{1}_{2}_{3}_{4}_{5}.sim'.format(i,j,year,month,day,time)
                
            args = {'outputfilepath': '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/jobfiles',
                    'filename' : 'Job_{0}_{1}_{2}_{3}_{4}_{5}.sh'.format(i,j,year,month,day,time),
                    'tempfilepath' : tmp,
                    'tempfilename' : 'tmp_{0}_{1}_{2}_{3}_{4}_{5}'.format(i,j,year,month,day,time),
                    'inputcodedir' : '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/scripts/py3ALH',
                    'script' : 'runJob.py',
                    'arguments' : '-scl {0} -pxl {1} -y {2} -m {3} -d {4} -t {5} -f {6} -a {7} -s {8} -R {9} -I {10} -M {11} -C {12} -b {13} --outdir {14} -T {15}'.format(i,j,
                    year,month,day,time,fresco,aai,os.path.join(simfilepath,simfilename),filenameRad,filenameIrr, 
                    meteofile, configFile,band,outputfilepath,tmp)
                    }
            
            createjobfiles(args)

            
    else:
        
        print('ERROR: File not found')
        break
