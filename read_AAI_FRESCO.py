# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 13:35:43 2017

@author: nanda
"""

import tables

#fresco          = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/FRESCO_AAI_data/G2A_OFFL_L2__FRESCO_20100810T082041_20100810T100159_19757_02_001002_20161115T024856.nc'    
#aai             = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/FRESCO_AAI_data/G2A_OFFL_L2__AER_AI_20100810T082041_20100810T100159_19757_02_001002_20161115T024836.nc'    
#scl = 108
#pxl = 20

def readFrescoAAI(fresco,aai,scl,pxl):

    with tables.open_file(fresco) as f:
        fresco_cld_fraction = f.root.PRODUCT.cloud_fraction[0,scl,pxl]
        
    with tables.open_file(aai) as g:
        AAI = g.root.PRODUCT.aerosol_index_340_380[0,scl,pxl]
        
    return fresco_cld_fraction, AAI