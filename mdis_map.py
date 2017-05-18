#!/usr/bin/env python

from pyhdf.SD import SD, SDC 
from scipy.misc import bytescale
from mpl_toolkits.basemap import Basemap, cm 
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib import colors
from matplotlib.collections import PolyCollection

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import pprint
import os, collections
import tables

# References:
# https://earthdata.nasa.gov/sites/default/files/field/document/MODIS_True_Color.pdf
# http://www.idlcoyote.com/ip_tips/brightmodis.html


def read(retrievalhdf):
    
    with tables.open_file(  retrievalhdf  ) as f:
        
        tau = f.root.parameters.aerosol_tau.read()        
        tau_error = f.root.parameters.precision_bias_aerosol_tau.read()
        
        alh_base = f.root.parameters.fit_interval_base_altitude.read()
        alh_top  = f.root.parameters.fit_interval_top_altitude.read()
        alh_error =  f.root.parameters.precision_bias_cloud_base_pressure.read()
        
        iterations_data = f.root.state_vector_during_iterations.state_vector.read()

                
        specifications_simdata = f.root.specifications_sim
        for specs in specifications_simdata:
            if specs.name == 'solar_zenith_angle':
                sza = specs.read()
            if specs.name == 'viewing_nadir_angle':
                vna = specs.read()
            if specs.name == 'solar_azimuth_angle':
                saa = specs.read()
            if specs.name == 'viewing_azimuth_angle':
                vaa = specs.read()

        inputspecs = [sza,vna,saa,vaa]
                
        wavelengths = f.root.sun_normalized_radiance.wavelengths.read()
                          
        number_of_iterations = float(f.root._f_getattr('number_of_iterations'))
        solnconv = f.root._f_getattr('solution_has_converged')
        stoppedBoundary = f.root._f_getattr('stopped_at_boundary')
        time_per_iteration = float(f.root.attachments.Config_in._f_getattr('time'))/float(number_of_iterations)
        total_time = float(f.root.attachments.Config_in._f_getattr('time'))
        deg_of_freedom = float(f.root._f_getattr('DFS'))
        cost_function = float(f.root._f_getattr('chi2_state_vector'))
        
        lat = f.root.external_data.latitude.read()
        lon = f.root.external_data.longitude.read()
        clat = [float(i) for i in f.root.external_data.cornerlatitude.read()[0].split('_')]
        clon = [float(i) for i in f.root.external_data.cornerlongitude.read()[0].split('_')]

        scl = f.root.external_data.AtrackNumber.read()
        pxl = f.root.external_data.XtrackNumber.read()
        

  
        RESULTS = {}
    
        RESULTS['latitude'] = lat[0]
        RESULTS['longitude'] = lon[0]
        RESULTS['cornerLat'] = clat
        RESULTS['cornerLon'] = clon
        RESULTS['scanline'] = scl
        RESULTS['pixel'] = pxl
        
        RESULTS['wavelengths'] = wavelengths
        RESULTS['DFS'] = deg_of_freedom
        RESULTS['costFunction'] = cost_function
        
        RESULTS['tau'] = tau 
        RESULTS['alh_base'] = alh_base
        RESULTS['alh_top'] = alh_top
        RESULTS['tau_error'] = tau_error
        RESULTS['alh_error'] = alh_error

        RESULTS['iterData'] = iterations_data
        RESULTS['inputspecs'] = inputspecs
        RESULTS['numberOfIterations'] =number_of_iterations
        
        if solnconv == 'true':
            RESULTS['solnConvergence'] = 1
        else:
            RESULTS['solnConvergence'] = 0
        if stoppedBoundary == 'true':
            RESULTS['stoppedboundary'] = 1
        else:
            RESULTS['stoppedboundary'] = 0    
            
        RESULTS['timePerIteration'] = time_per_iteration
        RESULTS['totalTIME'] = total_time
    
    return RESULTS

    
from matplotlib.patches import Polygon


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
#----------------------------------------------------------------------------------------#
# inputs

#file_name_myd021km = 'MYD021KM.A2013189.1350.006.2013190155358.hdf'
file_name_myd03 = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/MODIS021KM/MOD03.A2010220.0850.005.2010227212417.hdf'


hdffile = 'MOD021KM.A2010220.0850.005.2010228013232.hdf'
hdfDataSource = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/data/MODIS021KM'
hdffile = os.path.join(hdfDataSource, hdffile)

file_myd021km = SD(hdffile, SDC.READ)

#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_1.0_1.5'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/nofilter_bb'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/midP_5'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/midP_10'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_midP10'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_midP10_divPointsAdjusted_1pthreshold'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/pmid_50'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/pmid_10'
outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/o2bo2a_o2bcontinuum'

days = ['2010_08_08']
boundingbox = [60.0, 29.0, 52.5, 45.0]

resDict = {}
#for i in os.listdir(outputfilepath):
#    if i.endswith('__ALH.h5'):
#        filename = os.path.join(outputfilepath,i)
#        day = '_'.join(i.split('_')[0:3])
#        
#        if day in days:
#            if not day in resDict:
#                
#                resDict[day] = {'Converged':[],'NonConverged':[]}
#                
#            else:
#    
#                R = read(filename)
#                if R['solnConvergence'] == 1:
#                    
#                    data_converged = [  R['latitude'],R['longitude'],
#                              R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
#                              R['alh_error'][1], R['cornerLat'] , R['cornerLon'] ]
#    
#    
#                    resDict[day]['Converged'].append(data_converged)
#        
#                if R['solnConvergence'] == 0:
#                    
#                    data_nonconverged = [  R['latitude'],R['longitude'],
#                              R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
#                              R['alh_error'][1], R['cornerLat'] , R['cornerLon'] ]
#                
#                    resDict[day]['NonConverged'].append(data_nonconverged)


for i in os.listdir(outputfilepath):
#    if i.endswith('__ALH.h5'):
    if i.endswith('__ALH.h5'):

        filename = os.path.join(outputfilepath,i)
        day = '_'.join(i.split('_')[0:3])
        
        if day in days:
            if not day in resDict:
                
                resDict[day] = {'Converged':[],'NonConverged':[]}
                R = read(filename)
#                R = read_topPres(filename)
                
                lat = R['latitude']
                lon = R['longitude']
    
                if R['solnConvergence'] == 1:


                    if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                            and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):


                    
                        data_converged = [  R['latitude'],R['longitude'],
                                  R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
                                  R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
    #                              R['inputspecs'],
                                  R['costFunction'], R['cornerLat'] , R['cornerLon'] ]

#                    data_converged = [  R['latitude'],R['longitude'],
#                              R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
##                              R['inputspecs'],
#                              R['costFunction'], R['cornerLat'] , R['cornerLon'] ]
    
    
                        resDict[day]['Converged'].append(data_converged)
        
                if R['solnConvergence'] == 0:
                    

                    if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                            and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):

                        data_nonconverged = [  R['latitude'],R['longitude'],
                              R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#                              R['inputspecs'],
                              R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]

#                    data_nonconverged = [  R['latitude'],R['longitude'],
#                              R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
##                              R['inputspecs'],
#                              R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]
                
                        resDict[day]['NonConverged'].append(data_nonconverged)
                        
            else:
    
                R = read(filename)
#                R = read_topPres(filename)
                
                lat = R['latitude']
                lon = R['longitude']
    
                if R['solnConvergence'] == 1:


                    if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                            and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):


                    
                        data_converged = [  R['latitude'],R['longitude'],
                                  R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
                                  R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
    #                              R['inputspecs'],
                                  R['costFunction'], R['cornerLat'] , R['cornerLon'] ]

#                    data_converged = [  R['latitude'],R['longitude'],
#                              R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
##                              R['inputspecs'],
#                              R['costFunction'], R['cornerLat'] , R['cornerLon'] ]
    
    
                        resDict[day]['Converged'].append(data_converged)
        
                if R['solnConvergence'] == 0:
                    

                    if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                            and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):

                        data_nonconverged = [  R['latitude'],R['longitude'],
                              R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#                              R['inputspecs'],
                              R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]

#                    data_nonconverged = [  R['latitude'],R['longitude'],
#                              R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
##                              R['inputspecs'],
#                              R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]
                
                        resDict[day]['NonConverged'].append(data_nonconverged)

od = collections.OrderedDict(sorted(resDict.items()))

dConv = np.array([])
dnConv = np.array([])


#----------------------------------------------------------------------------------------#

selected_sds = file_myd021km.select('EV_250_Aggr1km_RefSB')

selected_sds_attributes = selected_sds.attributes()

for key, value in selected_sds_attributes.iteritems():
	#print key, value
	if key == 'reflectance_scales':
		reflectance_scales_250_Aggr1km_RefSB = np.asarray(value)		
	if key == 'reflectance_offsets':
		reflectance_offsets_250_Aggr1km_RefSB = np.asarray(value)	

sds_data_250_Aggr1km_RefSB = selected_sds.get()

data_shape = sds_data_250_Aggr1km_RefSB.shape

along_track = data_shape[1]
cross_trak = data_shape[2]

#----------------------------------------------------------------------------------------#

selected_sds = file_myd021km.select('EV_500_Aggr1km_RefSB')

selected_sds_attributes = selected_sds.attributes()

for key, value in selected_sds_attributes.iteritems():
	if key == 'reflectance_scales':
		reflectance_scales_500_Aggr1km_RefSB = np.asarray(value)	
	if key == 'reflectance_offsets':
		reflectance_offsets_500_Aggr1km_RefSB = np.asarray(value)	

sds_data_500_Aggr1km_RefSB = selected_sds.get()

#----------------------------------------------------------------------------------------#

file_myd03 = SD(file_name_myd03, SDC.READ)
#
selected_sds = file_myd03.select('Latitude')
myd03_lat = selected_sds.get()
#
selected_sds = file_myd03.select('Longitude')
myd03_long = selected_sds.get()

#----------------------------------------------------------------------------------------#

data_shape = sds_data_250_Aggr1km_RefSB.shape

along_track = data_shape[1]
cross_trak = data_shape[2]

z = np.zeros((along_track, cross_trak,3))

for i in np.arange(along_track):
    for j in np.arange(cross_trak): 
        z[i,j,0] = ( sds_data_250_Aggr1km_RefSB[0,i,j] - \
        reflectance_offsets_250_Aggr1km_RefSB[0] ) * \
        reflectance_scales_250_Aggr1km_RefSB[0] 

for i in np.arange(along_track):
    for j in np.arange(cross_trak): 
        z[i,j,1] = ( sds_data_500_Aggr1km_RefSB[1,i,j] - \
        reflectance_offsets_500_Aggr1km_RefSB[1] ) * \
        reflectance_scales_500_Aggr1km_RefSB[1]  

for i in np.arange(along_track):
    for j in np.arange(cross_trak): 
        z[i,j,2] = ( sds_data_500_Aggr1km_RefSB[0,i,j] - \
        reflectance_offsets_500_Aggr1km_RefSB[0] ) * \
        reflectance_scales_500_Aggr1km_RefSB[0] 

z[ z > 1 ] = 1.0
z[ z < 0 ] = 0.0

#----------------------------------------------------------------------------------------#
# Rough estimation of latitude and longitude at granule center (long_0, lat_0)

lat_min = myd03_lat[0,0]
lat_max = myd03_lat[along_track-1,cross_trak-1]

lat_0 = lat_min + (lat_max - lat_min) / 2.

long_min = min(myd03_long[0,0],myd03_long[along_track-1,cross_trak-1])
long_max = max(myd03_long[0,0],myd03_long[along_track-1,cross_trak-1])

lon_0 = long_min + (long_max - long_min) / 2.

#----------------------------------------------------------------------------------------#
# Orthographic Map Projection

fig = plt.figure()

ax = fig.add_subplot(111)

ax.patch.set_facecolor((0,0,0))

m1 = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
#m1 = Basemap(width=2400000,height=1500000,resolution='l',projection='laea',lat_0=56.,lon_0=45.)

xpt0, ypt0 = m1(lon_0,lat_0) 

xpt1, ypt1 = m1(myd03_long[0,0],myd03_lat[0,0]) 
xpt2, ypt2 = m1(myd03_long[0,cross_trak-1],myd03_lat[0,cross_trak-1]) 
xpt3, ypt3 = m1(myd03_long[along_track-1,cross_trak-1], \
                myd03_lat[along_track-1,cross_trak-1])
xpt4, ypt4 = m1(myd03_long[along_track-1,0],myd03_lat[along_track-1,0])

llx = min(xpt1,xpt2,xpt3,xpt4) - xpt0  # lower left
lly = min(ypt1,ypt2,ypt3,ypt4) - ypt0

urx = max(xpt1,xpt2,xpt3,xpt4) - xpt0  # upper right
ury = max(ypt1,ypt2,ypt3,ypt4) - ypt0


#m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',\
#    llcrnrx=llx,llcrnry=lly,urcrnrx=urx,urcrnry=ury)

m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',\
    llcrnrx=llx,llcrnry=lly,urcrnrx=urx,urcrnry=ury)

#m = Basemap(width=2400000,height=1500000,resolution='l',projection='laea',lat_0=56.,lon_0=45.)

#try:
#    
#    img = m.imshow(np.rot90(np.fliplr(rgb_projected)), origin='lower')
#    
#except Exception:
    
    
    
x_igrid, y_igrid = m(myd03_long,myd03_lat)

x_igrid = x_igrid - xpt0
y_igrid = y_igrid - ypt0

z_igrid_01 = np.zeros((along_track, cross_trak))
z_igrid_02 = np.zeros((along_track, cross_trak))
z_igrid_03 = np.zeros((along_track, cross_trak))

for i in np.arange(2030):
    for j in np.arange(1354): 
        z_igrid_01[i,j] = z[i,j,0]
        z_igrid_02[i,j] = z[i,j,1]
        z_igrid_03[i,j] = z[i,j,2]

x1_igrid = x_igrid.ravel()
y1_igrid = y_igrid.ravel()
z_igrid_01 = z_igrid_01.ravel()
z_igrid_02 = z_igrid_02.ravel()
z_igrid_03 = z_igrid_03.ravel()

xy1_igrid = np.vstack((x1_igrid, y1_igrid)).T
xi, yi = np.mgrid[llx:urx:1000j, lly:ury:1000j]

z_01 = griddata(xy1_igrid, z_igrid_01, (xi, yi), method='cubic')
z_02 = griddata(xy1_igrid, z_igrid_02, (xi, yi), method='cubic')
z_03 = griddata(xy1_igrid, z_igrid_03, (xi, yi), method='cubic')

rgb_projected = np.zeros((1000, 1000,3))
for i in np.arange(1000):
    for j in np.arange(1000): 
        rgb_projected[i,j,0] = z_01[i,j]
        rgb_projected[i,j,1] = z_02[i,j]
        rgb_projected[i,j,2] = z_03[i,j]

#rgb_projected[ z > 1 ] = 1.0
#rgb_projected[ z < 0 ] = 0.0
whereAreNaNs = np.isnan(rgb_projected);
rgb_projected[whereAreNaNs] = 0.75;

img = m.imshow(np.rot90(np.fliplr(rgb_projected)), origin='lower')

m.drawcoastlines()

m.drawparallels(np.arange(-90.,120.,5.), color='k', labels=[True,False,False,False])
m.drawmeridians(np.arange(0.,420.,5.), color='k', labels=[False,False,False,True])

ax.set_xlabel("", fontsize=10)
ax.set_ylabel("", fontsize=10)       


    
#boundingbox = [60.0, 29.0, 52.5, 45.0]
lats = [52.5, 52.5, 60., 60. ]
lons = [45.,29.,29.,45.]

x, y = m( lons, lats )
xy = zip(x,y)
poly = Polygon( xy,facecolor='none',edgecolor='red')
ax.add_patch(poly)

#lats = [53.79, 59.71, 59.71, 53.79 ]
#lons = [34.29,30.96,30.96,34.29]
#
#x, y = m( lons, lats )
#xy = zip(x,y)
#poly = Polygon( xy,facecolor='none',edgecolor='blue',linewidth=6.0)
#ax.add_patch(poly)



if resDict[od.keys()[0]]['Converged']:
    
    convergedSoln = resDict[od.keys()[0]]['Converged']    
    
    dConv = np.vstack([j[0:9] for j in convergedSoln])
    cLatdata_c = np.vstack([j[-2] for j in convergedSoln])
    cLondata_c = np.vstack([j[-1] for j in convergedSoln])


dataC = dConv[:,[0,1,3]]

conv = []
for xi,i in enumerate(dataC):

    lats = [ cLatdata_c[xi][0],cLatdata_c[xi][1],cLatdata_c[xi][2],cLatdata_c[xi][3] ]
    lons = [ cLondata_c[xi][0],cLondata_c[xi][1],cLondata_c[xi][2],cLondata_c[xi][3] ]
    
    mlats,mlons = m(lons,lats)
        
    conv.append([mlats,mlons,i[2]])


if resDict[od.keys()[0]]['NonConverged']:

    nconvergedSoln = resDict[od.keys()[0]]['NonConverged']    
        
#    dnConv = np.vstack([j[0:9] for j in nconvergedSoln])
    dnConv = np.vstack([j[0:9] for j in nconvergedSoln])
    cLatdata_nc = np.vstack([j[-2] for j in nconvergedSoln])
    cLondata_nc = np.vstack([j[-1] for j in nconvergedSoln]) 

datanC = dnConv[:,[0,1,3]]

#N = np.unique(np.append(dataC[:,2],datanC[:,2])).shape[0]
N = np.unique(np.append(datanC[:,2],[])).shape[0]
       
conv = np.array(conv)
convPoly = []
for i in conv:
    convPoly.append(zip(i[0],i[1]))
    


ConvColl = PolyCollection(convPoly,array=conv[:,2],edgecolor='none',cmap=mpl.cm.jet)
im = ax.add_collection(ConvColl)


nconv = []
for xi,i in enumerate(datanC):


    lats = [ cLatdata_nc[xi][0],cLatdata_nc[xi][1],cLatdata_nc[xi][2],cLatdata_nc[xi][3] ]
    
    lons = [ cLondata_nc[xi][0],cLondata_nc[xi][1],cLondata_nc[xi][2],cLondata_nc[xi][3] ]
    
    mlats,mlons = m(lons,lats)
        
    nconv.append([mlats,mlons,i[2]])

                
nconv = np.array(nconv)
nconvPoly = []
for i in nconv:
    nconvPoly.append(zip(i[0],i[1]))

nConvColl = PolyCollection(nconvPoly,edgecolor=(1,1,1,0.4),facecolor=(1,1,1,0),cmap=mpl.cm.jet)
#nConvColl = PolyCollection(nconvPoly,array=conv[:,2],edgecolor='none',cmap=discrete_cmap(N, 'jet'))
ax.add_collection(nConvColl)




#colorlimits = [0.,4.]
colorlimits = [0.,conv[:,2].mean()+2*conv[:,2].std()]
ConvColl.set_array(conv[:,2])
ConvColl.set_clim(colorlimits)
#nConvColl.set_array(nconv[:,2])
cbar = fig.colorbar(im)
#ax_alh.cax.toggle_label(True)
cbar.set_label(r'$h_{mid}$ [km]')
#cbar.set_label(r'$\tau$ [-]')

ax.set_title('08 August, 2010')
#fig.savefig('/usr/people/nanda/Dropbox/Work Files/Papers/Paper1/Figures/GOME_2_Russian_Fires_2010__PF_ALH.png',format='png',transparent=True)
