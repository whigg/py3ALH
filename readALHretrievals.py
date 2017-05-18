# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:21:32 2017

@author: nanda
"""

import tables, os

import mapPlots
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import collections
from matplotlib.font_manager import FontProperties

def read(retrievalhdf):
    
    with tables.open_file(  retrievalhdf,'r'  ) as f:
        
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
        solnconv = f.root._f_getattr('solution_has_converged').decode()
        stoppedBoundary = f.root._f_getattr('stopped_at_boundary').decode()
        time_per_iteration = float(f.root.attachments.Config_in._f_getattr('time'))/float(number_of_iterations)
        total_time = float(f.root.attachments.Config_in._f_getattr('time'))
        deg_of_freedom = float(f.root._f_getattr('DFS'))
        cost_function = float(f.root._f_getattr('chi2_state_vector'))
        
        lat = f.root.external_data.latitude.read()
        lon = f.root.external_data.longitude.read()
        clat = [float(i) for i in f.root.external_data.cornerlatitude.read()[0].decode().split('_')]
        clon = [float(i) for i in f.root.external_data.cornerlongitude.read()[0].decode().split('_')]

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
        RESULTS['numberOfIterations'] =int(number_of_iterations)

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

def read_topPres(retrievalhdf):
    
    with tables.open_file(  retrievalhdf  ) as f:
        
        tau = f.root.parameters.aerosol_tau.read()        
        tau_error = f.root.parameters.precision_bias_aerosol_tau.read()
        
        alh_top  = f.root.parameters.fit_interval_top_altitude.read()
        alh_error =  f.root.parameters.precision_bias_cloud_top_pressure.read()
        
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

#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/nofilter_bb'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/filtered_bb'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_1.0_1.5'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/fixedB'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/midP_5'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia/top_bot'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_midP10'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_midP10_divPointsAdjusted'
outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/GOME2_Aug_Russia_PF/doublePF_midP10_divPointsAdjusted_1pthreshold'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/pmid_50'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/pmid_10'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/onlyO2B'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/onlyO2B_ang0'
#outputfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/output/gome2_aug_pf_o2b/o2bo2a_o2bcontinuum'


jobfilepath = '/net/pc150230/nobackup/users/nanda/Projects/GOME2_orbits_LER_BRDF/jobfiles'
jobDict = {}


for i in os.listdir(jobfilepath):
    day = '_'.join(i.split('_')[3:6])
    if not day in jobDict:
        jobDict[day] = []
    
    else:
        jobDict[day].append([int(j) for j in i.split('_')[1:3] ] ) 

#days = ['2010_08_07','2010_08_08']
days = ['2010_08_08']
boundingbox = [60.0, 29.0, 52.5, 45.0]

resDict = {}
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
#
#                        data_converged = [  R['latitude'],R['longitude'],
#                                  R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                                  R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#    #                              R['inputspecs'],
#                                  R['costFunction'], R['cornerLat'] , R['cornerLon'] ]
#    
    
                        resDict[day]['Converged'].append(data_converged)
        
                if R['solnConvergence'] == 0:
                    

                    if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                            and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):

                        data_nonconverged = [  R['latitude'],R['longitude'],
                              R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#                              R['inputspecs'],
                              R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]

#                        data_nonconverged = [  R['latitude'],R['longitude'],
#                                  R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                                  R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#    #                              R['inputspecs'],
#                                  R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]
                
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
#
#                        data_converged = [  R['latitude'],R['longitude'],
#                                  R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                                  R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#    #                              R['inputspecs'],
#                                  R['costFunction'], R['cornerLat'] , R['cornerLon'] ]
#    
    
                        resDict[day]['Converged'].append(data_converged)
        
                if R['solnConvergence'] == 0:
                    

                    if (len(boundingbox) > 1) and ( (np.max([boundingbox[0],boundingbox[2]]) - lat > 0.) and (np.min([boundingbox[0],boundingbox[2]]) - lat < 0. ) ) \
                            and ( ( np.max([boundingbox[1],boundingbox[3]]) - lon > 0. ) and ( np.min([boundingbox[1],boundingbox[3]]) - lon < 0. ) ):

                        data_nonconverged = [  R['latitude'],R['longitude'],
                              R['tau'][2],0.5*(R['alh_base'][2]+R['alh_top'][2]), R['tau_error'][1], 
                              R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#                              R['inputspecs'],
                              R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]

#                        data_nonconverged = [  R['latitude'],R['longitude'],
#                                  R['tau'][2],R['alh_top'][2], R['tau_error'][1], 
#                                  R['alh_error'][1],R['scanline'][0],R['pixel'][0], R['numberOfIterations'],
#    #                              R['inputspecs'],
#                                  R['costFunction'] , R['cornerLat'] , R['cornerLon'] ]
                
                        resDict[day]['NonConverged'].append(data_nonconverged)

od = collections.OrderedDict(sorted(resDict.items()))
odkeys = list(od.keys())

""" Plot ALH data """
fig_alh = plt.figure()
grid_alh = ImageGrid(fig_alh, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,1),
                 axes_pad=0.65,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="2%",
                 cbar_pad=0.55,)


                 
fontP = FontProperties()
fontP.set_size('small')
                 
for i,ax_alh in enumerate(grid_alh):

    dConv = np.array([])
    dnConv = np.array([])
    
    ax_alh.set_title('{0} August, 2010'.format(int(odkeys[i].split('_')[2])), fontproperties=fontP )
    
    if resDict[odkeys[i]]['Converged']:
        
        convergedSoln = resDict[odkeys[i]]['Converged']    
        
        dConv = np.vstack([j[0:10] for j in convergedSoln])
        cLatdata_c = np.vstack([j[-2] for j in convergedSoln])
        cLondata_c = np.vstack([j[-1] for j in convergedSoln])

    if resDict[odkeys[i]]['NonConverged']:

        nconvergedSoln = resDict[odkeys[i]]['NonConverged']    
            
        dnConv = np.vstack([j[0:10] for j in nconvergedSoln])
        cLatdata_nc = np.vstack([j[-2] for j in nconvergedSoln])
        cLondata_nc = np.vstack([j[-1] for j in nconvergedSoln])    

    plt.sca(ax_alh)

    if dConv.any():
        
        im_alh = mapPlots.makesubplot_imageGrid(ax_alh,dConv[:,[0,1,3]],cLatdata_c,cLondata_c,[],True)
                
    if dnConv.any():
        
        mapPlots.makesubplot_imageGrid(ax_alh,dnConv[:,[0,1,3]],cLatdata_nc,cLondata_nc,[],False) 
        
        
ax_alh.cax.colorbar(im_alh)
ax_alh.cax.toggle_label(True)
ax_alh.cax.colorbar(im_alh).set_label_text(r'$h_{mid}$ [km]')
#fig_alh.tight_layout()
#fig_alh.savefig('/usr/people/nanda/Dropbox/Work Files/Papers/Paper1/Figures/GOME_2_Russian_Fires_2010_ALH.png',format='png',transparent=True)
#
##
# histograms

hist_alh = plt.figure()
hist_Grid_alh = ImageGrid(hist_alh, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,1),
                 axes_pad=0.65,share_all=True,aspect=False
                    )

                 
fontP = FontProperties()
fontP.set_size('small')
                 
for i,ax_alh_hist in enumerate(hist_Grid_alh):

    dConv = np.array([])
    
    ax_alh_hist.set_title('{0} August, 2010'.format(int(odkeys[i].split('_')[2])), fontproperties=fontP )
    
    if resDict[odkeys[i]]['Converged']:
        
        convergedSoln = resDict[odkeys[i]]['Converged']    
        
        dConv = np.vstack([j[0:4] for j in convergedSoln])

#    plt.sca(ax_alh_hist)

    if dConv.any():
        
        ax_alh_hist.hist(dConv[:,3],20,facecolor='gray', align='mid')
    ax_alh_hist.set_xlabel(r'$h_{mid}$ $[km]$')        
        
#hist_alh.savefig('/usr/people/nanda/Dropbox/Work Files/Papers/Paper1/Figures/GOME_2_Russian_Fires_20100808_ALH_hist_pmid10__PF.png',format='png',transparent=True)


""" colocation histogram """

conv_collocation = []
for i in od[odkeys[0]]['Converged']:
    
    if (i[0] > 53.79 and i[0]<59.71) and (i[1]>30.96 and i[1]<34.29):
        conv_collocation.append([i[0],i[1],i[2],i[3]])

conv_collocation = np.vstack(conv_collocation)

nconv_collocation = []
for i in od[odkeys[0]]['NonConverged']:
    
    if (i[0] > 53.79 and i[0]<59.71) and (i[1]>30.96 and i[1]<34.29):
        nconv_collocation.append([i[0],i[1],i[2],i[3]])

if nconv_collocation:
    nconv_collocation = np.vstack(nconv_collocation)


plt.figure()
plt.hist(conv_collocation[:,3])



""" Plot AOT data """

fig_aot = plt.figure()
grid_aot = ImageGrid(fig_aot, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,1),
                 axes_pad=0.65,
                 share_all=True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="2%",
                 cbar_pad=0.55,)


for i,ax_aot in enumerate(grid_aot):

    dConv = np.array([])
    dnConv = np.array([])
    
    ax_aot.set_title('{0} August, 2010'.format(int(odkeys[i].split('_')[2])), fontproperties=fontP )
    plt.sca(ax_aot)

    if resDict[odkeys[i]]['Converged']:
        
        convergedSoln = resDict[odkeys[i]]['Converged']    
        
        dConv = np.vstack([j[0:4] for j in convergedSoln])
        cLatdata_c = np.vstack([j[-2] for j in convergedSoln])
        cLondata_c = np.vstack([j[-1] for j in convergedSoln])

    if resDict[odkeys[i]]['NonConverged']:

        nconvergedSoln = resDict[odkeys[i]]['NonConverged']    
            
        dnConv = np.vstack([j[0:4] for j in nconvergedSoln])
        cLatdata_nc = np.vstack([j[-2] for j in nconvergedSoln])
        cLondata_nc = np.vstack([j[-1] for j in nconvergedSoln])   

    if dConv.any():
        
        im_aot = mapPlots.makesubplot_imageGrid(ax_aot,dConv[:,[0,1,2]],cLatdata_c,cLondata_c,[],True)
                
    if dnConv.any():
        
        mapPlots.makesubplot_imageGrid(ax_aot,dnConv[:,[0,1,2]],cLatdata_nc,cLondata_nc,[],False) 

ax_aot.cax.colorbar(im_aot)
ax_aot.cax.toggle_label(True)
ax_aot.cax.colorbar(im_aot).set_label_text(r'$\tau$ [-]')
#fig_aot.tight_layout()
##fig_aot.savefig('/usr/people/nanda/Dropbox/Work Files/Papers/Paper1/Figures/GOME_2_Russian_Fires_2010_AOT.png',format='png',transparent=True)
#
# histograms

hist_aot = plt.figure()
hist_Grid_aot = ImageGrid(hist_aot, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,1),
                 axes_pad=0.65,share_all=True,aspect=False
                    )

                 
fontP = FontProperties()
fontP.set_size('small')
                 
for i,ax_aot_hist in enumerate(hist_Grid_aot):

    dConv = np.array([])
    dnConv = np.array([])
    
    ax_aot_hist.set_title('{0} August, 2010'.format(int(odkeys[i].split('_')[2])), fontproperties=fontP )
    
    if resDict[odkeys[i]]['Converged']:
        
        convergedSoln = resDict[odkeys[i]]['Converged']    
        
        dConv = np.vstack([j[0:4] for j in convergedSoln])

#    plt.sca(ax_alh_hist)

    if dConv.any():
        
        ax_aot_hist.hist(dConv[:,2],40,facecolor='gray', align='mid')
    ax_aot_hist.set_xlabel(r'$\tau$ $[-]$')

#hist_aot.savefig('/usr/people/nanda/Dropbox/Work Files/Papers/Paper1/Figures/GOME_2_Russian_Fires_2010_AOT_hist_nofilter_pmid10__PF.png',format='png',transparent=True)
