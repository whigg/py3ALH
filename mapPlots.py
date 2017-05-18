# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 10:01:29 2017

@author: nanda
"""


from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PolyCollection
import matplotlib as mpl

handy_python_scripts = '/nobackup/users/nanda/Code/pythonscriptlibrary'

import numpy as np
import tables, sys
sys.path.insert(0,handy_python_scripts)
import matplotlib.pyplot as plt


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

def makeplot_scatter(data,unit,plot_title,colorlimits,splot):

    """ data = [lat lon data]"""
    
    fontP = FontProperties()
    fontP.set_size('small')

    lon_r = [-25,50]
    lat_r = [25,65]
    plt.subplot(splot)
    
#    m = Basemap(lat_0=0.,lon_0=0.,llcrnrlon=round(lon_r[0]),urcrnrlon=round(lon_r[1]),llcrnrlat=round(lat_r[0]),urcrnrlat=round(lat_r[1]),projection='merc',resolution='l')#,area_thresh=10000,rsphere=6371200.)
    m = Basemap(width=5000000,height=2000000,resolution='l',projection='stere',lat_ts=20,lat_0=55.,lon_0=40.)

    m.drawcoastlines(linewidth=1)
    x,y = m(data[:,1],data[:,0])

    im = m.scatter(x,y,c=data[:,2],s=20,lw=0,marker='o')
    m.bluemarble()
#    m.drawcountries(linewidth=1)
    m.drawparallels(np.arange(-90.,91.,5), labels=[False, True, True, False],fontsize=8)
    m.drawmeridians(np.arange(-180.,180.,5), labels=[True, False, False, True],fontsize=8)
    plt.title(plot_title,fontproperties=fontP)

    ax = plt.gca()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.25)
    
    plt.clim(colorlimits)

    plt.colorbar(im, orientation='horizontal',cax=cax, label=unit)

def makeplot_polygon(fig,data,cLatdata,cLondata,colorlimits,fill):
    
    ax = fig.add_subplot(111)
    m = Basemap(width=4200000,height=2100000,resolution='l',projection='stere',lat_ts=20,lat_0=55.,lon_0=40.)
    m.drawcoastlines()
    m.drawmapboundary()
    m.bluemarble()
    m.drawcountries(linewidth=1)
    m.drawparallels(np.arange(-90.,91.,5), labels=[False, True, True, False],fontsize=8)
    m.drawmeridians(np.arange(-180.,180.,5), labels=[True, False, False, True],fontsize=8)
    
    conv = []
    for xi,i in enumerate(data):
    
    
        lats = [ cLatdata[xi][0],cLatdata[xi][1],cLatdata[xi][2],cLatdata[xi][3] ]
        
        lons = [ cLondata[xi][0],cLondata[xi][1],cLondata[xi][2],cLondata[xi][3] ]
        
        mlats,mlons = m(lons,lats)
            
        conv.append([mlats,mlons,i[2]])
    
                    
    conv = np.array(conv)
    convPoly = []
    for i in conv:
        convPoly.append(zip(i[0],i[1]))
    

    if fill:
        ConvColl = PolyCollection(convPoly,array=conv[:,2],edgecolor='black',cmap=mpl.cm.jet)
        ax.add_collection(ConvColl)
    
        ConvColl.set_array(conv[:,2])
        ConvColl.set_clim(colorlimits)
        
    else:

        ConvColl = PolyCollection(convPoly,edgecolor=(1,1,1,1),facecolor=(1,1,1,0),cmap=mpl.cm.jet)
        ax.add_collection(ConvColl)
#        ConvColl.set_array(conv[:,2])
        
        

    
#    cb = plt.colorbar(ConvColl, shrink=0.5,orientation="horizontal")
    
    return ax,ConvColl

def makesubplot_polygon(fig,splot,data,cLatdata,cLondata,colorlimits,fill):
    
    ax = fig.add_subplot(splot)
    m = Basemap(width=4200000,height=2100000,resolution='l',projection='stere',lat_ts=20,lat_0=55.,lon_0=40.)
    m.drawcoastlines()
    m.drawmapboundary()
    m.bluemarble()
    m.drawcountries(linewidth=1)
    m.drawparallels(np.arange(-90.,91.,5), labels=[False, True, True, False],fontsize=8)
    m.drawmeridians(np.arange(-180.,180.,5), labels=[True, False, False, True],fontsize=8)
    
    conv = []
    for xi,i in enumerate(data):
    
    
        lats = [ cLatdata[xi][0],cLatdata[xi][1],cLatdata[xi][2],cLatdata[xi][3] ]
        
        lons = [ cLondata[xi][0],cLondata[xi][1],cLondata[xi][2],cLondata[xi][3] ]
        
        mlats,mlons = m(lons,lats)
            
        conv.append([mlats,mlons,i[2]])
    
                    
    conv = np.array(conv)
    convPoly = []
    for i in conv:
        convPoly.append(zip(i[0],i[1]))
    

    if fill:
        ConvColl = PolyCollection(convPoly,array=conv[:,2],edgecolor='black',cmap=mpl.cm.jet)
        ax.add_collection(ConvColl)
    
        ConvColl.set_array(conv[:,2])
        ConvColl.set_clim(colorlimits)
        
    else:

        ConvColl = PolyCollection(convPoly,edgecolor=(1,1,1,1),facecolor=(1,1,1,0),cmap=mpl.cm.jet)
        ax.add_collection(ConvColl)
        
    return ax,ConvColl

def makesubplot_imageGrid(ax,data,cLatdata,cLondata,colorlimits,fill):
    
    m = Basemap(width=2400000,height=1500000,resolution='l',projection='laea',lat_0=56.,lon_0=40.)
    m.drawcoastlines()
    m.drawmapboundary()
    m.bluemarble()
    m.drawcountries(linewidth=1)
    m.drawparallels(np.arange(-90.,91.,5.), labels=[False, True, True, False],fontsize=8)
    m.drawmeridians(np.arange(-180.,180.,10.), labels=[True, False, False, True],fontsize=8)
    
    conv = []
    for xi,i in enumerate(data):
    
    
        lats = [ cLatdata[xi][0],cLatdata[xi][1],cLatdata[xi][2],cLatdata[xi][3] ]
        
        lons = [ cLondata[xi][0],cLondata[xi][1],cLondata[xi][2],cLondata[xi][3] ]
        
        mlats,mlons = m(lons,lats)
            
        conv.append([mlats,mlons,i[2]])
    
                    
    conv = np.array(conv)
    convPoly = []
    for i in conv:
        convPoly.append(list(zip(i[0],i[1])))
    

    if fill:
        ConvColl = PolyCollection(convPoly,array=conv[:,2],edgecolor='none',cmap=mpl.cm.jet)
        im = ax.add_collection(ConvColl)
    
        ConvColl.set_array(conv[:,2])
        if colorlimits:
            ConvColl.set_clim(colorlimits)
        
    else:
        
        ConvColl = PolyCollection(convPoly,edgecolor=(1,1,1,1),facecolor=(1,1,1,0),cmap=mpl.cm.jet)
        im = ax.add_collection(ConvColl)

    return im