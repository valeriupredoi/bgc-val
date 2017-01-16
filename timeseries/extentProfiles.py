#!/usr/bin/ipython 

#
# Copyright 2017, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license. 

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk

"""
.. module:: extentProfiles
   :platform: Unix
   :synopsis: A tool for producing a depth profile/Transect of contours.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from shelve import open as shOpen
from netCDF4 import Dataset,num2date
import os
from glob import glob
import shutil
from matplotlib import pyplot, gridspec
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import Locator
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy import img_transform, feature as cfeature 

#Specific local code:
import UKESMpython as ukp
from bgcvaltools.pftnames import getLongName
import timeseriesTools as tst 
import timeseriesPlots as tsp 


try:	
	defcmap = pyplot.cm.jet
	defcmapstr = 'jet'
except:	
	defcmap = viridis
	defcmapstr = 'viridis'	




#####
# Plan:
#
#	1 determine transect
#	2 load data according to mask
#	3 plot data as a plot
#	4 add land mask.



class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    from: https://stackoverflow.com/questions/20470892/how-to-place-minor-ticks-on-symlog-scale
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in xrange(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))
                                  
                                  

def makeLonSafe(lon):
	while True:
		if 0.<lon<=360.:return lon
		if lon<=0.:lon+=360.
		if lon> 360.:lon-=360.	
		
def makeLonSafeArr(lon):
	if lon.ndim == 2:
	 for l,lon1 in enumerate(lon):
	  for ll,lon2 in enumerate(lon1):
	   lon[l,ll] = makeLonSafe(lon2)
	 return lon
	if lon.ndim == 1:
	 for l,lon1 in enumerate(lon):
	   lon[l] = makeLonSafe(lon1)
	 return lon
	 

def findClosest(arr, z,debug=False):
	d = 10000.
	best = -1
	for i,zz in enumerate(arr.squeeze()):
		#print i,z,zz,depth.shape
		d2 = abs(z-zz)
		if d2<d:
		   d=d2
		   best = i
	if debug: print 'Closest: in situ:', z,'index:', best, 'distance:',d,', closest model:',arr.shape, arr[best]
	return best	

def loadKeyFromFile(fn):
	return os.path.basename(fn).replace('u-ad371o_1y_','').replace('_ptrc_T.nc','')[:4]

def contourplot(
		zmin = 0.,
		zmax = 400.,
		oxcutoff = 80.,
		):
	contours = [zmin,oxcutoff,zmax]
	plotKey = 'Equator'
	files = glob('/data/euryale7/scratch/ledm/UKESM/MEDUSA/u-ad371/u-ad371o_1y_*_ptrc_T.nc')
	#####
	# Load plot details
	pd = {}
	pd['Data'] = {'label':'Data','c': ['k',],  'lw':[2.,],'ls':['-',]}		
	for i,fn in enumerate(files):
		lw =1
		key = loadKeyFromFile(fn)
		color = defcmap(float(i)/float(len(files)))
		label = loadKeyFromFile(fn)
		pd[key] = {'label':label,'c': [color,],  'lw':[lw,],'ls':['-',],}
		
	fig = pyplot.figure()
	fig.set_size_inches(14,8)
			
	ax = pyplot.subplot(111)			
			
	nc = Dataset(files[0],'r')
	lats 	= nc.variables['nav_lat'][:]
	lons 	= nc.variables['nav_lon'][:]	
	if plotKey=='Equator':	loc 	= findClosest(lats.mean(1), 0.,debug=True)	# equator	
	x 	= makeLonSafeArr(lons[loc,:])	
	z 	= nc.variables['deptht'][::-1]*-1.
    	newX, newZ = np.meshgrid(x,z)		

	nc.close()

	

	for fn in files:
		nc = Dataset(fn,'r')
		o2 	= nc.variables['OXY'][0,::-1,loc,:].squeeze()
		nc.close()
		
		key = loadKeyFromFile(fn)		
	    	o2 = np.ma.masked_where((newX>359.)+(newX<0.5)+o2.mask,o2)
		o2 = np.ma.masked_where(np.ma.array(o2).mask+(o2<1E-10) +(o2>1e10),o2).squeeze()					

		im = ax.contour(
			newX,newZ,o2,
			contours,
	 		colors=pd[key]['c'],
	 		linewidths=pd[key]['lw'],
	 		linestyles=pd[key]['ls'],
	 		)
	pyplot.pcolormesh(newX, newZ,np.ma.masked_where(o2.mask,o2),cmap='Blues_r',vmin=zmin,vmax=zmax)
	pyplot.colorbar(orientation='horizontal')	 		


	
	#####
	# Add legend:
	addLegend = True
	if addLegend:
		legendSize = len(pd.keys())+1
		ncols = int(legendSize/25)+1
		box = ax.get_position()
		ax.set_position([box.x0,
				  box.y0 ,
				  box.width*(1.-0.1*ncols), 
				  box.height ])
	
		for i in sorted(pd.keys()):
			pyplot.plot(
				[], [], 
				c    = pd[i]['c'][0], 
				lw   = pd[i]['lw'][0], 
				ls   = pd[i]['ls'][0], 
				label= pd[i]['label'],
				)
												
		legd = ax.legend(loc='center left', ncol=ncols,prop={'size':10},bbox_to_anchor=(1., 0.5))
		legd.draw_frame(False) 
		legd.get_frame().set_alpha(0.)	
			
	
	pyplot.axhline(y= -500.,c='k',ls='--')
	pyplot.axhline(y=-1000.,c='k',ls='--') 

	
	pyplot.xticks([0.,60.,120.,180.,240.,300.,360.])
    	ax.set_xlim([0.,360.])	
	ax.set_yscale('symlog',linthreshy=1000.)
	ax.set_ylim([z.min(),-1.])
	#yaxis = pyplot.gca().yaxis
	#yaxis.set_minor_locator(MinorSymLogLocator(0.01))	
	pyplot.yticks([-10.,-100.,-500.,-1000.,-2000.,-5000.],['10','100','500','1000','2000','5000'])		
	pyplot.xlabel('Longitude, degress E', )#ha='center', va='center')
	pyplot.ylabel('Depth, m', ha='center', va='center', rotation='vertical')	
	pyplot.show()
	
	
contourplot()




