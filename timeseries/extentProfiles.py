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
import paths	

try:	
	defcmap = pyplot.cm.jet
	defcmapstr = 'jet'
except:	
	defcmap = viridis
	defcmapstr = 'viridis'	


zonalCuts 	= ['Equator','10 N', '10 S',]
MeridionalCuts 	= ['Atlantic','Atlantic28W', 'Pacific135W']




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


def LoadZonalTransect(data, transectLat,lats,lons,depths):
	if dlats.ndim ==1:
		loc 	= findClosest(lats, transectLat)
		lons 	= makeLonSafeArr(lons)
		depths 	=depths[::-1]*-1.	    			

	    	newX, newZ = np.meshgrid(lons,depths)
	    	if data.ndim ==4:	dat = data[:,::-1,dloc,:].squeeze()
	    	return dat, newX,newY
	
	if dlats.ndim ==2:
		dloc 	= findClosest(dlats.mean(1), plotKeyDict[plotKey],debug=True)
		dx 	= makeLonSafeArr(dlons[dloc,:])				


def contourplot(
		jobID,
		name,
		modelfiles,
		datafile,
		modeldetails,
		datadetails,
		modelcoords,
		datacoords,		
		zmin = 0.,
		zmax = 400.,
		oxcutoff = 80.,
		plotKey = 'Equator',
		title	= '',
		filename = '',
		cbarlabel = '',
		):
		
	contours = [zmin,oxcutoff,zmax]
	plotKeyDict = {'Equator':0.,'10 N':10., '10 S':-10.,'Atlantic28W':-28., 'Pacific135W':-135.}
	zonalCuts 	= ['Equator','10 N', '10 S',]
	MeridionalCuts 	= ['Atlantic','Atlantic28W', 'Pacific135W']

	
	#####
	# Load plot details
	pd = {}
	pd['Data'] = {'label':'Data','c': ['k',],  'lw':[2.,],'ls':['-',]}		
	for i,fn in enumerate(modelfiles):
		lw =1
		key = loadKeyFromFile(fn)
		color = defcmap(float(i)/float(len(modelfiles)))
		label = loadKeyFromFile(fn)
		pd[key] = {'label':label,'c': [color,],  'lw':[lw,],'ls':['-',],}

	#####
	# Load data lats/lons, and data
	print "Loading:",datafile
	dnc = Dataset(datafile,'r')
	dlats 	= dnc.variables[datacoords['lat']][:]
	dlons 	= dnc.variables[datacoords['lon']][:]	
	do2 	= ukp.extractData(dnc, datadetails)
	
	if plotKey in zonalCuts:
		
		
		if dlats.ndim ==2:
			dloc 	= findClosest(dlats.mean(1), plotKeyDict[plotKey],debug=True)
			dx 	= makeLonSafeArr(dlons[dloc,:])				
		if dlats.ndim ==1:
			dloc 	= findClosest(dlats,         plotKeyDict[plotKey],debug=True)
			dx 	= makeLonSafeArr(dlons)	
	if plotKey in MeridionalCuts:
		if dlons.ndim ==2:
			dloc 	= findClosest(dlons.mean(0), plotKeyDict[plotKey],debug=True)
			dx 	= dlat[:,dloc]
		if dlons.ndim ==1:
			dloc 	= findClosest(dlons,         plotKeyDict[plotKey],debug=True)
			dx 	= dlats
	dz 	= dnc.variables[datacoords['z']][::-1]*-1.
    	dnewX, dnewZ = np.meshgrid(dx,dz)
    	

	if plotKey in zonalCuts:	
		if do2.ndim==4:	do2 = do2[0,::-1,dloc,:].squeeze()
	    	do2 	= np.ma.masked_where((dnewX>=359.2)+(dnewX<0.2)+do2.mask,do2)		
	if plotKey in MeridionalCuts	:
		if do2.ndim==4:	do2 = do2[0,::-1,:,dloc].squeeze()
				
	do2 	= np.ma.masked_where(np.ma.array(do2).mask+(do2<1E-10) +(do2>1e10),do2).squeeze()
	dnc.close()
	

	#####
	# Load model lats/lons
	print "Loading:",modelfiles[0]	
	nc = Dataset(modelfiles[0],'r')
	lats 	= nc.variables[modelcoords['lat']][:]
	lons 	= nc.variables[modelcoords['lon']][:]
	if plotKey in zonalCuts:	
		loc 	= findClosest(lats.mean(1), plotKeyDict[plotKey],debug=True)
		x 	= makeLonSafeArr(lons[loc,:])	
	if plotKey in MeridionalCuts:	
		loc 	= findClosest(lons.mean(0), plotKeyDict[plotKey],debug=True)
		x 	= lats[:,loc]
		
	z 	= nc.variables[modelcoords['z']][::-1]*-1.
    	newX, newZ = np.meshgrid(x,z)	
	nc.close()
	
			
	#####
	# Start making the figure	
	fig = pyplot.figure()
	fig.set_size_inches(10,6)
	ax = pyplot.subplot(111)			

	#####
	# Add data as a colormesh				
	cmesh = pyplot.pcolormesh(dnewX, dnewZ,np.ma.masked_where(do2.mask,do2),cmap='Blues_r',vmin=zmin,vmax=zmax)
	im = ax.contour(
		dnewX,dnewZ,do2,
		contours,
 		colors=pd['Data']['c'],
 		linewidths=pd['Data']['lw'],
 		linestyles=pd['Data']['ls'],
 		)	

	#####
	# Add model data as a colormesh
	for fn in modelfiles:
		nc = Dataset(fn,'r')
		
		o2 	= ukp.extractData(nc, modeldetails)		
		if plotKey in zonalCuts:	
			if o2.ndim==4:	o2 = o2[0,::-1,dloc,:].squeeze()
		    	o2 = np.ma.masked_where((newX>359.2)+(newX<0.2)+o2.mask,o2)			
		if plotKey in MeridionalCuts:
			if o2.ndim==4:	o2 = o2[0,::-1,:,dloc].squeeze()
		nc.close()
		
		key = loadKeyFromFile(fn)		

		o2 = np.ma.masked_where(np.ma.array(o2).mask+(o2<1E-10) +(o2>1e10),o2).squeeze()					

		im = ax.contour(
			newX,newZ,o2,
			contours,
	 		colors=pd[key]['c'],
	 		linewidths=pd[key]['lw'],
	 		linestyles=pd[key]['ls'],
	 		)

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
			
	c1 = pyplot.colorbar(cmesh, orientation='horizontal',)#pad=0.05,shrink=0.9 )
	if len(cbarlabel)>0: c1.set_label(cbarlabel)
		

	####
	# Sort out X axis:
	if plotKey in zonalCuts:	
		pyplot.xticks([0.,60.,120.,180.,240.,300.,360.])
	    	ax.set_xlim([0.,360.])	
		pyplot.xlabel('Longitude, degress E', )	    	
	if plotKey in MeridionalCuts:	
		pyplot.xlabel('Latitude, degress N', )		
	    	ax.set_xlim([-90.,90.])	
		pyplot.xticks([-90.,-60.,-30.,0.,30.,60.,90.])	    	
	####
	# Sort out Y axis:	
	pyplot.axhline(y= -500.,c='k',ls='--')
	pyplot.axhline(y=-1000.,c='k',ls='--') 		
	ax.set_yscale('symlog',linthreshy=1000.)
	ax.set_ylim([z.min(),-1.])
	pyplot.yticks([-10.,-100.,-500.,-1000.,-2000.,-5000.],['10','100','500','1000','2000','5000'])	
	pyplot.ylabel('Depth, m', ha='center', va='center', rotation='vertical')

	pyplot.title(title)
	
	#####
	#Save figure
	if filename=='':	filename=ukp.folder(['images',jobID,'OMZ'])+'-'.join([name,plotKey,'contour',str(int(oxcutoff))])+'.png'
	print "saving",filename
	pyplot.savefig(filename )		
	pyplot.close()		



def run():
	jobID = 'u-ad371'
	modelfiles 	= glob('/data/euryale7/scratch/ledm/UKESM/MEDUSA/'+jobID+'/'+jobID+'o_1y_*_ptrc_T.nc')	
	datafile 	= paths.WOAFolder_annual+'woa13_all_o00_01.nc'
	name 		= 'OMZExtent'

	medusaCoords 	= {'t':'time_counter', 	'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	woaCoords 	= {'t':'index_t', 	'z':'depth',  'lat': 'lat', 	'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	modeldetails 	= {'name': name, 'vars':['OXY',],  'convert': ukp.NoChange,'units':'mmol O2/m^3'}	
	datadetails  	= {'name': name, 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol O2/m^3'}	

	cbarlabel	= 'WOA Oxygen concentration, mmol O2/m^3'
	
	for plotKey in ['Pacific135W','Atlantic28W','Equator', '10 N', '10 S',]:
		oxcutoffs = [20.,50.,80.,]
		for oxcutoff in oxcutoffs:
			filename=ukp.folder(['images',jobID,'OMZ'])+'-'.join([name,plotKey,'contour',str(int(oxcutoff))])+'.png'
			title	= ' '.join([jobID, getLongName(plotKey),getLongName(name)+',',str(oxcutoff),modeldetails['units'] ])
			contourplot(
				jobID,
				name,
				modelfiles,datafile,
				modeldetails,datadetails,
				medusaCoords,woaCoords,
				plotKey = plotKey,
				oxcutoff = oxcutoff,
				filename = filename,
				title 	= title,
				cbarlabel = cbarlabel)

if __name__=="__main__":
	run() 
	print 'The end.'





