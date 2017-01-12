#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
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
.. module:: extentMaps
   :platform: Unix
   :synopsis: A tool for producing a map of contours.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from shelve import open as shOpen
from netCDF4 import Dataset,num2date
import os
import shutil
from matplotlib import pyplot, gridspec
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
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
	
	

def regrid(data,lat,lon):
        #print 'old regid shape:',data.shape, lat.shape,lon.shape               
    	nX = np.arange(-179.5,180.5,2.5)
    	nY = np.arange( -89.5, 90.5,2.5)
    	newLon, newLat = np.meshgrid(nX,nY)
    	
    	crojp1 = ccrs.PlateCarree(central_longitude=180.0, )#central_latitude=300.0)
    	crojp2 = ccrs.PlateCarree(central_longitude=180.0, )#central_latitude=300.0)

    	a = img_transform.regrid(data,
    			     source_x_coords=lon,
                             source_y_coords=lat,
                             source_cs=crojp1,
                             target_proj=crojp2,
                             target_x_points=newLon,
                             target_y_points=newLat
                             )
                             
        # print 'newregid shape:',a.shape                     
	return crojp2, a, newLon,newLat
	


def remask(data,lat,lon,regiddata, regridLat,regridLon, res=1.):
	#####
	# This is neccesairy as the img_tranform above doesn't see gaps in the data and fills it in.
	# We could do with re
	if len(data) == len(lat) ==len(lon): pass
	else:	assert "The data,lat,lon aren't the same length!"
	
	print data.shape,lat.shape,lon.shape,regridLat.shape,regridLon.shape
	newmask = np.zeros_like(regiddata)
	for (i,j), la in np.ndenumerate(regridLat):
		latdiff = np.abs(lat- la)
		londiff = np.abs(lon - regridLon[i,j])
		mdistance = np.min(np.power(latdiff,2.)+np.power(londiff,2.))
		if mdistance > res: newmask[i,j]+=1.
			
	return np.ma.masked_where( newmask>=1., regiddata)


			
def makeExtentPlot(
	modeldata,
	modellat,
	modellon,
	realdata,
	reallat,
	reallon,
	contours,
	filename,
	title='',
	labels='',
	zrange = '',
	colourmesh = True,
	contour	= True,
	):
	
	print modeldata.shape, modellat.shape,modellon.shape
	
	
	zmin = min([modeldata.min(),realdata.min()])
	zmax = max([modeldata.max(),realdata.max()])	

	
	if zrange in ['', 'auto',]:
		zrange = [zmin,zmax  ]
	
	if len(contours)==1:
		print "makeExtentPlot:\t adding min and max to contours."
		if zmin < contours[0] < zmax:
			contours = [zmin, contours[0],zmax]

	fig = pyplot.figure()
	fig.set_size_inches(14,8)
	ax = pyplot.subplot(111,projection=ccrs.PlateCarree(central_longitude=0., ))

	crojp2, mdregid, newmLon, newmLat  = regrid(modeldata,modellat,modellon)
	crojp2, rdregid, newdLon, newdLat  = regrid(realdata,reallat,reallon)

	mdregid = remask(modeldata,modellat,modellon,mdregid, newmLat,newmLon)
	
#	mdregid = np.ma.masked_where((mdregid < modeldata.min()) + (mdregid > modeldata.max()), mdregid)
#	rdregid = np.ma.masked_where((rdregid < realdata.min())  + (rdregid > realdata.max()) , rdregid)	
	
#	print 'overall:',zmin,zmax, 'model', [modeldata.min(),modeldata.max()], 'data', [realdata.min(),realdata.max()]
#	print 'regridded:', 'model', [mdregid.min(),mdregid.max()], 'data', [rdregid.min(),rdregid.max()]	
	
 	if colourmesh:
 		pyplot.pcolormesh(newmLon, newmLat,mdregid,transform=ccrs.PlateCarree(),vmin=zrange[0],vmax=zrange[1])
		pyplot.colorbar()
		
	if contour:
	 	ax.contour(newmLon,newmLat,mdregid,contours,colors=['darkblue',],linewidths=[1.5,],linestyles=['-',],transform=ccrs.PlateCarree(),zorder=1)
 		#ax.contour(newdLon,newdLat,rdregid,contours,colors=['black',],   linewidths=[1.5,],linestyles=['-',],transform=ccrs.PlateCarree(),zorder=1)
 		
	ax.add_feature(cfeature.LAND,  facecolor='white',zorder=2)
	ax.coastlines(lw=0.5,zorder=3)
	pyplot.title(title)

	print "saving",filename
	pyplot.savefig(filename )		
	pyplot.close()

	
def interannualExtendMap(
		modeldata,
		modellat,
		modellon,
		realdata,
		reallat,
		reallon,
		contours,
		filename,
		title='',
		labels='',
		zrange = '',
		#colourmesh = True,
		showdata	= True,
		addLegend	= True,
		):	
	keys = sorted(	modeldata.keys())
		
	zmin = [realdata.min(),]
	zmax = [realdata.max(),]
	for key in keys:
		zmin.append(np.ma.min(modeldata[key]))
		zmax.append(np.ma.max(modeldata[key]))
	zmin = np.ma.min(zmin)
	zmax = np.ma.max(zmax)	

	####
	# Add plot details
	pd = {}
	pd['Data'] = {'label':'Data','c': ['k',],  'lw':[2,],'ls':['-',]}		
	for i,key in enumerate(keys):
		lw =1
		color = defcmap(float(i)/float(len(keys)))
		label = key
		pd[key] = {'label':label,'c': [color,],  'lw':[lw,],'ls':['-',],}
		
		
	
	if zrange in ['', 'auto',]:
		zrange = [zmin,zmax  ]	
	
	if len(contours)==1:
		print "makeExtentPlot:\t adding min and max to contours."
		if zmin < contours[0] < zmax:
			contours = [zmin, contours[0],zmax]

	fig = pyplot.figure()
	fig.set_size_inches(14,8)
	ax = pyplot.subplot(111,projection=ccrs.PlateCarree(central_longitude=0., ))

	crojp2, rdregid, newdLon, newdLat  = regrid(realdata,reallat,reallon)

 	#####
 	# Add model contours	
	for key in keys:
		crojp2, mdregid, newmLon, newmLat  = regrid(modeldata[key],modellat[key],modellon[key])
		mdregid = remask(modeldata[key],modellat[key],modellon[key],mdregid, newmLat,newmLon)
	 	ax.contour(
	 		newmLon,newmLat,mdregid,
	 		contours,
	 		colors=pd[key]['c'],
	 		linewidths=pd[key]['lw'],
	 		linestyles=pd[key]['ls'],
	 		transform=ccrs.PlateCarree(),
	 		zorder=1,
	 		)
 	#####
 	# Add data
	if showdata:	
	 	ax.contour(
	 		newdLon,newdLat,rdregid,
	 		contours,
	 		colors=pd['Data']['c'],
	 		linewidths=pd['Data']['lw'],
	 		linestyles=pd['Data']['ls'],
	 		transform=ccrs.PlateCarree(),
	 		zorder=1,
	 		)
	#####
	# Add land:	 	
	ax.add_feature(cfeature.LAND,  facecolor='white',zorder=2)
	ax.coastlines(lw=0.5,zorder=3)
	pyplot.title(title)


	#####
	# Add legend:
	if addLegend:
		legendSize = len(pd.keys())+1
		ncols = int(legendSize/5)+1
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
		


	#####
	# Saving image:
	print "saving",filename
	pyplot.savefig(filename )		
	pyplot.close()
	

		
	

class extentMaps:
  def __init__(self,
  		modelFiles, 
		dataFile,
		dataType	= '',
		modelcoords 	= '',
		modeldetails 	= '',
		datacoords 	= '',
		datadetails 	= '',								
		datasource	= '',
		model 		= '',
		jobID		= '',
		layers	 	= '',
		regions	 	= '',			
#		metrics	 	= '',
                contours	= '',	
		zrange 		= '',
		workingDir	= '',
		imageDir	= '',						
		grid		= '',
		gridFile	= '',
		debug		= True,
		):
	#####
	#	This is the class that does most of the legwork.
	#	First we save all the initialisation settings as class attributes.
		
	
	if debug: print "timeseriesAnalysis:\t init."	
	self.modelFiles 	= modelFiles 		
	self.dataFile		= dataFile
	self.dataType		= dataType
	self.modelcoords 	= modelcoords		
	self.modeldetails 	= modeldetails
	self.datacoords 	= datacoords
	self.datadetails 	= datadetails						
	self.datasource		= datasource
	self.model 		= model
	self.jobID		= jobID
	self.layers	 	= layers
	self.regions	 	= regions			
	self.contours		= contours
	self.zrange		= zrange	
	#self.metrics	 	= metrics						
	self.grid		= grid
	self.gridFile		= gridFile
	self.workingDir		= workingDir
  	self.imageDir 		= imageDir
	self.debug		= debug
	
		
  	self.shelvefn 		= ukp.folder(self.workingDir)+'_'.join([self.jobID,self.dataType,])+'_contour.shelve'
	self.shelvefn_insitu	= ukp.folder(self.workingDir)+'_'.join([self.jobID,self.dataType,])+'_contour_insitu.shelve'

	#####
	# Run everything
 	self.run()
	

  def run(self,):
  	
  	
 # 	nc = Dataset(self.modelFiles[0],'r')
 	print self.modelFiles
 	print self.dataFile

	
	dataDL = tst.DataLoader(self.dataFile,'',self.datacoords,self.datadetails, regions = self.regions, layers = self.layers,)

	for l in self.layers:		
	    for r in self.regions:
	    	modeldata = {}
	    	modellat = {}
	    	modellon = {}
	    	
	    	
	    	realdata = dataDL.load[(r,l,)]
	    	reallat = dataDL.load[(r,l,'lat')]
	    	reallon = dataDL.load[(r,l,'lon')]
		    		    	
		for mfile in self.modelFiles:
			nc = Dataset(mfile,'r')
			ts = tst.getTimes(nc,self.modelcoords)
			meantime = int(np.mean(ts))
			print "\ttime:",meantime
		
			modelDL = tst.DataLoader(mfile,nc,self.modelcoords,self.modeldetails, regions = self.regions, layers = self.layers,)	
	
	    
			modeldata[meantime] = modelDL.load[(r,l)]
			modellat[meantime] = modelDL.load[(r,l,'lat')]
			modellon[meantime] = modelDL.load[(r,l,'lon')]
			nc.close()
		    		
		filename = ukp.folder(self.imageDir)+'_'.join([self.jobID,self.dataType,l,r])+'.png'	    		
	    	title = ' '.join([getLongName(na) for na in [self.jobID, self.dataType, l, str(meantime) ]])
    		interannualExtendMap(
    			modeldata, modellat, modellon,
			realdata, reallat, reallon,
			self.contours,
			filename,
			title= title,
			labels='',
			zrange = self.zrange,
			)
	return
	#assert 0    				
	

	for mfile in self.modelFiles:
		nc = Dataset(mfile,'r')
		ts = tst.getTimes(nc,self.modelcoords)
		meantime = int(np.mean(ts))
		print "\ttime:",meantime
		
		modelDL = tst.DataLoader(mfile,nc,self.modelcoords,self.modeldetails, regions = self.regions, layers = self.layers,)	
	
		for l in self.layers:		
		    for r in self.regions:
		    
			modeldata = modelDL.load[(r,l)]
			modellat = modelDL.load[(r,l,'lat')]
			modellon = modelDL.load[(r,l,'lon')]
		
		    	realdata = dataDL.load[(r,l,)]
		    	reallat = dataDL.load[(r,l,'lat')]
		    	reallon = dataDL.load[(r,l,'lon')]
		    	
		    		
	    		for mesh in [1,]:
			    	
			    	if mesh: filename = ukp.folder(self.imageDir)+'_'.join([self.jobID,self.dataType,l,r,str(meantime), 'mesh',])+'.png'
			    	else:	 filename = ukp.folder(self.imageDir)+'_'.join([self.jobID,self.dataType,l,r,str(meantime)])+'.png'	    		
			    	
		    		makeExtentPlot(	modeldata, modellat, modellon,
	    				realdata, reallat, reallon,
	    				self.contours,
	    				filename,
	    				title=' '.join([getLongName(na) for na in [self.jobID, self.dataType, l, str(meantime) ]]),
	    				labels='',
	    				zrange = self.zrange,
	    				colourmesh = mesh,
	    				)
  	
  	  	
  	
  	




