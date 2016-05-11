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

from matplotlib import pyplot, gridspec

from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import cartopy
import numpy as np 
from scipy import interpolate 

import timeseriesTools as tst

try:	defcmap = pyplot.cm.viridis
except:	
	from bgcvaltools.viridis import viridis
	defcmap = viridis
	
	

def trafficlights(ax,xlims, bands,labels=[],drawlegend=True):
	if len(bands) != 6: print "wrong number of bands to traffic light."

	ax.fill_between(xlims,bands[0], bands[1] ,color='r', 		alpha = 0.2)		
	ax.fill_between(xlims,bands[1] ,bands[2] ,color='DarkOrange', 	alpha = 0.2)				
	ax.fill_between(xlims,bands[2] ,bands[3] ,color='g', 		alpha = 0.2)
	ax.fill_between(xlims,bands[3] ,bands[4] ,color='DarkOrange', 	alpha = 0.2)
	ax.fill_between(xlims,bands[4] ,bands[5] ,color='r', 		alpha = 0.2)	
	
	if len(labels)==5:
		patch4 = mpatches.Patch(color='r', 		alpha = 0.2,label=labels[4])
		patch3 = mpatches.Patch(color='DarkOrange', 	alpha = 0.2,label=labels[3])
		patch2 = mpatches.Patch(color='g', 		alpha = 0.2,label=labels[2])
		patch1 = mpatches.Patch(color='DarkOrange', 	alpha = 0.2,label=labels[1])
		patch0 = mpatches.Patch(color='r', 		alpha = 0.2,label=labels[0])		
		handles = [patch4,patch3,patch2,patch1,patch0,]
		if drawlegend: pyplot.legend(handles=handles)
	return ax

def drawgreyband(ax,xaxis, bands,labels=[]):
	ax.fill_between(xaxis,bands[0], bands[1] ,color='k', 	alpha = 0.05)
	return ax

	
	
def trafficlightsPlots(
		times, 			# model times (in floats)
		arr,			# model time series
		dataslice,		# in situ data distribution
		title 	='',
		filename='',
		greyband = False,
	):

	xlims= [times[0],times[-1]]
	
	fig = pyplot.figure()
	
	ax = fig.add_subplot(211)	
	pyplot.plot(times,arr)
	pyplot.xlim(xlims)	
	pyplot.title(title)	

	if len(dataslice):
		#pyplot.axhline(y=np.ma.mean(dataslice),c='k',ls='-',lw=2,alpha=0.5)
		pyplot.axhline(y=np.ma.median(dataslice),c='k',ls='-',lw=1,)#alpha=0.5)	
		pc1 = np.array([np.percentile(dataslice,25.) for i in xlims]) 
		pc2 = np.array([np.percentile(dataslice,35.) for i in xlims])
		pc3 = np.array([np.percentile(dataslice,45.) for i in xlims])
		pc4 = np.array([np.percentile(dataslice,55.) for i in xlims])
		pc5 = np.array([np.percentile(dataslice,65.) for i in xlims])
		pc6 = np.array([np.percentile(dataslice,75.) for i in xlims])
		labels = ['25-35 pc','35-45 pc','45-55 pc','55-65 pc','65-75 pc',]
		ax = trafficlights(ax,xlims, [pc1,pc2,pc3,pc4,pc5,pc6],labels=labels)
		
		pcmin 	= np.array([dataslice.min() for i in xlims]) 
		pcmax 	= np.array([dataslice.min() for i in xlims]) 	
		if greyband:
			ax  	= drawgreyband(ax,xlims, [pcmin,pc1],)
			ax  	= drawgreyband(ax,xlims, [pc6,pcmax],)
				
	
	ax = fig.add_subplot(212)
	newt,cusum = tst.calcCuSum(times,arr)
	
	if len(dataslice):
		pyplot.plot(newt,cusum/np.ma.median(dataslice))
		pyplot.title('Cumulative sum / in situ media')		
	else:
		pyplot.plot(newt,cusum)		
		pyplot.title('Cumulative sum')
		
	pyplot.xlim(xlims)	
	pyplot.axhline(y=0.,c='k',ls='-',lw=2,alpha=0.5)
	
	if len(dataslice):
		m25 = np.array([-0.25 for i in xlims]) 
		m15 = np.array([-0.15 for i in xlims])
		m05 = np.array([-0.05 for i in xlims])
		p05 = np.array([0.05  for i in xlims])		
		p15 = np.array([0.15 for i in xlims])
		p25 = np.array([0.25 for i in xlims]) 
		ax = trafficlights(ax,xlims, [m25,m15,m05,p05,p15,p25])
	
	print "UKESMpython:\ttrafficlightsPlots:\tSaving:" , filename
	pyplot.savefig(filename )
	pyplot.close()	


def percentilesPlot(
		timesDict,		# model times dict
		modeldataDict,		# model data dictionairy
		dataslice,		# in situ data distribution
		title 	='',
		filename='',
		units = '',		
		greyband= True,
		dolog=True
	):
	
	metrics = sorted(timesDict.keys())

	#####
	# Determine the x axis and y axis corners.	
	mint = np.ma.min(timesDict[metrics[0]])
	maxt = np.ma.max(timesDict[metrics[0]])
	miny = np.ma.min(modeldataDict['min'])
	maxy = np.ma.max(modeldataDict['max'])
	for m in metrics:
		if np.ma.min(timesDict[m]) < mint: mint = np.ma.min(timesDict[m])
		if np.ma.max(timesDict[m]) > maxt: maxt = np.ma.max(timesDict[m])
		if np.ma.min(modeldataDict[m]) < miny: miny = np.ma.min(modeldataDict[m])
		if np.ma.max(modeldataDict[m]) > maxy: maxy = np.ma.max(modeldataDict[m])
	if len(dataslice):	
		try:	dataslice = dataslice.compressed()
		except:	pass
		if np.ma.min(dataslice) < miny: miny = np.ma.min(dataslice)
		if np.ma.max(dataslice) > maxy: maxy = np.ma.max(dataslice)			
	xlims= [mint,maxt]
	ylims= [miny,maxy]

	if dolog and ylims[0] *ylims[-1] <=0.:
		print "zero (or less) in ylims", ylims,
		smallestNonZeroModel = getSmallestAboveZero(modeldataDict['min'])
		smallestNonZeroData = getSmallestAboveZero(dataslice)
		ylims[0] = getSmallestAboveZero([smallestNonZeroModel, smallestNonZeroData])
		print "new ylims:",ylims
		#assert 0
		
	if np.ma.max(ylims)/np.ma.min(ylims)>20.:	
		#####
		# log
		ylims[0] = ylims[0]*0.6
		ylims[1] = ylims[1]*1.4	
	else:	#####
		# not log
		diff = abs(ylims[1] - ylims[0])
		ylims[0] = ylims[0]-diff/20.
		ylims[1] = ylims[1]+diff/20.	
			
	#####
	# Make, resize and divide the figure into an uneven grid/
	fig = pyplot.figure()
	fig.set_size_inches(10,6)
	gs = gridspec.GridSpec(1, 2, width_ratios=[1, 12], wspace=0.0, hspace=0.0) 		
	
	#####
	# Data plot to the left.
	axd = pyplot.subplot(gs[0])	
	if len(dataslice):

		
		pc1 = np.clip(np.array([np.percentile(dataslice,20.) for i in xlims]),ylims[0],ylims[1])
		pc2 = np.clip(np.array([np.percentile(dataslice,30.) for i in xlims]),ylims[0],ylims[1])
		pc3 = np.clip(np.array([np.percentile(dataslice,40.) for i in xlims]),ylims[0],ylims[1])
		pc4 = np.clip(np.array([np.percentile(dataslice,60.) for i in xlims]),ylims[0],ylims[1])
		pc5 = np.clip(np.array([np.percentile(dataslice,70.) for i in xlims]),ylims[0],ylims[1])
		pc6 = np.clip(np.array([np.percentile(dataslice,80.) for i in xlims]),ylims[0],ylims[1])
		labels = ['20-30 pc','30-40 pc','40-60 pc','60-70 pc','70-80 pc',]
		axd = trafficlights(axd,xlims, [pc1,pc2,pc3,pc4,pc5,pc6],labels=labels,drawlegend=False)

		axd.axhline(y=np.ma.mean(dataslice),  c='k',ls='--',lw=1.5,)
		axd.axhline(y=np.ma.median(dataslice),c='k',ls='-', lw=1.5,)#alpha=0.5)	
				
		pcmin 	= np.clip(np.array([dataslice.min() for i in xlims]),ylims[0],ylims[1])
		pcmax 	= np.clip(np.array([dataslice.max() for i in xlims]),ylims[0],ylims[1])
		if greyband:
			axd  	= drawgreyband(axd,xlims, [pcmin,pc1],)
			axd  	= drawgreyband(axd,xlims, [pc6,pcmax],)
		pyplot.title('Data')	
		pyplot.ylabel(units)
		
			
	#####
	# Model plot to the right.	
	axm = pyplot.subplot(gs[1])	

	pyplot.plot(timesDict[m],np.clip(modeldataDict['mean'],  ylims[0],ylims[1]),'k--',lw=1,label='mean')
	pyplot.plot(timesDict[m],np.clip(modeldataDict['median'],ylims[0],ylims[1]),'k-', lw=1,label='median')	
			
	legend = pyplot.legend(loc='lower center',  numpoints = 1, ncol=2, prop={'size':12}) 
	legend.draw_frame(False) 
	legend.get_frame().set_alpha(0.)
		
	
	labels = ['20-30 pc','30-40 pc','40-60 pc','60-70 pc','70-80 pc',]	
	trafficarray = [np.clip(modeldataDict[m],ylims[0],ylims[1]) for m in ['20pc','30pc','40pc','60pc','70pc','80pc']]
	axm = trafficlights(axm,timesDict['min'], trafficarray,labels=labels,drawlegend=False)		

	if greyband:
		axm  	= drawgreyband(axm,timesDict['min'], [np.clip(modeldataDict[m],ylims[0],ylims[1]) for m in ['min','20pc']],)
		axm  	= drawgreyband(axm,timesDict['min'], [np.clip(modeldataDict[m],ylims[0],ylims[1]) for m in ['80pc','max']],)				


	#yticks = list(axm.get_xticks())
	#ytickLabels = list(axm.get_yticklabels())

	
	axm.set_xlim(xlims)	
	axm.set_title(title)	


					
	axm.set_ylim(ylims)	
	axd.set_ylim(ylims)

	if np.ma.max(ylims)/np.ma.min(ylims)>20.:
		axm.set_yscale('log')	
		axd.set_yscale('log')

	axm.get_yaxis().set_ticklabels([])			
	axd.get_xaxis().set_ticks([])	
	
		
	ytickLabels	= ['min', '20pc','30pc','40pc','60pc','70pc','80pc','max']
	yticks = np.clip([ modeldataDict[m][-1]	for m in ytickLabels],ylims[0],ylims[1])
	axm.set_yticks(yticks)
	axm.set_yticklabels(ytickLabels,fontsize=10)
	axm.yaxis.tick_right()		

	
	print "UKESMpython:\tpercentilesPlot:\tSaving:" , filename
	pyplot.savefig(filename )
	pyplot.close()	


def getSmallestAboveZero(arr):
	arr = np.array(arr)
	return np.ma.masked_where(arr<=0.,arr).min()
	



def trafficlightsPlot(
		times, 			# model times (in floats)
		arr,			# model time series
		dataslice,		# in situ data distribution
		metric  = '',
		title 	='',
		filename='',
		units = '',
		greyband = False		
	):

	if len(times) ==0 or len(arr) == 0:
		print "trafficlightsPlot:\tWARNING:\tdata or time arrays are empty.",len(times),len(arr),title
		return
	if np.ma.is_masked(arr):
		print "trafficlightsPlot:\tWARNING:\tdata arrays is masked",len(times),len(arr),title
		return
				
	xlims= [times[0],times[-1]]
	
	fig = pyplot.figure()
	
	ax = fig.add_subplot(111)	
	pyplot.plot(times,arr)
	pyplot.xlim(xlims)	
	pyplot.title(title)	

	if len(dataslice) and metric != 'sum':
		
		#pyplot.axhline(y=np.ma.mean(dataslice),c='k',ls='-',lw=2,alpha=0.5)
		pyplot.axhline(y=np.ma.median(dataslice),c='k',ls='-',lw=1,)#alpha=0.5)	
		pcmin 	= np.array([dataslice.min() for i in xlims]) 
		pc1 = np.array([np.percentile(dataslice,20.) for i in xlims]) 
		pc2 = np.array([np.percentile(dataslice,30.) for i in xlims])
		pc3 = np.array([np.percentile(dataslice,45.) for i in xlims])
		pc4 = np.array([np.percentile(dataslice,60.) for i in xlims])
		pc5 = np.array([np.percentile(dataslice,70.) for i in xlims])
		pc6 = np.array([np.percentile(dataslice,80.) for i in xlims])
		pcmax 	= np.array([dataslice.max() for i in xlims])
						
		labels = ['20-30 pc','30-40 pc','40-60 pc','60-70 pc','70-80 pc',]
		pcs = [pc1,pc2,pc3,pc4,pc5,pc6]
		ax = trafficlights(ax,xlims, pcs ,labels=labels)
		if greyband:
			ax  	= drawgreyband(ax,xlims, [pcmin,pc1],)
			ax  	= drawgreyband(ax,xlims, [pc6,pcmax],)
				
	if len(dataslice) and metric == 'sum':
		
		pyplot.axhline(y=np.ma.sum(dataslice),c='k',ls='-',lw=1,label ='data'+str(np.ma.sum(dataslice)))#alpha=0.5)		
		pyplot.legend()
	
	print "UKESMpython:\tscatterPlot:\tSaving:" , filename
	pyplot.savefig(filename )
	pyplot.close()	


def regrid(data,lat,lon):
    	nX = np.arange(-179.5,180.5,0.25)
    	nY = np.arange( -89.5, 90.5,0.25)
    	if lat.ndim ==1:
    		oldLon, oldLat = np.meshgrid(lon,lat) 
    	else:
    		   oldLon, oldLat = lon,lat		
    	newLon, newLat = np.meshgrid(nX,nY)
    	
    	crojp1 = cartopy.crs.PlateCarree(central_longitude=0.0, ) #central_latitude=300.0)
    	crojp2 = cartopy.crs.PlateCarree(central_longitude=0.0, ) #central_latitude=300.0)

    	a = cartopy.img_transform.regrid(data,
    			     source_x_coords=oldLon,
                             source_y_coords=oldLat,
                             source_cs=crojp1,
                             target_proj=crojp2,
                             target_x_points=newLon,
                             target_y_points=newLat
                             )
       # print 'newregid shape:',a.shape                     
	return crojp2, a, newLon,newLat
	
	
	
def makemapplot(fig,ax,lons,lats,data,title, zrange=[-100,100],lon0=0.,drawCbar=True,cbarlabel='',doLog=False,):

	
	lons = np.array(lons)
	lats = np.array(lats)
	data = np.ma.array(data)	
	
	
	
	if doLog and zrange[0]*zrange[1] <=0.:
		print "makemapplot: \tMasking"
		data = np.ma.masked_less_equal(ma.array(data), 0.)
	
	print data.min(),lats.min(),lons.min(), data.shape,lats.shape,lons.shape
	
	if data.ndim ==1:
		if doLog:
			im = ax.scatter(lons, lats,c=data, lw=0,marker='s', transform=cartopy.crs.PlateCarree(),norm=LogNorm(),vmin=zrange[0],vmax=zrange[1])
		else:	
			im = ax.scatter(lons, lats,c=data, lw=0,marker='s',transform=cartopy.crs.PlateCarree(),vmin=zrange[0],vmax=zrange[1])
	else:
		crojp2, data, newLon,newLat = regrid(data,lats,lons)

		if doLog:
			im = ax.pcolormesh(newLon, newLat,data, transform=cartopy.crs.PlateCarree(),norm=LogNorm(vmin=zrange[0],vmax=zrange[1]),)
		else:	
			im = ax.pcolormesh(newLon, newLat,data, transform=cartopy.crs.PlateCarree(),vmin=zrange[0],vmax=zrange[1])
	
	ax.add_feature(cartopy.feature.LAND,  facecolor='0.85')	


	
	if drawCbar:
	    c1 = fig.colorbar(im,pad=0.05,shrink=0.75)
	    if len(cbarlabel)>0: c1.set_label(cbarlabel)

	pyplot.title(title)


	ax.set_axis_off()
	pyplot.axis('off')
	ax.axis('off')
		
	
		
	return fig, ax
	
def mapPlotSingle(lons1, lats1, data1,filename,titles=['',],lon0=0.,drawCbar=True,cbarlabel='',doLog=False,dpi=100,):#**kwargs):

	fig = pyplot.figure()
	fig.set_size_inches(10,6)
	
	lons1 = np.array(lons1)
	lats1 = np.array(lats1)
	data1 = np.ma.array(data1)
	
	rbmi = data1.min()
	rbma = data1.max()
	
	if rbmi * rbma >0. and rbma/rbmi > 100.: doLog=True
	ax1 = pyplot.subplot(111,projection=cartopy.crs.PlateCarree(central_longitude=0.0, ))
		
	fig,ax1 = makemapplot(fig,ax1,lons1,lats1,data1,titles[0], zrange=[rbmi,rbma],lon0=0.,drawCbar=True,cbarlabel='',doLog=doLog,)
	ax1.set_extent([-180.,180.,-90.,90.])
	print "mapPlotSingle.py:\tSaving:" , filename
	pyplot.savefig(filename ,dpi=dpi)		
	pyplot.close()
	

def mapPlotPair(lons1, lats1, data1,lons2,lats2,data2,filename,titles=['',''],lon0=0.,drawCbar=True,cbarlabel='',doLog=False,dpi=100,):#**kwargs):

	fig = pyplot.figure()

	fig.set_size_inches(10,10)

	
	lons1 = np.array(lons1)
	lats1 = np.array(lats1)
	data1 = np.ma.array(data1)
	
	lons2 = np.array(lons2)
	lats2 = np.array(lats2)
	data2 = np.ma.array(data2)
	
	rbmi = min([data1.min(),data2.min()])
	rbma = max([data1.max(),data2.max()])		
	
	if rbmi * rbma >0. and rbma/rbmi > 100.: doLog=True
	ax1 = pyplot.subplot(211,projection=cartopy.crs.PlateCarree(central_longitude=0.0, ))
		
	fig,ax1 = makemapplot(fig,ax1,lons1,lats1,data1,titles[0], zrange=[rbmi,rbma],lon0=0.,drawCbar=True,cbarlabel='',doLog=doLog,)
	ax1.set_extent([-180.,180.,-90.,90.])	

	ax2 = pyplot.subplot(212,projection=cartopy.crs.PlateCarree(central_longitude=0.0, ))	
	try:fig,ax2 = makemapplot(fig,ax2,lons2,lats2,data2,titles[1], zrange=[rbmi,rbma],lon0=0.,drawCbar=True,cbarlabel='',doLog=doLog,)
	except: 
		mapPlotSingle(lons1, lats1, data1,filename,titles=titles,lon0=lon0,drawCbar=drawCbar,cbarlabel=cbarlabel,doLog=doLog,dpi=dpi)
		return
	ax2.set_extent([-180.,180.,-90.,90.])	
		
	print "mapPlotPair: \tSaving:" , filename
	pyplot.savefig(filename ,dpi=dpi)		
	pyplot.close()
		

def hovmoellerAxis(fig,ax,title,xaxis,yaxis,data,vmin='',vmax='',cmap = defcmap ,debug = False):
	yaxis = np.array(yaxis)
	if yaxis.min()*yaxis.max() <=0.:
		if yaxis.mean()<0:yaxis = np.clip(yaxis,-10000.,-0.1)
	if debug:	print "hovmoellerAxis:\txaxis:",xaxis,"\tyaxis:",yaxis,"\tdata:",data
	if debug:	print "hovmoellerAxis:\txaxis:",xaxis.shape,"\tyaxis:",yaxis.shape,"\tdata:",data.shape	
	
	if vmin==vmax=='':	p=pyplot.pcolormesh(xaxis,yaxis,data,cmap = cmap)
	else:			p=pyplot.pcolormesh(xaxis,yaxis,data,vmin=vmin,vmax=vmax,cmap = cmap)
	pyplot.title(title)
	#ax.set_yscale("log", nonposy='clip')
	ax.set_yscale('symlog')
	
def zaxisfromCC(arr):
	"""
		This function creates a 1D array for the depth axis.
		It assumes that the z axis given is cell centered (cc) depth, it adds a zero to the start
		and takes the mid point between the cc depths.
		It then converts the array into negative values.
	"""
	arr = np.array(arr)
	zarr = list((arr[1:]+arr[:-1])/2.)
	zarr.insert(0,0.)
	zarr.append(arr[-1])	
	#print arr,zarr
	return -1.*np.abs(np.array(zarr))
		
def taxisfromCC(arr):
	"""
		This function creates a 1D array for the time axis.
		It assumes that the x axis given is cell centered (cc), 
		it adds a starting time to the start and takes the mid point between the cc depths.
		It then converts the array into negative values.
	"""
	arr = np.array(arr)
	diff = np.mean(np.abs(arr[1:]-arr[:-1]))/2.
	zarr = list((arr[1:]+arr[:-1])/2.)
	zarr.insert(0,arr[0]-diff)
	zarr.append(arr[-1]+diff)	
	#print arr,zarr
	return np.array(zarr)
	

def hovmoellerPlot(modeldata,dataslice,filename, modelZcoords = {}, dataZcoords= {}, title = '',dpi=100,diff=True):	
	#####
	# creating model data dictionairies
	md = []
	times_cc = []
	yaxis_cc = []
	
	for l in modelZcoords:
		if l not in modeldata.keys():continue
		yaxis_cc.append(modelZcoords[l])

		times_cc = sorted(modeldata[l].keys())
		try:	md.append(np.array([modeldata[l][t][0] for t in times_cc]))
		except:	md.append(np.array([modeldata[l][t] for t in times_cc]))
		
	md = np.ma.array(md)#
	md = md.squeeze()
	md = np.ma.masked_where(np.ma.masked_invalid(md).mask + md.mask, md)
	times = taxisfromCC(np.array(times_cc))
	yaxis = zaxisfromCC(yaxis_cc)
	print "hovmoellerPlot:", title, md.shape,md.mean(),times.shape,yaxis.shape #, times, yaxis
	
	#####
	# creating data data dictionairies
	dd = []
	dxaxis= []
	dyaxis_cc = []
	
	for l in dataZcoords:
		if l not in dataslice.keys():continue
		dyaxis_cc.append(dataZcoords[l])
		print 'preparing data for hov:',l,dataZcoords[l], dataslice[l]
		dd.append([np.ma.array(dataslice[l]),])
	dd = np.ma.array(dd)#.squeeze()
	dd = np.ma.masked_where(np.ma.masked_invalid(dd).mask + dd.mask, dd)	
	dyaxis_cc = np.array(dyaxis_cc)
		
	
	#####
	# A hack to defuse an unusal bug that occurs when all layers of the data set are masked.
	while len(dd.shape)> 2:
		if dd.shape[-1] == 1:
			dd = dd[:,:,0]
		if len(dd.shape) >2 and dd.shape[-1] !=1: 
			print "Something very strange is happenning with this array:", dd
			assert False
		

	dxaxis = np.array([0,1,])
	dyaxis = zaxisfromCC(dyaxis_cc)
	print "hovmoellerPlot: - data:", title, dd.shape,dd.mean(),dxaxis.shape,dyaxis.shape
	
 	if len(dd.squeeze().compressed())!=0 and diff: 
		print "hovmoellerPlot:\tWARNNG: Data entire masked, can not take difference of entirely masked data."
		return
	
	
	######
	# Locate min/max colours
	rbmi = min([md.min(),dd.min(),])
	rbma = max([md.max(),dd.max(),])
	
	######
	# Locate min/max depths
	zmi = min([dyaxis.min(),yaxis.min(),])
	zma = max([dyaxis.max(),yaxis.max(),])	
	
	

	if diff:
		#####
		# diff:
		# If the diff key is true, this subtracts the data from the model, it also changes the colour scheme.
		
		if len(dd.squeeze().compressed())!=0:
			print "hovmoellerPlot:\tWARNNG: Data entire masked, can not take difference of entirely masked data."
			return
		# First step is to perform the interpollartion.		
		f = interpolate.interp1d(dyaxis_cc, dd.squeeze(), kind='linear')
		newData = f(np.ma.clip(yaxis_cc,dyaxis_cc.min(),dyaxis_cc.max()))	

		###
		# subtract data from model:		
		for i,t in enumerate(times_cc):
			
			md[:,i] = md[:,i] -  newData
	
		####
		# change plot ranges, title, colorscale
		ax2max = max([abs(md.max()),abs(md.min()),])
		ax2min = - ax2max
		cmapax2 = pyplot.cm.bwr
		title = 'Model - Data: '+title		
	else:
		ax2max	= rbma
		ax2min	= rbmi
		cmapax2 = defcmap
		title = 'Model: '+title
	#####
	# Start drawing
	# Grid spec allows you to make un-even shaped subplots.
	# here we want the in situ data to be much narrower.
	fig = pyplot.figure()	
	fig.set_size_inches(10,6)
	gs = gridspec.GridSpec(1, 2, width_ratios=[1, 12]) 

	#####
	# Data  subplot
	
	ax1 = pyplot.subplot(gs[0])
	if len(dd.squeeze().compressed())!=0:
		hovmoellerAxis(fig,ax1,'Data',dxaxis,dyaxis,dd,vmin=rbmi,vmax=rbma)
		if diff: pyplot.colorbar()		
	pyplot.ylim([zmi,zma])
	ax1.get_xaxis().set_ticks([])
	pyplot.ylabel('Depth')


	#####
	# model subplot
	ax2 = pyplot.subplot(gs[1])
	hovmoellerAxis(fig,ax2,title,times,yaxis,md,vmin=ax2min,vmax=ax2max,cmap = cmapax2)
	pyplot.xlim([times.min(),times.max()])
	pyplot.ylim([zmi,zma])	
	pyplot.colorbar()
	ax2.yaxis.set_ticklabels([])
	pyplot.xlabel('Year')

	pyplot.tight_layout()			
	print "hovmoellerPlot.py: \tSaving:" , filename
	pyplot.savefig(filename ,dpi=dpi)		
	pyplot.close()	
	
	
	
	
			
