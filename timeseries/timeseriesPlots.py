#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
#
# This file is part of the ukesm-validation library.
#
# ukesm-validation is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license. 

# ukesm-validation is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with ukesm-validation.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk

from matplotlib import pyplot
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import cartopy
import numpy as np 

import timeseriesTools as tst

def trafficlights(ax,xlims, bands,labels=[]):
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
		pyplot.legend(handles=handles)
	return ax
	
	
def trafficlightsPlots(
		times, 			# model times (in floats)
		arr,			# model time series
		dataslice,		# in situ data distribution
		title 	='',
		filename='',
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
	
	print "UKESMpython:\tscatterPlot:\tSaving:" , filename
	pyplot.savefig(filename )
	pyplot.close()	


def trafficlightsPlot(
		times, 			# model times (in floats)
		arr,			# model time series
		dataslice,		# in situ data distribution
		metric  = '',
		title 	='',
		filename='',
		units = '',
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
		pc1 = np.array([np.percentile(dataslice,20.) for i in xlims]) 
		pc2 = np.array([np.percentile(dataslice,30.) for i in xlims])
		pc3 = np.array([np.percentile(dataslice,45.) for i in xlims])
		pc4 = np.array([np.percentile(dataslice,60.) for i in xlims])
		pc5 = np.array([np.percentile(dataslice,70.) for i in xlims])
		pc6 = np.array([np.percentile(dataslice,80.) for i in xlims])
		labels = ['20-30 pc','30-40 pc','40-60 pc','60-70 pc','70-80 pc',]
		pcs = [pc1,pc2,pc3,pc4,pc5,pc6]
		ax = trafficlights(ax,xlims, pcs ,labels=labels)
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
			im = ax.scatter(lons, lats,c=data, lw=0,marker='s', transform=cartopy.crs.PlateCarree(),norm=LogNorm())#vmin=zrange[0],vmax=zrange[1]),)
		else:	
			im = ax.scatter(lons, lats,c=data, lw=0,marker='s',transform=cartopy.crs.PlateCarree(),)#vmin=zrange[0],vmax=zrange[1])	
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

	ax2 = pyplot.subplot(212,projection=cartopy.crs.PlateCarree(central_longitude=0.0, ))	
	fig,ax2 = makemapplot(fig,ax2,lons2,lats2,data2,titles[1], zrange=[rbmi,rbma],lon0=0.,drawCbar=True,cbarlabel='',doLog=doLog,)
	
		
	print "timeseriespots.py:\tmapPlotPair: \tSaving:" , filename
	pyplot.savefig(filename ,dpi=dpi)		
	pyplot.close()
		
		
