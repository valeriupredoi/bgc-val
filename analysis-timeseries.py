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
#

#Standard Python modules:
from sys import argv,exit
from os.path import exists
from calendar import month_name
from socket import gethostname

from glob import glob
from shelve import open as shOpen
from netCDF4 import Dataset,num2date
from matplotlib import pyplot
import numpy as np


#Specific local code:
import UKESMpython as ukp


def getMeanSurfaceChl(fn,coords,details):
	# this needs to be combined with the tools from elsewhere.
	nc = Dataset(fn,'r')
	ts = getTimes(nc,coords)
	d = getHorizontalSlice(nc,coords,details,depthLevel='Surface')
	
	return ts.mean(), d.mean()
	
def getMedianSurfaceChl(fn,coords,details):
	# this needs to be combined with the tools from elsewhere.
	nc = Dataset(fn,'r')
	ts = getTimes(nc,coords)
	d = getHorizontalSlice(nc,coords,details,depthLevel='Surface')
	
	return ts.mean(), np.ma.median(d)	
	
def getTimes(nc, coords):
	dtimes = num2date(nc.variables[coords['t']][:], nc.variables[coords['t']].units,calendar=coords['cal'])[:]
	ts = np.array([float(dt.year) + dt.dayofyr/365. for dt in dtimes])
	return ts
	
def getHorizontalSlice(nc,coords,details,depthLevel):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
		
	if depthLevel in ['Surface','100m','200m','500m','1000m','2000m',]:	
		if depthLevel == 'Surface':	z = 0.
		if depthLevel == '100m': 	z = 100.			
		if depthLevel == '200m': 	z = 200.
		if depthLevel == '500m': 	z = 500.
		if depthLevel == '1000m': 	z = 1000.
		if depthLevel == '2000m': 	z = 2000.
	else:
		assert 0
	k =  ukp.getORCAdepth(z,nc.variables[coords['z']][:],debug=True)
	d = ukp.extractData(nc,details)[:,k,:,:]
	return d



def calcCuSum(times,arr):
	newt,cusum = [],[]
	
	c = 0.
	for i,ti in enumerate(times):
		if i==0:continue
		t0 = times[i-1]
		t= t0 + (ti-t0)/2.
		c += arr[i] - arr[i-1]
		newt.append(t)
		cusum.append(c)

	return newt,cusum
	
		

def analysis_timeseries(jobID = "u-ab671",
			clean = 0,
			):

	"""
		The role of this code is to produce time series analysis.
		
	"""	

	
	#doIntPP = True
	#if doIntPP:
	# code plan:
	# for each metric.

	#	for each model file:
	#	for r in regions:	
	#		iterate over a list of analysis	
	#		extract mean, depth int, etc for each region.
	#		do it monthly, then annual.
	#		Save it into a 
	#	load data
	#	for r in regions:
	#		iterate over a list of analysis
	#		produce a distribution of values.
	#		do it monthly, then annual.	
	#		

	#jobID = "u-ab671"
	files = sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
	
	shelvefn = ukp.folder('shelves/timeseries/'+jobID)+'shelvefn.shelve'
	
	
	modelCoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	dataCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	modelDetails = 	{'name': 'CHL', 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}
	datadetails  =  {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
	
	
			
	####
	# load and calculate the model info
	try:
		if clean: 
			print "User requested clean run. Wiping old data."
			assert 0		
		sh = shOpen(shelvefn)
		readFiles 	= sh['readFiles']
		data 		= sh['data']
		sh.close()
		print "Opened shelve:", shelvefn, '\tread', len(readFiles)
	except:
		readFiles = []
		data = {}
		print "Could not open shelve:", shelvefn, '\tread', len(readFiles)	

	openedFiles = 0					
	for fn in files:
		if fn in readFiles:continue
		print "loading ",fn
		#func = getMeanSurfaceChl
			
		t,m = getMeanSurfaceChl(fn,modelCoords,modelDetails)
		data[t] = m
		readFiles.append(fn)
		openedFiles+=1
	
	if openedFiles:
		print "Saving shelve:", shelvefn, '\tread', len(readFiles)				
		sh = shOpen(shelvefn)
		sh['readFiles']	= readFiles
		sh['data'] 	= data
		sh.close()
	
	# load and calculate data info.
	
	MAREDATFn 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"+"MarEDat20121001Pigments.nc"	
	dataslice = getHorizontalSlice(MAREDATFn,dataCoords,datadetails,depthLevel='Surface')

	med = np.ma.median(dataslice)
	print "data mean:",dataslice.mean(),"data std:",dataslice.std(), 'median:',med
		
	filename = ukp.folder('images/timeseries')+jobID+'meansurfacechl.png'
	

		
	times = sorted(data.keys())
	arr = [data[t] for t in times]
	xlims= [times[0],times[-1]]
	
	fig = pyplot.figure()
	ax = fig.add_subplot(211)	
	pyplot.plot(times,arr)

	pyplot.axhline(y=med,c='k',ls='-',lw=1,)#alpha=0.5)

#	ax.fill_between(xlims,med-dataslice.std() ,med+dataslice.std() ,color='g', alpha = 0.2)
		
	
	pyplot.title('Times series')
	
	pyplot.xlim(xlims)	
	
	
	ax = fig.add_subplot(212)
	newt,cusum = calcCuSum(times,arr)
	
	pyplot.plot(newt,cusum/np.mean(arr))
	pyplot.title('normalised Cumulative sum')
	pyplot.xlim(xlims)	
	pyplot.axhline(y=0.,c='k',ls='-',lw=2,alpha=0.5)
	

	lll = np.array([-0.25 for i in xlims]) 
	l40 = np.array([-0.15 for i in xlims])
	l50 = np.array([-0.05 for i in xlims])
	l60 = np.array([0.05  for i in xlims])		
	l70 = np.array([0.15 for i in xlims])
	lul = np.array([0.25 for i in xlims]) 
	ax.fill_between(xlims,lll, l40 ,color='r', alpha = 0.2)		
	ax.fill_between(xlims,l40 ,l50 ,color='DarkOrange', alpha = 0.2)				
	ax.fill_between(xlims,l50 ,l60 ,color='g', alpha = 0.2)
	ax.fill_between(xlims,l60 ,l70 ,color='DarkOrange', alpha = 0.2)
	ax.fill_between(xlims,l70 ,lul ,color='r', alpha = 0.2)		
	
		
	print "UKESMpython:\tscatterPlot:\tSaving:" , filename
	pyplot.savefig(filename )
	pyplot.close()	
	


if __name__=="__main__":	
	analysis_timeseries(jobID = "u-ab671")		
	analysis_timeseries(jobID = "u-ab749")			
	
	
	
	
	 
