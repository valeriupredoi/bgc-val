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

import numpy as np
from netCDF4 import Dataset,num2date

#Specific local code:
import UKESMpython as ukp




def getTimeAndData(fn, coords,details, layer='Surface',region = 'Global'):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	
	ts = getTimes(nc,coords)
	
	if layer in ['Surface','100m','200m','500m','1000m','2000m',]:
		data = getHorizontalSlice(nc,coords,details,layer)
	if layer == 'depthIntergrated':
		data = getDepthIntegrated(nc,coords,details,layer)
	
	if region !='Global':
		print "Not ready for non-global cuts"
	return ts,data		



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
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	dtimes = num2date(nc.variables[coords['t']][:], nc.variables[coords['t']].units,calendar=coords['cal'])[:]
	ts = np.array([float(dt.year) + dt.dayofyr/365. for dt in dtimes])
	return ts

def loadData(nc,details):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	return ukp.extractData(nc,details)[:]
	
def getHorizontalSlice(nc,coords,details,depthLevel,data = ''):
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
	if data =='': return ukp.extractData(nc,details)[:,k,:,:]
	return data[:,k,:,:]

def getDepthIntegrated(nc,coords,details,depthLevel):

	return []



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
	
	
	
