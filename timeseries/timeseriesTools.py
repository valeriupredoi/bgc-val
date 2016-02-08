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
import os 

#Specific local code:
import UKESMpython as ukp
from convertToOneDNC import convertToOneDNC



#def getTimeAndData(fn, coords,details, layer='Surface',region = 'Global'):
#	###
#	# mostly unused
#	if type(nc) == type('filename'):
#		nc = Dataset(nc,'r')
#	
#	ts = getTimes(nc,coords)
#	
#	if layer in ['Surface','100m','200m','500m','1000m','2000m',]:
#		data = getHorizontalSlice(nc,coords,details,layer)
#	if layer == 'depthIntergrated':
#		data = getDepthIntegrated(nc,coords,details,layer)
#	if layer == 'Surface - 1000m':
#		data = getHorizontalSlice(nc,coords,details,'Surface')
#		data = data - getHorizontalSlice(nc,coords,details,'1000m',)
#	if region !='Global':
#		print "Not ready for non-global cuts"
#	return ts,data		
#


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


def ApplyDepthSlice(arr,k):
	if arr.ndim == 4: return arr[:,k,:,:]
	if arr.ndim == 3: return arr[k,:,:]
	if arr.ndim == 2: return arr	
	
def getHorizontalSlice(nc,coords,details,layer,data = ''):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	
	if coords['z'] == '' or coords['z'] not in nc.variables.keys():
		print "getHorizontalSlice:\tNo depth field in",details['name']
		if data =='': 	data = ukp.extractData(nc,details)
		return data
				
	if len(nc.variables[coords['z']][:]) ==1 and layer in ['Surface',]:
		print "getHorizontalSlice:\tNo depth field only 1 value",details['name']	
		if data =='': 	data = ukp.extractData(nc,details)
		return ApplyDepthSlice(data, 0)
	
	if layer in ['Surface','100m','200m','300m','500m','1000m','2000m',]:	
	
		if layer == 'Surface':	z = 0.
		if layer == '100m': 	z = 100.			
		if layer == '200m': 	z = 200.
		if layer == '300m': 	z = 300.		
		if layer == '500m': 	z = 500.
		if layer == '1000m': 	z = 1000.
		if layer == '2000m': 	z = 2000.
		k =  ukp.getORCAdepth(z,nc.variables[coords['z']][:],debug=True)
		if data =='': 	data = ukp.extractData(nc,details)
		print "getHorizontalSlice:\tSpecific depth field requested",details['name'], layer,[k],nc.variables[coords['z']][k], data.shape
		return ApplyDepthSlice(data, k)
			
	elif layer in  ['Surface - 1000m', 'Surface - 300m']:
		if layer == 'Surface - 300m':  	z = 300.
		if layer == 'Surface - 1000m': 	z = 1000.
		k_surf =  ukp.getORCAdepth(0., nc.variables[coords['z']][:],debug=True)
		k_low  =  ukp.getORCAdepth(z , nc.variables[coords['z']][:],debug=True)
		print "getHorizontalSlice:\t",layer,"surface:",k_surf,'-->',k_low
		if data =='': 
			return ApplyDepthSlice(ukp.extractData(nc,details),k_surf) - ApplyDepthSlice(ukp.extractData(nc,details),k_low)
		return ApplyDepthSlice(data, k_surf) - ApplyDepthSlice(data, k_low)

	elif layer.lower() == 'depthint':
		assert 0		
#		return getDepthIntegrated(nc,coords,details,layer)
	else:
		assert 0



def applyRegionMask(nc,coords,details, region, layer = '',data = ''):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	if data == '': data = ukp.extractData(nc,details)
	if region in ['Global','All']:return data

	#xt = 	
	#xz = 
	#xy =
	#xx = 
	#xd = 
	m = ukp.makeMask(details['name'],region, xt,xz,xy,xx,xd)
	return np.ma.masked_where(m,data).compressed()
	
	
	
	
class DataLoader:
  def __init__(self,fn,nc,coords,details, regions = ['Global',], layers = ['Surface',],data = ''):
  	self.fn = fn
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')  
  	self.nc 	= nc
  	self.coords 	= coords
  	self.details 	= details
  	self.regions 	= regions
  	self.layers 	= layers
	if data == '': data = ukp.extractData(nc,self.details)  	
  	self.Fulldata 	= data
	
	for r in self.regions:
		if r in ['Global','All']:  continue
		tmpname = ukp.folder('tmp/')+os.path.basename(self.fn).replace('.nc','-'+self.details['name']+'.nc')
		if ukp.shouldIMakeFile(self.fn,tmpname,):
			c = convertToOneDNC(self.fn, tmpname, variables=self.details['vars'], debug=True)
		self.loadXYZArrays(tmpname)
		break

	self.run()
	
  def run(self):
  	self.load = {}
  	for region in self.regions:
  	    for layer in self.layers:  	
  	    	if region in  ['Global','All']:
   			self.load[(region,layer)] =  getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata)
		else:
   			self.load[(region,layer)] =  self.applyRegionMask(region,layer)
		
		
  def loadXYZArrays(self,fn1d):
  	"""	This code loads the lat, lon, depth arrays.
  		For masking.
  	"""
  	nc1d = Dataset(fn1d,'r')
  	
	self.t = np.ma.array(nc1d.variables[self.coords['t']][:])
	self.z = np.ma.array(nc1d.variables[self.coords['z']][:])
	self.z_index = np.ma.array(nc1d.variables['index_z'][:])
	#lat and lon
	self.y = np.ma.array(nc1d.variables[self.coords['lat']][:])
	self.x = np.ma.array(nc1d.variables[self.coords['lon']][:])
	
	self.d = np.ma.array(ukp.extractData(nc1d,self.details)[:])

	self.m = self.d.mask + self.z.mask + self.t.mask + self.y.mask + self.x.mask
  	nc1d.close()
  	return
  
  def applyRegionMask(self,region,layer):
  	"""	This code loads the lat, lon, depth arrays.
  		For masking.
  	""" 
  	mask = self.m
  	
	if layer in ['Surface','100m','200m','300m','500m','1000m','2000m',]:	
		if layer == 'Surface':	z = 0.
		if layer == '100m': 	z = 100.			
		if layer == '200m': 	z = 200.
		if layer == '300m': 	z = 300.		
		if layer == '500m': 	z = 500.
		if layer == '1000m': 	z = 1000.
		if layer == '2000m': 	z = 2000.
		k =  ukp.getORCAdepth(z,self.nc.variables[self.coords['z']][:],debug=True)
		mask += np.ma.masked_where(self.z_index!=k,self.z).mask
		z = np.ma.masked_where(mask,self.z).compressed()
		y = np.ma.masked_where(mask,self.y).compressed()
		x = np.ma.masked_where(mask,self.x).compressed()
		t = np.ma.masked_where(mask,self.t).compressed()
		d = np.ma.masked_where(mask,self.d).compressed()
		m = ukp.makeMask(self.details['name'],region, t,z,y,x,d)
		return np.ma.masked_where(m,d).compressed()
	else:
  		assert 0
  		
  		
  	
  	return 	

def getDepthIntegrated(nc,coords,details,layer,data):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	if data == '': data = ukp.extractData(nc,details)
	
	
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
	
	
	
