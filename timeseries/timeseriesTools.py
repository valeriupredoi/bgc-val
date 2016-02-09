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
	
def ApplyDepthrange(arr,k1,k2):
	if arr.ndim == 4: return arr[:,k1:k2,:,:]
	if arr.ndim == 3: return arr[k1:k2,:,:]
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
		k =  ukp.getORCAdepth(z,nc.variables[coords['z']][:],debug=False)
		if data =='': 	data = ukp.extractData(nc,details)
		print "getHorizontalSlice:\tSpecific depth field requested",details['name'], layer,[k],nc.variables[coords['z']][k], data.shape
		return ApplyDepthSlice(data, k)
			
	elif layer in  ['Surface - 1000m', 'Surface - 300m']:
		if layer == 'Surface - 300m':  	z = 300.
		if layer == 'Surface - 1000m': 	z = 1000.
		k_surf =  ukp.getORCAdepth(0., nc.variables[coords['z']][:],debug=False)
		k_low  =  ukp.getORCAdepth(z , nc.variables[coords['z']][:],debug=False)
		print "getHorizontalSlice:\t",layer,"surface:",k_surf,'-->',k_low
		if data =='': 
			return ApplyDepthSlice(ukp.extractData(nc,details),k_surf) - ApplyDepthSlice(ukp.extractData(nc,details),k_low)
		return ApplyDepthSlice(data, k_surf) - ApplyDepthSlice(data, k_low)
		
	elif layer in  ['Surface to 100m', 'Surface to 300m']:
		if layer == 'Surface to 300m':  z = 300.
		if layer == 'Surface to 100m': 	z = 100.
		k_surf =  ukp.getORCAdepth(0., nc.variables[coords['z']][:],debug=False)
		k_low  =  ukp.getORCAdepth(z , nc.variables[coords['z']][:],debug=False)
		print "getHorizontalSlice:\t",layer,"surface:",k_surf,'-->',k_low
		if data =='': 
			return ApplyDepthrange(ukp.extractData(nc,details),k_surf,k_low)
		return ApplyDepthrange(data, k_surf,k_low)
	elif layer.lower() == 'depthint':
		assert 0		
#		return getDepthIntegrated(nc,coords,details,layer)
	else:
		assert 0




	
class DataLoader:
  def __init__(self,fn,nc,coords,details, regions = ['Global',], layers = ['Surface',],data = ''):
  	self.fn = fn
	if type(nc) == type('filename'):
		nc = Dataset(fn,'r')  
  	self.nc 	= nc
  	self.coords 	= coords
  	self.details 	= details
  	self.regions 	= regions
  	self.layers 	= layers
	if data == '': data = ukp.extractData(nc,self.details)
  	self.Fulldata 	= data
  	self.__lay__ 	= ''
	self.run()
	
  def run(self):
  	self.load = {}
    	for layer in self.layers:  	
  	    for region in self.regions:
  	    	if region in  ['Global','All']:
  	    		dat = self.__getlayerDat__(layer)
  	    		
   			self.load[(region,layer)] =  dat
   			print "Loading Global Data", (region,layer), dat.min(),dat.max(),dat.mean()
   			#getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata)
		else:
   			self.load[(region,layer)] =  self.createDataArray(region,layer)
   			
  def __getlayerDat__(self,layer):
  	""" Minimise quick load and save to minuimise disk-reading time.
  	"""
  	if self.__lay__ == layer:
  		return self.__layDat__
  	else:
  		 self.__layDat__ = np.ma.array(getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata))
  		 self.__lay__ = layer
  		 return self.__layDat__
  		 
  	
  def createDataArray(self,region,layer):
  	print 'createDataArray',self.details['name'],region,layer
  	
  	lat = self.nc.variables[self.coords['lat']][:]
  	lon = ukp.makeLonSafeArr(self.nc.variables[self.coords['lon']][:]) # makes sure it's between +/-180
  	
  	dat = self.__getlayerDat__(layer)
  	#dat = np.ma.array(getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata))
  	
  	arr = []
  	if region == 'SouthernHemisphere':	regionlims  = {'lat_min':-90.,'lat_max': 0.,'lon_min':-360.,'lon_max':360.}
  	if region == 'NorthernHemisphere':	regionlims  = {'lat_min':  0.,'lat_max':90.,'lon_min':-360.,'lon_max':360.}  	
  	
	if region == 'NorthAtlanticOcean': 	regionlims  = {'lat_min':  10.,'lat_max': 60.,'lon_min':-80.,'lon_max':10.}
	if region == 'SouthAtlanticOcean': 	regionlims  = {'lat_min': -50.,'lat_max':-10.,'lon_min':-65.,'lon_max':20.}	

	if region == 'EquatorialAtlanticOcean': regionlims  = {'lat_min': -15.,'lat_max': 15.,'lon_min':-65.,'lon_max':20.}	
	
	if dat.ndim == 2:	dat =dat[None,:,:]
	
	dims =   self.nc.variables[self.details['vars'][0]].dimensions
  	print 'createDataArray',self.details['name'],region,layer, dims	
	
	if dims[-1].lower() in ['lon','longitude','lonbnd','nav_lon','x'] and dims[-2].lower() in ['lat','latitude','latbnd','nav_lat','y']:
	    # dims order is [t,z,y,x] or [t,y,x] or [z,y,x]
  	    print 'createDataArray',self.details['name'],region,layer, "Sensible dimsions order:",dims		
 	    for index,v in ukp.maenumerate(dat):
  			try:	(t,z,y,x) 	= index
  			except: (t,y,x) 	= index 
 			
  			try:	la = lat[y,x]  			
  			except:	la = lat[y]
  			
  			if la<regionlims['lat_min']:continue
  			if la>regionlims['lat_max']:continue  	
  			
  			try:	lo = lon[y,x]  			
  			except:	lo = lon[x]  			
  			
  			if lo<regionlims['lon_min']:continue
  			if lo>regionlims['lon_max']:continue
  			
  			#if np.ma.is_masked(v):continue  
  			#if v > 1E20:continue    						
  			arr.append(v)

	elif dims[-2].lower() in ['lon','longitude','lonbnd','nav_lon','x'] and dims[-1].lower() in ['lat','latitude','latbnd','nav_lat','y']:
  	    print 'createDataArray',self.details['name'],region,layer, "Ridiculous dimsions order:",dims				
 	    for index,v in ukp.maenumerate(dat):
  			try:	(t,z,x,y) 	= index
  			except: (t,x,y) 	= index 
 			
  			try:	la = lat[y,x]  			
  			except:	la = lat[y]
  			
  			if la<regionlims['lat_min']:continue
  			if la>regionlims['lat_max']:continue  	
  			
  			try:	lo = lon[y,x]  			
  			except:	lo = lon[x]  			
  			
  			if lo<regionlims['lon_min']:continue
  			if lo>regionlims['lon_max']:continue
  			
  			#if np.ma.is_masked(v):continue  
  			#if v > 1E20:continue    						
  			arr.append(v)
  	else:
  		print "Unknown dimensions order", dims
  		assert False	
  			
  	arr = np.ma.masked_invalid(np.ma.array(arr))
  	arr = np.ma.masked_where((arr>1E20) + arr.mask,arr).compressed()
	#print arr.min(), arr.max(),arr.mean(),arr.shape  	
  	#if arr.mean()>40. or arr.max()>100.:
  	#	print arr.min(), arr.max(),arr.mean(),arr.shape
  	#	assert 0
  	return arr
  	
  	
  	
  
 		
 		
  	

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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
class DataLoader_1Darrays:
  def __init__(self,fn,nc,coords,details, regions = ['Global',], layers = ['Surface',],data = ''):
  	assert False
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
		k =  ukp.getORCAdepth(z,self.nc.variables[self.coords['z']][:],debug=False)
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
  		
  		
