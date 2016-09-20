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
	if arr.ndim == 1: return arr	
	return arr
	
def ApplyDepthrange(arr,k1,k2):
	if arr.ndim == 4: return arr[:,k1:k2,:,:]
	if arr.ndim == 3: return arr[k1:k2,:,:]
	if arr.ndim == 2: return arr
		
def getHorizontalSlice(nc,coords,details,layer,data = ''):
	if type(nc) == type('filename'):
		nc = Dataset(nc,'r')
	
	#####
	# In the case that there is no depth field provided, or the depth field is not the netcdf.
	# We just attempt to extract the data. 
	# This is useful
	if coords['z'] == '' or coords['z'] not in nc.variables.keys():
		print "getHorizontalSlice:\tNo depth field in",details['name']
		if data =='': 	data = ukp.extractData(nc,details)
		return data
				
	####
	#
	if len(nc.variables[coords['z']][:]) ==1 and layer in ['Surface',]:
		print "getHorizontalSlice:\tNo depth field only 1 value",details['name']	
		if data =='': 	data = ukp.extractData(nc,details)
		return ApplyDepthSlice(data, 0)
	if layer in ['layerless',]:
		print "getHorizontalSlice:\tNo layer data requested", layer
		if data =='': 	data = ukp.extractData(nc,details)
		return data

	#####
	# This is when there is only one dimension in the data/model file.
	if len(nc.dimensions.keys()) == 1 and layer in ['Surface','layerless']:
		print "getHorizontalSlice:\tOne D file",details['name']	
		if data =='': 	data = ukp.extractData(nc,details)
		data = np.ma.masked_where(nc.variables[coords['z']][:]>0,data)
		return data
		#return ApplyDepthSlice(data, 0)
		
	
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
	elif layer == 'depthint':
		print "getHorizontalSlice\t:ERROR:\tDepth in should be done in the extractData phase by passing a function in the details dictionary."
		assert 0

	if type(layer) in [type(0),np.int64,np.int,]:
		k = layer
		try:	z = nc.variables[coords['z']][k]
		except:	return []
		if data =='': 	data = ukp.extractData(nc,details)				
		print "getHorizontalSlice:\tSpecific depth level requested",details['name'], layer,nc.variables[coords['z']][k], data.shape	
		return ApplyDepthSlice(data, k)
			
	if layer in nc.variables[coords['z']][:]:
		z = layer
		k =  ukp.getORCAdepth(z,nc.variables[coords['z']][:],debug=False)
		if data =='': 	data = ukp.extractData(nc,details)		
		print "getHorizontalSlice:\tSpecific depth requested",details['name'], layer,nc.variables[coords['z']][k], data.shape	
		return ApplyDepthSlice(data, k)
			
	print "getHorizontalSlice\t:ERROR:\tunrecoginised layer instructions",layer, coords,type(layer)
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
  	self.name	= self.details['name']
	if data == '': data = ukp.extractData(nc,self.details)
  	self.Fulldata 	= data
  	self.__lay__ 	= -999.
	self.run()
	
  def run(self):
  	self.load = {}
   	try:	depth = {i:z for i,z in enumerate(self.nc.variables[self.coords[z]][:])} 
   	except: depths = {}
   	
   	 	
    	for layer in self.layers: 
  	    if type(layer) == type(1):
                
 		if layer not in depths.keys():
                    for region in self.regions:
			a = np.ma.array([-999,],mask=[True,])
		
   			self.load[(region,layer)] =  a 
   			self.load[(region,layer,'t')] =  a
	   		self.load[(region,layer,'z')] =  a
   			self.load[(region,layer,'lat')] =  a 
   			self.load[(region,layer,'lon')] =  a
   	            return	   		
   			   			   				
  	    for region in self.regions:
  	    
  	    	#if region in  ['Global','All']:
  	    	#	dat = self.__getlayerDat__(layer)
  	    		
   		#	self.load[(region,layer)] =  dat
   		#	print "Loading Global Data", (region,layer), dat.min(),dat.max(),dat.mean()
   		#	#getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata)
		#else:
		arr, arr_t, arr_z, arr_lat, arr_lon 	= self.createDataArray(region,layer)
   		self.load[(region,layer)] =  arr 
   		self.load[(region,layer,'t')] =  arr_t
   		self.load[(region,layer,'z')] =  arr_z
   		self.load[(region,layer,'lat')] =  arr_lat 
   		self.load[(region,layer,'lon')] =  arr_lon
   			   			   			   			
   			
  		print "DataLoader:\tLoaded",self.name, 'in',
  		print '{:<24} layer: {:<8}'.format(region,layer),
  		print '\tdata length:',len(self.load[(region,layer)]), '\tmean:',np.ma.mean(self.load[(region,layer)])
  		
  def __getlayerDat__(self,layer):
  	""" Minimise quick load and save to minimise disk-reading time.
  	"""

  	if self.__lay__ == layer:
  		return self.__layDat__
  	else:
  		 self.__layDat__ = np.ma.array(getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata))
  		 print "DataLoader:\tgetlayerDat:",self.name,layer
  		 self.__lay__ = layer
  		 return self.__layDat__
  		 
  	
  def createDataArray(self,region,layer):
  	"""	
  		This creates a set of 1D arrays of the dat and 4D coordinates for the required region.
  		The leg work is done in UKESMpython.py: makeMask()
  	"""
  	
  	#print 'DataLoader:\tcreateDataArray:\t',self.details['name'],region,layer
  	
  	self.createOneDDataArray(layer)
  	  	
  	m = ukp.makeMask(self.details['name'],region, 
  				self.oneDData['arr_t'],
  				self.oneDData['arr_z'],
  				self.oneDData['arr_lat'],
  				self.oneDData['arr_lon'],
  				self.oneDData['arr'],)

  	return 	np.ma.masked_where(m,self.oneDData['arr']),\
  		np.ma.masked_where(m,self.oneDData['arr_t']),\
  		np.ma.masked_where(m,self.oneDData['arr_z']),\
  		np.ma.masked_where(m,self.oneDData['arr_lat']),\
  		np.ma.masked_where(m,self.oneDData['arr_lon'])
  		  		  		  		
  	
  def createOneDDataArray(self,layer):
  	""" 	This is a relatively simple routine that takes a layer and makes a series of 1D arrays containing points.
  		These output 1D arrays are then passed to UKESMpython.py's makemasks toolkit.
  	"""

	#####
	# load lat, lon and data.
  	lat = self.nc.variables[self.coords['lat']][:]
  	lon = ukp.makeLonSafeArr(self.nc.variables[self.coords['lon']][:]) # makes sure it's between +/-180

	dims =   self.nc.variables[self.details['vars'][0]].dimensions
  	
  	dat = self.__getlayerDat__(layer)
	if dat.ndim == 2:	dat =dat[None,:,:]

	try:	l = len(dat)
	except:
		dat = np.ma.array([dat,])
		l = len(dat)
		print "createOneDDataArray: \tWarning:\tdata was a single float:",self.name,dat, l, dat.shape,dat.ndim
	if l == 0:
		a = np.ma.array([-999,],mask=[True,])
		return a,a,a,a,a

		
	#####
	# Create Temporary Output Arrays.
  	arr 	= []
  	arr_lat = []
  	arr_lon = []
  	arr_t 	= []  	  		
  	arr_z 	= []  	


	latnames = ['lat','latitude','latbnd','nav_lat','y',u'lat',]
	lonnames = ['lon','longitude','lonbnd','nav_lon','x',u'lon',]
	

	####
	# Different data has differnt shapes, and order or dimensions, this takes those differences into account.
	if dat.ndim > 2:
	  if dims[-2].lower() in latnames and dims[-1].lower() in lonnames:
	  
  	    #print 'createDataArray',self.details['name'],layer, "Sensible dimsions order:",dims		
 	    for index,v in ukp.maenumerate(dat):

  			try:	(t,z,y,x) 	= index
  			except: 
  				(t,y,x) 	= index 
  				z = 0
 			
  			try:	la = lat[y,x]  			
  			except:	la = lat[y]
  			try:	lo = lon[y,x]  			
  			except:	lo = lon[x]  			
  			
  			arr.append(v)
  			arr_t.append(t)
  			arr_z.append(z)
  			arr_lat.append(la)
  			arr_lon.append(lo)
  			  			  			  			
	  elif dims[-2].lower() in lonnames and dims[-1].lower() in latnames:
	  
  	    #print 'createDataArray',self.details['name'],layer, "Ridiculous dimsions order:",dims				
 	    for index,v in ukp.maenumerate(dat):
  			try:	(t,z,x,y) 	= index
  			except: 
  				(t,x,y) 	= index 
				z = 0	  				
 			
  			try:	la = lat[y,x]  			
  			except:	la = lat[y]
  			try:	lo = lon[y,x]  			
  			except:	lo = lon[x]  			
  			
  			arr.append(v)
  			arr_t.append(t)
  			arr_z.append(z)
  			arr_lat.append(la)
  			arr_lon.append(lo)
  			  			
  	  else:
  		print "Unknown dimensions order", dims
  		assert False	
  		  			
  	elif dat.ndim == 1:
   	  if dims[0] == 'index':
  	    print 'createDataArray',self.name,layer, "1 D data:",dims,dat.shape
 	    for i,v in enumerate(dat):
  			la = lat[i]  			
  			lo = lon[i]  			
  			
  			arr.append(v)  	
  			arr_t.append(0)
  			arr_z.append(0)
  			arr_lat.append(la)
  			arr_lon.append(lo)
  	  elif len(dat)==1:
  	    	print 'createDataArray',self.name,layer, "single point data:",dims,dat.shape
  		arr = dat  	
  		arr_t = [0,]
  		arr_z = [0.,]
  		arr_lat = [0.,]
  		arr_lon = [0.,]
  	    		  	
  	  else:
  		print "Unknown dimensions order", dims
  		assert False	  	
  	else:
  		print "Unknown dimensions order", dims
  		assert False
  			

  	arr = np.ma.masked_invalid(np.ma.array(arr))
  	mask = np.ma.masked_where((arr>1E20) + arr.mask,arr).mask
  	
  	self.oneDData={}
  	self.oneDData['arr_lat'] = np.ma.masked_where(mask,arr_lat).compressed()
  	self.oneDData['arr_lon'] = np.ma.masked_where(mask,arr_lon).compressed()
  	self.oneDData['arr_z']   = np.ma.masked_where(mask,arr_z  ).compressed() 	
  	self.oneDData['arr_t']   = np.ma.masked_where(mask,arr_t  ).compressed()
  	self.oneDData['arr']     = np.ma.masked_where(mask,arr    ).compressed()

	
  	#print 'createDataArray:',arr.min(),arr.mean(),arr.max(),arr, len(arr)  	
  	#print 'createDataArray:',region, arr_lat.min(),arr_lat.mean(),arr_lat.max(),arr_lat, len(arr_lat)  	  	
  	#print 'createDataArray:',region, arr_lon.min(),arr_lon.mean(),arr_lon.max(),arr_lon, len(arr_lon)  	  	  	
  	#return arr, arr_t,arr_z,arr_lat,arr_lon
  	
  	
  	
  
 		
 		
  	




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
		if r in ['Global','All',]: 
			print "Loading global file 1D. "
			assert 0 
			continue
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
  	    	#if region in  ['Global','All']:
   		#	self.load[(region,layer)] =  getHorizontalSlice(self.nc,self.coords,self.details,layer,data = self.Fulldata)
		#else:
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
	imask = np.ma.masked_invalid(self.d).mask
	
	self.m = self.d.mask + self.z.mask + self.t.mask + self.y.mask + self.x.mask + imask
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
  		
  		
