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
#


from netCDF4 import Dataset
import UKESMpython as ukp
import numpy as np
from changeNC import changeNC, AutoVivification


""" 	This code makes a mask netcdf for the regions written below.
	this code is needed for profileAnalysis.py
"""

def makeMaskNC(outFile, regions, grid):

	if grid in ['eORCA1',]:
		orcaGridfn 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'	

	#####
	# load mask and coordinates.
	ncmesh = Dataset(orcaGridfn,'r')
	landmask = ncmesh.variables['tmask'][:]
	lats 	 = ncmesh.variables['nav_lat'][:]
	lons 	 = ncmesh.variables['nav_lon'][:]
	depths 	 = ncmesh.variables['nav_lev'][:]
	floattype = ncmesh.variables['e3t_0'][:].dtype
	inttype = ncmesh.variables['umask'][:].dtype
	ncmesh.close()

	landmask = np.ma.masked_where(landmask==0,landmask)

	#####
	# Create Temporary Output Arrays.
	arr 	= []
	arr_lat = []
	arr_lon = []
	arr_t 	= []  	  		
	arr_z 	= []  
	arr_y 	= []  
	arr_x 	= []  

	######
	# Make 1D arrays of coordinates
	print 'Make 1D arrays of coordinates'
	for index,v in ukp.maenumerate(landmask):
		(z,y,x) 	= index 
		t = 0
		
		la = lats[y,x]  			
		lo = lons[y,x]  			
	
		arr.append(v)
		arr_t.append(t)
		arr_z.append(z)
		arr_y.append(y)
		arr_x.append(x)	
		arr_lat.append(la)
		arr_lon.append(lo)

	arr_t 	= np.array(arr_t)
	arr_z	= np.array(arr_z)
	arr_lat	= np.array(arr_lat)
	arr_lon = np.array(arr_lon)
	arr 	= np.array(arr)

	######
	# Calculate 1D mask
	oneDmasks = {}	
	for r in regions:
		print 'Calculate 1D mask',r
		mask = ukp.makeMask('mask name', r, 
	  				arr_t,
	  				arr_z,
	  				arr_lat,
	  				arr_lon,
	  				arr)
		oneDmasks[r] = np.ma.masked_where((mask==1)+ (arr==0), arr )
	######
	# Convert 1D mask to 3D
	threeDmasks={}
	for r in regions:
		print 'Convert 1D mask to 3D',r
		mask = np.zeros_like(landmask)
		for i,m in enumerate(oneDmasks[r]):
			mask[arr_z[i],arr_y[i],arr_x[i]] = m
		threeDmasks[r] = mask

	plotting = 1
	if plotting:
		from matplotlib import pyplot
		for r in regions:	
			pyplot.pcolormesh(threeDmasks[r].sum(0),cmap='jet')
			pyplot.colorbar()
			pyplot.savefig('images/mask_'+r+'.png' )
			print "Saved", r, 'map image'
			pyplot.close()

	av = AutoVivification()
	for r in regions:
		av['newVar'][r]['name']		= r
		av['newVar'][r]['long_name']	= r+ ' mask'
		av['newVar'][r]['units']	= ''
		av['newVar'][r]['newDims']	= ('z','y','x')
		av['newVar'][r]['dtype']	= inttype
		av['newVar'][r]['newData']	= threeDmasks[r][:]	

	removes = [u'e1f', u'e1t', u'e1u', u'e1v', u'e2f', u'e2t', u'e2u', u'e2v', u'ff', u'fmask', u'fmaskutil', u'gdepu', u'gdepv', u'glamf', u'glamt', u'glamu', u'glamv', u'gphif', u'gphit', u'gphiu', u'gphiv', u'isfdraft',  u'misf',   u'tmaskutil',u'umaskutil',  u'vmaskutil', u'e3t', u'e3u', u'e3v', u'e3w', u'e3t_0', u'e3w_0', u'gdept,', u'gdepw', u'gdept_0', u'gdepw_0']
	for rem in removes:
		av[rem]['name']='False'
	
	c = changeNC(orcaGridfn, outfile, av)

def main():
	regions	= ['Global', 
		'ignoreInlandSeas',
	  	'SouthernOcean','ArcticOcean',
		'Equator10', 'Remainder',
		'NorthernSubpolarAtlantic','NorthernSubpolarPacific',
		]
	fn = 	'data/eORCA1_masks.nc'
	grid = 'eORCA'
	makeMaskNC(fn, regions, grid)
	
if __name__=="__main__":
	main()


