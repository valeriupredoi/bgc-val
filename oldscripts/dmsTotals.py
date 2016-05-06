#!/usr/bin/ipython 
from matplotlib import pyplot
from calendar import month_name
from glob import glob
from shelve import open as shOpen
from bgcvaltools.StatsDiagram import rmsds
from UKESMpython import folder, makeMask
from itertools import product
import os
import numpy as np
from ncdfView import ncdfView
from pyproj import Proj
from shapely.geometry import shape
from changeNC import changeNC,AutoVivification

purple = [125./256., 38./256., 205./256.]

oceans = ['SouthPacificOcean',  'ArcticOcean',
	  'AntarcticOcean',   'NorthAtlanticOcean','SouthAtlanticOcean',
	  'NorthPacificOcean','IndianOcean',] 
months = [month_name[i] for i in xrange(1,13)]	  

seasons = ['JFM','AMJ','JAS','OND',]

dmsmetrics  = ['dms_and','dms_ara','dms_hal','dms_sim']
dmspmetrics = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']	

###
#	The strategy of this code:
#		load paired netcdfs.
#		apply mask
#		calculate volume
#		save calculated volume as a shelve.
#		make plot from volume. like dmsPatternStats.py
#	

def calculateArea():
	# makes a netcdf with the areas of each pixel in the flat map.
	nc = ncdfView('data/Flat1deg.nc',Quiet=True)
	lat = nc('lat')[:]
	lon = nc('lon')[:]	 
	areas = np.zeros_like(nc('nav_lat')[:])
	
	for i,l in enumerate(lat): #(u'lat', u'lon')
	
		co = {"type": "Polygon", "coordinates": [
		    [(1., l+0.5), #('lon', 'lat')
		     (1., l-0.5),
		     (0., l-0.5),
		     (0., l+0.5)]]}
		clon, clat = zip(*co['coordinates'][0])

		pa = Proj("+proj=aea +lat_1="+str(l-0.5)+" +lat_2="+str(l+0.5)+"+lat_0="+str(l)+" +lon_0=0.5")		
		x, y = pa(clon, clat)
		cop = {"type": "Polygon", "coordinates": [zip(x, y)]}

		area = shape(cop).area  
		print i, l , area
	
		for j,lo in enumerate(lon):
			areas[i,j]=area
	av = AutoVivification()
	av['newVar']['area']['name'] ='area'
	av['newVar']['area']['newDims']	= (u'lat', u'lon')
	av['newVar']['area']['dtype'] = nc.variables['lat'].dtype
	av['newVar']['area']['long_name'] = 'area'
	av['newVar']['area']['units'] = 'm^2'
	av['newVar']['area']['newData'] = areas
	
	c = changeNC('data/Flat1deg.nc', 'data/Flat1deg-Areas.nc', av)
	
	
def run():
	model = 'dms_p_sim'
	fold = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-iMarNet-2526/'+model+'/'
	ncDfn = fold+'Data_dms_p_sim_Surface_MEDUSA-iMarNet-2526_1D.nc'
	ncMfn = fold+'Model_'+model+'_Surface_MEDUSA-iMarNet-2526_1D.nc'
	
	ncGfn = 'data/Flat1deg-Areas.nc'


def main():
	#kys = [o + 'Months' for o in oceans]
	#kys.extend([o + 'Seasons' for o in oceans])
	#kys = [o + 'Seasons' for o in oceans]

	#kys.extend(['seasons','OceansSeasons',])# 'OceanMonths','months', 'oceans',
	#for xkeys in kys: #[,]:
	#    for dmsmodels in ['dmspmetrics','dmsmetrics', ]:
	run()

if __name__=="__main__":
	calculateArea()
	#main()		
		
	print 'The end.'



