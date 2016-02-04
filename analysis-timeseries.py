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



import numpy as np


#Specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis#trafficlightsPlots, getMeanSurfaceChl,getHorizontalSlice




		

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

	

		#1. DIN at surface, 100m, 300m and 1000m (vs. "nitrate" in WOA)
		#2. SIL at surface, 100m, 300m and 1000m (vs. "silicate" in WOA)
		#3. OXY at surface, 100m, 300m and 1000m (vs. "dissolved oxygen" in WOA; note the units aren't the same)
		#4. DIC at surface, 100m, 300m and 1000m (vs. pre-industrial version in GLODAP [*])
		#5. ALK at surface, 100m, 300m and 1000m (vs. "total alkalinity" in GLODAP [*])

		#Here I reckon simply reproducing the format (geographical maps of fields and differences from observations) 
		#	and analysis (regional averages and average differences from observations) would be a good start.


		#6. (CHN+CHD) at surface (compare with an appropriate chlorophyll product)
		#7. ((PRN + PRD) *  6.625 * 12.011 / 1e3) for integrated primary production (compare with the OSU data products) 
	#####
	# Location of data files.
	annual = 0
	if gethostname().find('pmpc')>-1:	
		print "analysis-JASMIN.py:\tBeing run at PML on ",gethostname()
		
		MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
		NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	
		if annual:	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/annual/"
		else:		WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
		MAREDATFolder 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"
		GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
		MLDFolder	= "/data/euryale7/scratch/ledm/IFREMER-MLD/"
		workDir		= "/data/euryale7/scratch/ledm/ukesm_postProcessed/"
		imgDir		= ukp.folder('images')
		
	if gethostname().find('ceda.ac.uk')>-1:
		print "analysis-JASMIN.py:\tBeing run at CEDA on ",gethostname()
			
		esmvalFolder = "/group_workspaces/jasmin/esmeval/example_data/bgc/"
		
		#####
		# Location of model files.	
		MEDUSAFolder_pref	= ukp.folder(esmvalFolder+"MEDUSA/")
		NEMOFolder_pref		= ukp.folder(esmvalFolder+"MEDUSA/")
		
		#####
		# Location of data files.
		if annual:	WOAFolder 	= ukp.folder(esmvalFolder+"WOA/annual")
		else:		WOAFolder 	= ukp.folder(esmvalFolder+"WOA/")
		MAREDATFolder 	= ukp.folder(esmvalFolder+"/MAREDAT/")
		WOAFolder 	= ukp.folder(esmvalFolder+"WOA/")
		GEOTRACESFolder = ukp.folder(esmvalFolder+"GEOTRACES/GEOTRACES_PostProccessed/")
		TakahashiFolder = ukp.folder(esmvalFolder+"Takahashi2009_pCO2/")
		MLDFolder  	= ukp.folder(esmvalFolder+"IFREMER-MLD/")
	
		# Directory for output files:
		workDir 	= ukp.folder(esmvalFolder+"ukesm_postProcessed/")
		imgDir		= ukp.folder('images')	
		
		
		


	doChl		= True
	doN		= True
	doSi		= True
	doO2		= True
	doIntPP		= False
	
	medusaCoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	
		
	av = ukp.AutoVivification()	
	if doChl:
		name = 'Chlorophyll'
		av[name]['modelFiles']  	= sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= maredatCoords
	
		av[name]['modeldetails'] 	= {'name': 'Chlorophylla', 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}
		av[name]['datadetails']  	= {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
	
		av[name]['layers'] 		= ['Surface','100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= ['Global',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'MAREDAT'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'


	if doN:
		name = 'Nitrate'
		av[name]['modelFiles']  	= sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  WOAFolder+'nitrate_monthly_1deg.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'nitrate', 'vars':['DIN',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'nitrate', 'vars':['n_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		= ['Surface','100m','200m','Surface - 300m',]
		av[name]['regions'] 		= ['Global',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'

	if doSi:

		name = 'Silicate'
		av[name]['modelFiles']  	= sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  WOAFolder+'silicate_monthly_1deg.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'silicate', 'vars':['SIL',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'silicate', 'vars':['i_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		= ['Surface','100m','200m','Surface - 300m',]
		av[name]['regions'] 		= ['Global',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
		
	
	if doO2:
		name = 'Oxygen'
		av[name]['modelFiles']  	= sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  WOAFolder+'oxygen-woa13.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'oxygen', 'vars':['OXY',], 'convert': ukp.NoChange,}	
		av[name]['datadetails']  	= {'name': 'oxygen', 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol/m^3'}
	
		av[name]['layers'] 		= ['Surface','100m','200m','Surface - 300m',]
		av[name]['regions'] 		= ['Global',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
			
			

	if doIntPP:
		#def NoChange(nc,keys):	return nc.variables[keys[0]][:]		
		
		#def depthInt(nc,keys):
			#sum(nc.variables[keys[0]] * 5.	
		name = 'IntegratedPrimaryProduction'
		av[name]['modelFiles']  	= sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"+"MarEDat20121001Pigments.nc"	
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= maredatCoords
	
		av[name]['modeldetails'] 	= {'name': 'Chlorophylla', 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}
		av[name]['datadetails']  	= {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
	
		av[name]['layers'] 		= ['Surface','100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= ['Global',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'MAREDAT'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
		
		
	shelvedir 	= ukp.folder('shelves/timeseries/'+jobID)
	imagedir	 = ukp.folder('images/timeseries/'+jobID)
	
	shelves = {}
	shelves_insitu={}
	for name in av.keys():
		tsa = timeseriesAnalysis(
			av[name]['modelFiles'], 
			av[name]['dataFile'],
			dataType	= name,
  			modelcoords 	= av[name]['modelcoords'],
  			modeldetails 	= av[name]['modeldetails'],
  			datacoords 	= av[name]['datacoords'],
  			datadetails 	= av[name]['datadetails'],								
			datasource	= av[name]['datasource'],
			model 		= av[name]['model'],
			jobID		= jobID,
			layers	 	= av[name]['layers'],
			regions	 	= av[name]['regions'],			
			metrics	 	= av[name]['metrics'],
			workingDir	= shelvedir,
			imageDir	= imagedir,					
			grid		= av[name]['grid'],
			gridFile	= av[name]['gridFile'],
			clean 		= clean,
		)
		shelves[name] = tsa.shelvefn
		shelves_insitu[name] = tsa.shelvefn_insitu

if __name__=="__main__":	
	analysis_timeseries(jobID = "u-ab671")		
	analysis_timeseries(jobID = "u-ab749")			
	
	
	
	
	 
