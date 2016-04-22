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


#####	
# Load Standard Python modules:
from sys import argv,exit
from os.path import exists
from calendar import month_name
from socket import gethostname
from netCDF4 import Dataset
from glob import glob
from scipy.interpolate import interp1d
import numpy as np
import os

#####	
# Load specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis




		

def analysis_timeseries(jobID = "u-ab671",
			clean = 0,
			annual = True,
			strictFileCheck = True,
			):

	"""
		The role of this code is to produce time series analysis.
		The jobID is the monsoon/UM job id and it looks for files with a specific format
		
		The clean flag allows you to start the analysis without loading previous data.
		
		The annual flag means that we look at annual (True) or monthly (False) data.
		
		The strictFileCheck switch checks that the data/model netcdf files exist.
			It fails if the switch is on and the files no not exist.
	"""	

	#####
	# Switches:
	# These are some booleans that allow us to choose which analysis to run. 
	# I hope that they have sensible enough names.
	
	#####	
	# BGC switches:
	doChl		= 0#True
	doN		= True
	doSi		= True
	doO2		= True
	doAlk		= True
	doDIC		= True
	doOMZ		= 0#True
	doAirSeaFlux	= 0#True	
	doIntPP_Lester	= 0#True
	doIntPP_OSU	= 0#True
	doExportRatio   = True
	
	#####	
	# Physics switches:
	doT		= True
	doS		= True
	doMLD		= 0#True
		

	#####
	# Location of images directory
	# the imagedir is where the analysis images will be saved.
	imagedir	 = ukp.folder('images/timeseries/'+jobID)
	
	#####
	# Location of shelves folder
	# The shelve directory is where the intermediate processing files are saved in python's shelve format.
	# This allows us to put away a python open to be re-opened later.
	# This means that we can interupt the analysis without loosing lots of data and processing time, 
	# or we can append new simulation years to the end of the analysis without starting from scratch each time.
	shelvedir 	= ukp.folder('shelves/timeseries/'+jobID)

	

	#####
	# Location of data files.
	# The first thing that this function does is to check which machine it is being run. 
	# This is we can run the same code on multiple machines withouht having to make many copies of this file.
	# So far, this has been run on the following machines:
	#	PML
	#	JASMIN
	#	Charybdis (Julien's machine at NOCS)
	#
	# Feel free to add other macihines onto this list, if need be.
	machinelocation = ''
	
	#####
	# PML
	if gethostname().find('pmpc')>-1:	
		print "analysis-timeseries.py:\tBeing run at PML on ",gethostname()
		machinelocation = 'PML'
		MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
		NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	
		if annual:	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/annual/"
		else:		WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
		
		ObsFolder = "/data/euryale7/backup/ledm/Observations/"
		MAREDATFolder 	= ObsFolder+"/MAREDAT/MAREDAT/"
		GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
		MLDFolder	= ObsFolder+"/IFREMER-MLD/"
		LesterFolder	= ObsFolder+"/LestersReportData/"
		GlodapDir	= ObsFolder+"/GLODAP/"
		GLODAPv2Dir	= ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
		
		eORCAgrid 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
				
	#####
	# JASMIN		
	if gethostname().find('ceda.ac.uk')>-1:
		print "analysis-timeseries.py:\tBeing run at CEDA on ",gethostname()
		machinelocation = 'JASMIN'	
				
		ObsFolder 	= "/group_workspaces/jasmin/esmeval/example_data/bgc/"
		esmvalFolder 	= "/group_workspaces/jasmin/esmeval/data/"		
		#####
		# Location of model files.	
		MEDUSAFolder_pref	= ukp.folder(esmvalFolder)
		NEMOFolder_pref		= ukp.folder(esmvalFolder)
		
		#####
		# Location of data files.
		if annual:	WOAFolder 	= ukp.folder(esmvalFolder+"WOA/annual")
		else:		WOAFolder 	= ukp.folder(esmvalFolder+"WOA/")
		#MAREDATFolder 	= ukp.folder(esmvalFolder+"/MAREDAT/")
		#WOAFolder 	= ukp.folder(esmvalFolder+"WOA/")
		#GEOTRACESFolder = ukp.folder(esmvalFolder+"GEOTRACES/GEOTRACES_PostProccessed/")
		#TakahashiFolder = ukp.folder(esmvalFolder+"Takahashi2009_pCO2/")
		#MLDFolder  	= ukp.folder(esmvalFolder+"IFREMER-MLD/")
		
		ObsFolder = esmvalFolder
		MAREDATFolder 	= ObsFolder+"/MAREDAT/MAREDAT/"
		GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
		MLDFolder	= ObsFolder+"/IFREMER-MLD/"
		LesterFolder	= ObsFolder+"/LestersReportData/"
		GlodapDir	= ObsFolder+"/GLODAP/"
		GLODAPv2Dir	= ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
		
		
		eORCAgrid 	= '/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'

	#####
	# NOC		
	if gethostname().find('charybdis')>-1:	
		print "analysis-timeseries.py:\tBeing run at NOC on ",gethostname()
		machinelocation = 'NOC'
		
		MEDUSAFolder_pref	= "/home/jpp1m13/Documents/WORKING/UKESM/Compar_Atm_forcings/netcdf_files/"
		NEMOFolder_pref		="/home/jpp1m13/Documents/WORKING/UKESM/Compar_Atm_forcings/netcdf_files/" 
	
		if annual:	WOAFolder 	= "/home/jpp1m13/Documents/WORKING/UKESM/Compar_Atm_forcings/netcdf_files/"
		else:		WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
		MAREDATFolder 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"
		GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
		MLDFolder	= "/data/euryale7/scratch/ledm/IFREMER-MLD/"
		eORCAgrid 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
		GlodapDir	= "/data/euryale7/backup/ledm/Observations/GLODAP/"		

	#####
	# Unable to find location of files/data.	
	if not machinelocation:
		print "analysis-timeseries.py:\tFATAL:\tWas unable to determine location of host: ",gethostname()
		assert False
		


	#####
	# Coordinate dictionairy
	# These are python dictionairies, one for each data source and model.
	# This is because each data provider seems to use a different set of standard names for dimensions and time.
	# The 'tdict' field is short for "time-dictionary". 
	#	This is a dictionary who's indices are the values on the netcdf time dimension.
	#	The tdict indices point to a month number in python numbering (ie January = 0)
	# 	An example would be, if a netcdf uses the middle day of the month as it's time value:
	#		tdict = {15:0, 45:1 ...}	
	
	
	medusaCoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	takahashiCoords	= {'t':'index_t', 'z':'index_z',  'lat': 'LAT', 'lon': 'LON', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	glodapCoords	= {'t':'index_t', 'z':'depth',  'lat': 'latitude', 'lon': 'longitude', 'cal': 'standard','tdict':[] }
	glodapv2Coords	= {'t':'time',    'z':'Pressure','lat':'lat',      'lon':'lon',        'cal': '',        'tdict':{0:0,} }
	mldCoords	= {'t':'index_t', 'z':'index_z','lat':'lat','lon':'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}



	#####
	# Some lists of region.
	# This are pre-made lists of regions that can be investigated.
	# Note that each analysis below can be given its own set of regions.	
#	regions 	= ['Global','NorthAtlanticOcean','SouthAtlanticOcean',]
#	oldRegions	= ['Global','SouthernHemisphere','NorthernHemisphere',
#			  'NorthAtlanticOcean','SouthAtlanticOcean','EquatorialAtlanticOcean',
#			  'Atlantic','Arctic','nino3','nino3.4','atl_spg','ne_atl','persian']
	#shortRegions 	= ['Global','SouthernHemisphere','NorthernHemisphere',]
  	#AndyRegions 	= ['SouthernOcean','NorthernSubpolarAtlantic','NorthernSubpolarPacific','Arctic','SouthRemainder','NorthRemainder', 'Global']
 
	# need to add:  ,'SouthRemainder','NorthRemainder', 
 	debugRegions	= ['Global','Equator10', 'Remainder','ArcticOcean','NorthernSubpolarAtlantic','NorthernSubpolarPacific','ignoreInlandSeas','SouthernOcean',]
 	
  	allRegions 	= debugRegions	
  	keyRegions 	= debugRegions	  	
  	
	#alllayers = range(40)
	alllayers = [0,2,5,10,15,20,25,30,35,40,45,50,55,60,70,]#80,]
	alllayers.append('Surface')  	
  	
  	
  	#####
  	# The analysis settings:
  	# Below here is a list of analysis settings.
  	# The settings are passed to timeseriesAnalysis using a nested dictionary (called an autovivification, here).
  	# 
  	# These analysis were switched on or off at the start of the function.
  	# Each analysis requires:
  	#	model files
  	#	data files
  	#	model and data coordinate dictionaries, (defines above)
  	#	model and data details (a set of instructions of what to analyse:
  	#		name: 		field name
  	#		vars:		variable names in the netcdf
  	#		convert: 	a function to manipuate the data (ie change units, or add two fields together.
  	#				There are some standard ones in UKESMPython.py, but you can write your own here.
  	#		units: 		the units after the convert function has been applied.
  	#		layers:		which vertical layers to look at (ie, surface, 100m etc...)
  	#		regions:	which regions to look at. Can be speficied here, or use a pre-defined list (from above)
  	#		metrics:	what metric to look at:  mean, median or sum
  	#		model and data source: 	the name of source of the model/data (for plotting)
  	#		model grid: 	the model grid, usually eORCA1
  	#		the model grid file: 	the file path for the model mesh file (contains cell area/volume/masks, etc)
  	#
  	#	Note that the analysis can be run with just the model, it doesn't require a data file.
  	#	If so, just set to data file to an empty string:
  	#		av[name]['dataFile']  = ''
  	
	av = ukp.AutoVivification()
	if doChl:
		name = 'Chlorophyll'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= maredatCoords
	
		av[name]['modeldetails'] 	= {'name': 'Chlorophylla', 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}
		av[name]['datadetails']  	= {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
	
		av[name]['layers'] 		= ['Surface',]#'100m','200m','Surface to 100m', 'Surface to 300m']#'Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= allRegions #['SouthernHemisphere','NorthernHemisphere','Global',	'NorthAtlanticOcean','SouthAtlanticOcean','EquatorialAtlanticOcean']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'MAREDAT'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid


	if doN:
		name = 'Nitrate'
		if annual:
			av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av[name]['dataFile'] 		=  WOAFolder+'/woa13_all_n00_01.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'nitrate', 'vars':['DIN',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'nitrate', 'vars':['n_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions
		
		#av[name]['layers'] 		= ['Surface','300m',]#'1000m',]#'Surface - 300m',]'100m',
		#av[name]['regions'] 		= allRegions#['Global',]#'NorthAtlanticOcean','SouthAtlanticOcean',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

	if doSi:

		name = 'Silicate'
		if annual:	
			av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av[name]['dataFile'] 		=  WOAFolder+'woa13_all_i00_01.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'silicate', 'vars':['SIL',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'silicate', 'vars':['i_an',], 'convert': ukp.NoChange,}
		
		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions
			
		#av[name]['layers'] 		=  ['Surface','100m','300m','1000m',]
		#av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
		
	
	if doO2:
		name = 'Oxygen'
		if annual:
			av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av[name]['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'oxygen', 'vars':['OXY',], 'convert': ukp.NoChange,}	
		av[name]['datadetails']  	= {'name': 'oxygen', 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol/m^3'}

		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions
			
#		av[name]['layers'] 		= ['Surface','100m','300m','1000m',]#
	#	av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid


	if doOMZ:
		# Here we calculate the volume of the OMZ, where the O2 concentration is below 20.
		nc = Dataset(eORCAgrid,'r')

	

		assert 0 
		# not ready yet
		name = 'Oxygen'
		if annual:
			av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av[name]['dataFile'] 		=  WOAFolder+'oxygen-woa13.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'oxygen', 'vars':['OXY',], 'convert': ukp.NoChange,}	
		av[name]['datadetails']  	= {'name': 'oxygen', 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol/m^3'}
	
		av[name]['layers'] 		= ['Surface',] #'100m','300m','1000m',]
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['sum', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

	
	if doDIC:
	
		name = 'DIC'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		= GLODAPv2Dir+'GLODAPv2.tco2.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= glodapv2Coords
	
		av[name]['modeldetails'] 	= {'name': 'DIC', 'vars':['DIC',],  'convert': ukp.NoChange,'units':'mmol-C/m3'}
		av[name]['datadetails']  	= {'name': 'DIC', 'vars':['tco2',], 'convert': ukp.NoChange,'units':'micro-mol kg-1'}
	
		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'GLODAP'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid		

	if doAlk:
		def convertmeqm3TOumolkg(nc,keys):
			return nc.variables[keys[0]][:]* 1.027
		
		name = 'Alkalinity'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  GlodapDir+'Alk.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= glodapCoords
	
		av[name]['modeldetails'] 	= {'name': 'Alkalinity', 'vars':['ALK',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'Alkalinity', 'vars':['Alk',], 'convert': convertmeqm3TOumolkg,}
	
	#	av[name]['layers'] 		=  ['Surface','100m','300m','1000m',]
	#	av[name]['regions'] 		= keyRegions
		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions		
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'GLODAP'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
					


	if doAirSeaFlux:
	
		nc = Dataset(eORCAgrid,'r')
		area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
		nc.close()
		def eOrcaTotal(nc,keys):
			factor = 365. / 1.E9
			arr = nc.variables[keys[0]][:].squeeze()
			if arr.ndim ==3:
				for i in np.arange(arr.shape[0]):
					arr[i] = arr[i]*area
			elif arr.ndim ==2: arr = arr*area
			else: assert 0
			return arr * factor
					
		def takaTotal(nc,keys):
			
			factor = 1.E12
			arr = nc.variables[keys[0]][:].squeeze()
			area = nc.variables[keys[1]][:].squeeze()
			arr = arr*area
			print arr.sum(), arr.sum()*factor
			return arr * factor
			# area 10^6 km^2			
			# flux:  10^15 g Carbon month^-1. (GT)/m2/month

			
		name = 'AirSeaFluxCO2'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_diad_T.nc"))
		av[name]['dataFile'] 		=  TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= takahashiCoords
	
		av[name]['modeldetails'] 	= {'name': 'AirSeaFluxCO2', 'vars':['CO2FLUX',], 'convert': eOrcaTotal,'units':'ton C/yr'}
		av[name]['datadetails']  	= {'name': 'AirSeaFluxCO2', 'vars':['TFLUXSW06','AREA_MKM2'], 'convert': takaTotal,'units':'ton C/yr'}
	
		av[name]['layers'] 		= ['Surface',]
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['sum',]

		av[name]['datasource'] 		= ''
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
					

										
					

	if doIntPP_Lester:
		#def NoChange(nc,keys):	return nc.variables[keys[0]][:]		
		
		def medusadepthInt(nc,keys):
			return (nc.variables[keys[0]][:]+ nc.variables[keys[1]][:])* 6.625 * 12.011 / 1000.	
		name = 'IntegratedPrimaryProduction_1x1'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_diad_T.nc"))
		av[name]['dataFile'] 		= LesterFolder+"/PPint_1deg.nc"

				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= maredatCoords
	
		av[name]['modeldetails'] 	= {'name': 'IntPP', 'vars':['PRN' ,'PRD'], 'convert': medusadepthInt,'units':'mg C/m^3'}
		#av[name]['datadetails']  	= {'name': 'IntPP', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
		av[name]['datadetails']  	= {'name': 'IntPP', 'vars':['PPint',], 'convert': ukp.div1000,'units':'[ug/L/d'}		

	
		av[name]['layers'] 		= ['Surface',]#'100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', 'sum' ]

		av[name]['datasource'] 		= 'MAREDAT'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
		
		
	if doIntPP_OSU:
		nc = Dataset(eORCAgrid,'r')
		area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
		nc.close()
		def medusadepthInt(nc,keys):
			#	 mmolN/m2/d        [mg C /m2/d]   [mgC/m2/yr] [gC/m2/yr]     Gt/m2/yr
			factor = 1.		* 6.625 * 12.011 * 365.	      / 1000.   /     1E15
			arr = (nc.variables[keys[0]][:]+ nc.variables[keys[1]][:]).squeeze()*factor
			if arr.ndim ==3:
				for i in np.arange(arr.shape[0]):
					arr[i] = arr[i]*area
			elif arr.ndim ==2: arr = arr*area
			else: assert 0
			return arr
			
		
				
		name = 'IntegratedPrimaryProduction_OSU'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_diad_T.nc"))
		av[name]['dataFile'] 		= "/data/euryale7/scratch/ledm/OSU/standard_VGPM.SeaWIFS.global.nc"

		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= glodapCoords



		nc = Dataset(av[name]['dataFile'] ,'r')
		lats = nc.variables['latitude'][:]
		osuareas = np.zeros((1080, 2160))
		osuarea = (111100. / 6.)**2. # area of a pixel at equator. in m2
		for a in np.arange(1080):osuareas[a] = np.ones((2160,))*osuarea*np.cos(np.deg2rad(lats[a]))
		
		
		def osuconvert(nc,keys):
			arr = nc.variables[keys[0]][:,:,:] 
			tlen = arr.shape[0]
			
			arr  = arr.sum(0)/tlen * 365.	/ 1000. /     1E15
			if arr.ndim ==3:
				for i in np.arange(arr.shape[0]):
					arr[i] = arr[i]*osuarea
			elif arr.ndim ==2: arr = arr*osuarea
			else: assert 0
			return arr
						
			

	
		av[name]['modeldetails'] 	= {'name': 'IntPP', 'vars':['PRN' ,'PRD'], 'convert': medusadepthInt,'units':'gC/yr'}
		#av[name]['datadetails']  	= {'name': 'IntPP', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
		av[name]['datadetails']  	= {'name': 'IntPP', 'vars':['NPP',], 'convert': osuconvert,'units':'gC/yr'}

	
		av[name]['layers'] 		= ['Surface',]#'100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['sum', ]

		av[name]['datasource'] 		= 'OSU'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid


	if doExportRatio:
		
		def calcExportRatio(nc,keys):
			a = (nc.variables['SDT__100'][:] +nc.variables['FDT__100'][:])/ (nc.variables['PRD'][:] +nc.variables['PRN'][:] )
			a = np.ma.masked_where(a>1.01, a)
			return 	a
			
		name = 'exportRatio'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_diad_T.nc"))
		av[name]['dataFile'] 		= ""

				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= maredatCoords
	
		av[name]['modeldetails'] 	= {'name': 'exportRatio', 'vars':['SDT__100','FDT__100' ,'PRD','PRN',], 'convert': calcExportRatio,'units':''}
		av[name]['datadetails']  	= {'name':'','units':'',}
	
		av[name]['layers'] 		= ['Surface',]#'100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median',]

		av[name]['datasource'] 		= ''
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid		


		
	if doT:
		name = 'Temperature'
		if annual:		
			av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
			av[name]['dataFile'] 		= WOAFolder+'woa13_decav_t00_01v2.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'temperature', 'vars':['votemper',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'temperature', 'vars':['t_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions	
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'NEMO'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

					
	if doS:
		name = 'salinity'
		if annual:
			av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
			av[name]['dataFile'] 		= WOAFolder+'woa13_decav_s00_01v2.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'salinity', 'vars':['vosaline',], 'convert': ukp.NoChange,}	
		av[name]['datadetails']  	= {'name': 'salinity', 'vars':['s_an',], 'convert': ukp.NoChange,}

		av[name]['layers'] 		=  alllayers
		av[name]['regions'] 		= debugRegions		
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'NEMO'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

	if doMLD:
		
		def mldapplymask(nc,keys):
			mld = nc.variables[keys[0]][:]
			return np.ma.masked_where((np.tile(nc.variables[keys[1]][:],(12,1,1))==0.)+mld.mask+(mld==1.E9),mld)	
		nc = Dataset(eORCAgrid,'r')
		depth = nc.variables['nav_lev'][:]#
		nc.close()
		
	#	depth10 = x
	#	ndepth = range(0,50,1)
	#	ndepth.extend(range(50,80,2))
	#	ndepth.extend(range(80,110,3))
	#	ndepth.extend(range(110,200,5))
	#	ndepth.extend(range(200,500,10))
	#	ndepth.extend(range(500,1000,50))
	#	ndepth.extend(range(1000,2000,100))
	#	ndepth.extend(range(2000,6000,200))

		def calcMLD(nc,keys):
			#mlds = np.arange(nc.zeros_like(nc.variables[keys[0]][)

			temp = nc.variables[keys[0]][:,:,:,:]
			f_out = interp1d(depth[7:9],temp[7:9], axis=1)
			tcrit = 0.2
			t10m =  f_out(10.) 
			t10m = np.ma.masked_where(t10m>1E20, t10m) - tcrit
			
			#nc.variables[keys[0]][:,depth10,:,:]
			# linear regression to extrapolate below this level to find the first? 
			f_out = interp1d(temp, depth, axis=1)
			#t_out = f_out(newdepth)

			

			
		name = 'MLD'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
		av[name]['dataFile'] 		= MLDFolder+"mld_DT02_c1m_reg2.0.nc"
			#MLD_DT02 = depth where (T = T_10m +/- 0.2 degC)
 	
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= mldCoords
	
		#av[name]['modeldetails'] 	= {'name': 'mld', 'vars':['somxl010',],   'convert': ukp.NoChange,'units':'m'}	
		av[name]['modeldetails'] 	= {'name': 'mld', 'vars':['votemper',],   'convert': calcMLD,'units':'m'}	
		av[name]['datadetails']  	= {'name': 'mld', 'vars':['mld','mask',], 'convert': mldapplymask,'units':'m'}
	
		av[name]['layers'] 		= ['Surface',]#'Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'IFREMER'
		av[name]['model']		= 'NEMO'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
		
			
				

  	#####
  	# Calling timeseriesAnalysis
	# This is where the above settings is passed to timeseriesAnalysis, for the actual work to begin.
	# We loop over all fiels in the first layer dictionary in the autovificiation, av.
	#	
	# Once the timeseriesAnalysis has completed, we save all the output shelves in a dictionairy.
	# At the moment, this dictioary is not used, but we could for instance open the shelve to highlight specific data,
	#	(ie, andy asked to produce a table showing the final year of data.
	
	shelves = {}
	shelves_insitu={}
	for name in av.keys():
		print "------------------------------------------------------------------"	
		print "analysis-Timeseries.py:\tBeginning to call timeseriesAnalysis for ", name

		if len(av[name]['modelFiles']) == 0:
			print "analysis-Timeseries.py:\tWARNING:\tmodel files are not found:",av[name]['modelFiles'] 
			if strictFileCheck: assert 0		

		modelfilesexists = [os.path.exists(f) for f in av[name]['modelFiles'] ]
		if False in modelfilesexists:
			print "analysis-Timeseries.py:\tWARNING:\tnot model files do not all exist:",av[name]['modelFiles'] 
			if strictFileCheck: assert 0
			
			
		if av[name]['dataFile']!='':
		   if not os.path.exists(av[name]['dataFile']):
			print "analysis-Timeseries.py:\tWARNING:\tdata file is not found:",av[name]['dataFile']
			if strictFileCheck: assert 0
			
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
	#analysis_timeseries(jobID = "u-ab671")		
	analysis_timeseries(jobID = "u-ab749",)#clean=1)			
	#analysis_timeseries(jobID = "u-ab963")			
	
	
	
	
	 
