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
from netCDF4 import Dataset
from glob import glob



from scipy.interpolate import interp1d
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
		
		MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
		NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	
		if annual:	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/annual/"
		else:		WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
		MAREDATFolder 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"
		GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
		MLDFolder	= "/data/euryale7/scratch/ledm/IFREMER-MLD/"
		workDir		= "/data/euryale7/scratch/ledm/ukesm_postProcessed/"
		imgDir		= ukp.folder('images')
		eORCAgrid 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
		GlodapDir	= "/data/euryale7/backup/ledm/Observations/GLODAP/"
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
		
		
		


	doChl		= 0#True
	doN		= True
	doSi		= True
	doO2		= True
	doAlk		= True
	doDIC		= 0#True
	doAirSeaFlux	= 0#True	
	doIntPP_Lester	= 0#True
	doIntPP_OSU	= 0#True

	doT		= 0#True
	doS		= 0#True
	doMLD		= 0#True
	
	medusaCoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	takahashiCoords	= {'t':'index_t', 'z':'index_z',  'lat': 'LAT', 'lon': 'LON', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	glodapCoords	= {'t':'index_t', 'z':'depth',  'lat': 'latitude', 'lon': 'longitude', 'cal': 'standard','tdict':[] }
	mldCoords	= {'t':'index_t', 'z':'index_z','lat':'lat','lon':'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	#regions 	= ['Global','NorthAtlanticOcean','SouthAtlanticOcean',]
	allRegions	= ['Global','SouthernHemisphere','NorthernHemisphere',
			  'NorthAtlanticOcean','SouthAtlanticOcean','EquatorialAtlanticOcean',
			  'Atlantic','Arctic','nino3','nino3.4','atl_spg','ne_atl','persian']
	keyRegions 	=  ['Global','SouthernHemisphere','NorthernHemisphere',]
	
	av = ukp.AutoVivification()
		
	if doChl:
		name = 'Chlorophyll'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
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
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  WOAFolder+'nitrate_monthly_1deg.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'nitrate', 'vars':['DIN',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'nitrate', 'vars':['n_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		= ['Surface','300m',]#'1000m',]#'Surface - 300m',]'100m',
		av[name]['regions'] 		= allRegions#['Global',]#'NorthAtlanticOcean','SouthAtlanticOcean',]#'NorthAtlantic']
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

	if doSi:

		name = 'Silicate'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  WOAFolder+'silicate_monthly_1deg.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'silicate', 'vars':['SIL',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'silicate', 'vars':['i_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		=  ['Surface','100m','300m','1000m',]
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
		
	
	if doO2:
		name = 'Oxygen'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  WOAFolder+'oxygen-woa13.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'oxygen', 'vars':['OXY',], 'convert': ukp.NoChange,}	
		av[name]['datadetails']  	= {'name': 'oxygen', 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol/m^3'}
	
		av[name]['layers'] 		= ['Surface','100m','300m','1000m',]
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid
			

	if doAlk:
		name = 'Alkalinity'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
		av[name]['dataFile'] 		=  GlodapDir+'Alk.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= glodapCoords
	
		av[name]['modeldetails'] 	= {'name': 'Alkalinity', 'vars':['ALK',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'Alkalinity', 'vars':['Alk',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		=  ['Surface','100m','300m','1000m',]
		av[name]['regions'] 		= keyRegions
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
		av[name]['dataFile'] 		= "/data/euryale7/backup/ledm/Observations/LestersReportData/PPint_1deg.nc"

#	#		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
#	#		av['intpp']['Data']['Vars'] 	= ['PPint',]

				
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
		
	if doT:
		name = 'Temperature'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
		av[name]['dataFile'] 		= WOAFolder+'temperature_monthly_1deg.nc'	
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'temperature', 'vars':['votemper',], 'convert': ukp.NoChange,}
		av[name]['datadetails']  	= {'name': 'temperature', 'vars':['t_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		= ['Surface',]#'Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'NEMO'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

					
	if doS:
		name = 'salinity'
		av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
		av[name]['dataFile'] 		= WOAFolder+'salinity_monthly_1deg.nc'	
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': 'salinity', 'vars':['vosaline',], 'convert': ukp.NoChange,}	
		av[name]['datadetails']  	= {'name': 'salinity', 'vars':['s_an',], 'convert': ukp.NoChange,}
	
		av[name]['layers'] 		= ['Surface',]#'Surface - 1000m','Surface - 300m',]#'depthint']
		av[name]['regions'] 		= keyRegions
		av[name]['metrics']		= ['mean','median', ]

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'NEMO'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= eORCAgrid

	if doMLD:
		
		def mldapplymask(nc,keys):
			mld = nc.variables[keys[0]][:]
			return np.ma.masked_where((np.tile(nc.variables[keys[1]][:],(12,1,1))==0.)+mld.mask+(mld==1.E9),mld)	
		nc = Datasaet(eORCAgrid,'r')
		depth = nc('nav_lev')[:]#
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
	#analysis_timeseries(jobID = "u-ab671")		
	analysis_timeseries(jobID = "u-ab749")			
	
	
	
	
	 
