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
import os,sys
from getpass import getuser

#####	
# Load specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis
from timeseries import profileAnalysis
from timeseries import timeseriesPlots as tsp 
from pftnames import getLongName


def timeseries_compare():
	### strategy here is a simple wrapper.
	# It's a little cheat-y, as I'm copying straight from analysis_timeseries.py
	
	
	#jobs = ['u-af123', 'u-af139','u-af420','u-af421', ]
	jobs = ['u-af725', 'u-af728','u-af420','u-af421', 'u-af586','u-af730','u-ae748']	
	colours = {'u-af725':'red', 'u-af728':'orange','u-af420':'blue','u-af421':'purple', 'u-af586':'green','u-af730':'pink','u-ae748':'black'}
	annual = True
	strictFileCheck = False

	analysisKeys = []
	#                        analysisKeys.append('N')                        # WOA Nitrate
	#                        analysisKeys.append('Si')                       # WOA Siliate
	#                        analysisKeys.append('O2')                       # WOA Oxygen
	#                        analysisKeys.append('Alk')                      # Glodap Alkalinity
	#                        analysisKeys.append('DIC')                      # Globap tCO2
	#                        analysisKeys.append('AirSeaFlux')               # work in progress

	#                        analysisKeys.append('PP_OSU')                   # OSU Integrated primpary production                    
	#                        analysisKeys.append('LocalExportRatio')         # Export ratio (no data)
	#                        analysisKeys.append('GlobalExportRatio')        # Export ratio (no data)
	#                        analysisKeys.append('TotalOMZVolume')           # Total OMZ Volume
	#                        analysisKeys.append('OMZThickness')             # Total OMZ Volume
	#                        analysisKeys.append('Iron')                     # Iron

		                #####   
		                # Physics switches:
	#                        analysisKeys.append('T')                        # WOA Temperature
	#                        analysisKeys.append('S')                        # WOA Salinity
	#			analysisKeys.append('MLD')			# iFERMER Mixed Layer Depth - work in prgress                        
	


	analysisKeys.append('DrakePassageTransport')	# DrakePassageTransport		
	analysisKeys.append('TotalAirSeaFlux')          # work in progress              
	analysisKeys.append('IntPP_OSU')                # OSU Integrated primpary production    
	analysisKeys.append('GlobalExportRatio')
	analysisKeys.append('TotalIceArea')		# TotalIceArea	
	analysisKeys.append('NorthernTotalIceArea')	# North TotalIceArea
	analysisKeys.append('SouthernTotalIceArea')	# South TotalIceArea

	analysisKeys.append('N')                        # WOA Nitrate
	analysisKeys.append('Si')                       # WOA Siliate
	analysisKeys.append('O2')                       # WOA Oxygen
	analysisKeys.append('Iron')
	analysisKeys.append('Alk')	
	analysisKeys.append('DIC')
	
	layerList 	= ['Surface',]
	metricList 	= ['mean',]
  	regionList	= ['Global',]
	
	#####
	# JASMIN		
	if gethostname().find('ceda.ac.uk')>-1:
		print "analysis-timeseries.py:\tBeing run at CEDA on ",gethostname()
		machinelocation = 'JASMIN'	
				
		ObsFolder 	= "/group_workspaces/jasmin/esmeval/example_data/bgc/"
		esmvalFolder 	= "/group_workspaces/jasmin2/ukesm/BGC_data/"

			
		#####
		# Location of model files.	
		MEDUSAFolder_pref	= ukp.folder(esmvalFolder)
		NEMOFolder_pref		= ukp.folder(esmvalFolder)
		
		#####
		# Location of data files.
		if annual:	WOAFolder 	= ukp.folder(ObsFolder+"WOA/annual")
		else:		WOAFolder 	= ukp.folder(ObsFolder+"WOA/")
		
		MAREDATFolder 	= ObsFolder+"/MAREDAT/MAREDAT/"
		GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
		MLDFolder	= ObsFolder+"/IFREMER-MLD/"
		iMarNetFolder	= ObsFolder+"/LestersReportData/"
		GlodapDir	= ObsFolder+"/GLODAP/"
		GLODAPv2Dir	= ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
		OSUDir		= ObsFolder+"OSU/"
		CCIDir		= ObsFolder+"CCI/"
		icFold		= ObsFolder+"/InitialConditions/"

		# New eORCA1 grid		
		orcaGridfn 	= '/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'
	else:
		print "This only runs on jasmin."
		#assert 0		


  	




	
	

	#####
	# Coordinate dictionairy
	# These are python dictionairies, one for each data source and model.
	# This is because each data provider seems to use a different set of standard names for dimensions and time.
	# The 'tdict' field is short for "time-dictionary". 
	#	This is a dictionary who's indices are the values on the netcdf time dimension.
	#	The tdict indices point to a month number in python numbering (ie January = 0)
	# 	An example would be, if a netcdf uses the middle day of the month as it's time value:
	#		tdict = {15:0, 45:1 ...}	


	medusaCoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	icCoords 	= {'t':'time_counter', 'z':'nav_lev', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.	
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	takahashiCoords	= {'t':'index_t', 'z':'index_z','lat': 'LAT', 'lon': 'LON', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	glodapCoords	= {'t':'index_t', 'z':'depth',  'lat': 'latitude', 'lon': 'longitude', 'cal': 'standard','tdict':[] }
	glodapv2Coords	= {'t':'time',    'z':'Pressure','lat':'lat',      'lon':'lon',        'cal': '',        'tdict':{0:0,} }
	mldCoords	= {'t':'index_t', 'z':'index_z','lat':'lat','lon': 'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	cciCoords	= {'t':'index_t', 'z':'index_z','lat': 'lat',      'lon': 'lon', 'cal': 'standard','tdict':['ZeroToZero'] }



			



	dataD = {}		
	modeldataD = {}
		
	for jobID in jobs:
	
		#####
		# Location of images directory
		# the imagedir is where the analysis images will be saved.
		imagedir	 = ukp.folder('images/'+jobID+'/timeseries')
		shelvedir 	= ukp.folder("/group_workspaces/jasmin2/ukesm/BGC_data/"+getuser()+"/shelves/timeseries/"+jobID)
		
		
		def listModelDataFiles(jobID, filekey, datafolder, annual):
			if annual:
				return sorted(glob(datafolder+jobID+"/"+jobID+"o_1y_*_"+filekey+".nc"))
			else:
				return sorted(glob(datafolder+jobID+"/"+jobID+"o_1m_*_"+filekey+".nc"))

		av = ukp.AutoVivification()							
		if 'DrakePassageTransport' in analysisKeys:
			name = 'DrakePassageTransport'
			####
			# Note that this will only work with the eORCA1grid.
		
			# coordinates of Drake Passage
			LON=219
			LAT0=79
			LAT1=109
		
			nc = Dataset(orcaGridfn,'r')
			e2u = nc.variables['e2u'][LAT0:LAT1,LON]
			umask = nc.variables['umask'][:,LAT0:LAT1,LON]
			nc.close()			

			def drake(nc,keys):
				e3u = nc.variables['e3u'][0,:,LAT0:LAT1,LON]
				velo = nc.variables['vozocrtx'][0,:,LAT0:LAT1,LON]
				return np.sum(velo*e3u*e2u*umask)*1.e-6			
								
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_U', MEDUSAFolder_pref, annual)
			av[name]['dataFile'] 	= ''
			
			av[name]['modelcoords'] = medusaCoords 	
			av[name]['datacoords'] 	= medusaCoords

			av[name]['modeldetails']= {'name': name, 'vars':['e3u','vozocrtx',], 'convert': drake,'units':'Sv'}		    	
				
			av[name]['regions'] 		=  ['regionless',]		
			av[name]['datadetails']  	= {'name':'','units':'',}
			av[name]['layers'] 		=  ['layerless',]		
			av[name]['metrics']		= ['metricless',]
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 1



		if 'NorthernTotalIceArea' in analysisKeys or 'SouthernTotalIceArea' in analysisKeys or 'TotalIceArea' in analysisKeys:
		    for name in ['NorthernTotalIceArea','SouthernTotalIceArea','TotalIceArea']:
		    	if name not in analysisKeys:continue
		
			nc = Dataset(orcaGridfn,'r')
			area = nc.variables['e2t'][:] * nc.variables['e1t'][:]			
			tmask = nc.variables['tmask'][0,:,:]
			lat = nc.variables['nav_lat'][:,:]
			nc.close()			

			def calcTotalIceArea(nc,keys):	#Global
				arr = nc.variables[keys[0]][:].squeeze() * area
				return np.ma.masked_where(tmask==0,arr).sum()/1E12
			
			def calcTotalIceAreaN(nc,keys): # North
				arr = nc.variables[keys[0]][:].squeeze() * area
				return np.ma.masked_where((tmask==0)+(lat<0.),arr).sum()/1E12

			def calcTotalIceAreaS(nc,keys): # South
				arr = nc.variables[keys[0]][:].squeeze() * area
				return np.ma.masked_where((tmask==0)+(lat>0.),arr).sum()/1E12
								
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', MEDUSAFolder_pref, annual)												
			av[name]['dataFile'] 		= ''
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= medusaCoords

		    	if name in ['NorthernTotalIceArea',]:
				av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov',], 'convert': calcTotalIceAreaN,'units':'1E6 km^2'}		    	
			#	av[name]['regions'] 		=  ['NorthHemisphere',]	
				
		    	if name in ['SouthernTotalIceArea',]:
				av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov',], 'convert': calcTotalIceAreaS,'units':'1E6 km^2'}		    	
			#	av[name]['regions'] 		=  ['SouthHemisphere',]				

		    	if name in ['TotalIceArea',]:
				av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov',], 'convert': calcTotalIceArea,'units':'1E6 km^2'}		    	
			#	av[name]['regions'] 		=  ['Global',]					
			av[name]['regions'] 		=  ['regionless',]		
			
			av[name]['datadetails']  	= {'name':'','units':'',}
			#av[name]['layers'] 		=  ['Surface',]
			av[name]['layers'] 		=  ['layerless',]		
			av[name]['metrics']		= ['metricless',]
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'CICE'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 1
		
		
				
		if 'GlobalExportRatio' in analysisKeys:
		
			def calcExportRatio(nc,keys):
				a = (nc.variables['SDT__100'][:] +nc.variables['FDT__100'][:]).sum()/ (nc.variables['PRD'][:] +nc.variables['PRN'][:] ).sum()
				#a = np.ma.masked_where(a>1.01, a)
				return 	a
			
			name = 'ExportRatio'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', MEDUSAFolder_pref, annual)								
				
			av[name]['dataFile'] 		= ""
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= maredatCoords
			av[name]['modeldetails'] 	= {'name': name, 'vars':['SDT__100','FDT__100' ,'PRD','PRN',], 'convert': calcExportRatio,'units':''}
			av[name]['datadetails']  	= {'name':'','units':'',}
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]
			av[name]['metrics']		= ['metricless',]
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'MEDUSA'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn		
			av[name]['Dimensions']		= 1
		
		
		#####
		# Total 
		if 'IntPP_OSU' in analysisKeys:
			noOSU = True
			nc = Dataset(orcaGridfn,'r')
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
				return arr.sum()
			
			name = 'TotalIntegratedPrimaryProduction'
			if annual:
				av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', MEDUSAFolder_pref, annual)							
				if noOSU:	av[name]['dataFile']            = ''
				else:		av[name]['dataFile'] 		= OSUDir +"/standard_VGPM.SeaWIFS.global.average.nc"
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= glodapCoords

		        av[name]['modeldetails']        = {'name': 'IntPP', 'vars':['PRN' ,'PRD'], 'convert': medusadepthInt,'units':'Gt/yr'}
			if noOSU: 
			        av[name]['datadetails']         = {'name': '', 'units':''}

			else:
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
					return arr.sum()
			       	av[name]['datadetails']         = {'name': 'IntPP', 'vars':['NPP',], 'convert': osuconvert,'units':'Gt/yr'}

			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]
			av[name]['metrics']		= ['metricless',]
			if noOSU:	av[name]['datasource']          = ''
			else:		av[name]['datasource'] 		= 'OSU'
			av[name]['model']		= 'MEDUSA'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn		
			av[name]['Dimensions']		= 1
		

		if 'TotalAirSeaFlux' in analysisKeys:
			name = 'TotalAirSeaFluxCO2'	
			nc = Dataset(orcaGridfn,'r')
			area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
			nc.close()
		
			def eOrcaTotal(nc,keys):
				factor =  365.25 * 12./1000. / 1.E15
				arr = nc.variables['CO2FLUX'][:].squeeze() * factor	# mmolC/m2/d
				if arr.ndim ==3:
					for i in np.arange(arr.shape[0]):
						arr[i] = arr[i]*area
				elif arr.ndim ==2: arr = arr*area
				else: assert 0
				return arr.sum()
					
			def takaTotal(nc,keys):
				arr = nc.variables['TFLUXSW06'][:].squeeze()	# 10^12 g Carbon year^-1
				arr = -1.E12* arr /1.E15#/ 365.				#g Carbon/day
				#area = nc.variables['AREA_MKM2'][:].squeeze() *1E12	# 10^6 km^2
				#fluxperarea = arr/area
				return arr.sum()
				# area 10^6 km^2
				# flux:  10^15 g Carbon month^-1. (GT)/m2/month

			


			av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', MEDUSAFolder_pref, annual)								
			if annual:
				av[name]['dataFile'] 		=  TakahashiFolder+'takahashi_2009_Anual_sumflux_2006c_noHead.nc'		
			else:	
				av[name]['dataFile'] 		=  TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'				
				print "Air Sea Flux CO2 monthly not implemented"
				assert 0
				#av[name]['dataFile'] 		=  TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
				
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= takahashiCoords
			av[name]['modeldetails'] 	= {'name': 'AirSeaFluxCO2', 'vars':['CO2FLUX',], 'convert': eOrcaTotal,'units':'Pg C/yr'}
			av[name]['datadetails']  	= {'name': 'AirSeaFluxCO2', 'vars':['TFLUXSW06','AREA_MKM2'], 'convert': takaTotal,'units':'Pg C/yr'}
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]
			av[name]['metrics']		= ['metricless',]
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'MEDUSA'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 2					
										
			noTaka = True
			if noTaka:
				av[name]['datadetails'] =  {'name': '',	'units':''}
				av[name]['dataFile']	= ''
				av[name]['datasource']  = ''
							

		if  'Iron' in analysisKeys:
			name = 'Iron'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', MEDUSAFolder_pref, annual)								
				
			av[name]['dataFile'] 		= icFold+"/UKESM_fields_1860_eORCA1_small.nc"
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= icCoords
			av[name]['modeldetails']	= {'name': name, 'vars':['FER',], 'convert': ukp.mul1000, 'units':'umolFe/m3'}
			av[name]['datadetails']  	= {'name': name, 'vars':['FER',], 'convert': ukp.mul1000, 'units':'umolFe/m3'}
			av[name]['layers'] 		= layerList
			av[name]['regions'] 		= regionList
			av[name]['metrics']		= metricList
			av[name]['datasource'] 		= 'InititialCondition'
			av[name]['model']		= 'MEDUSA'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3

		if 'N' in analysisKeys:
			name = 'Nitrate'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', MEDUSAFolder_pref, annual)				
			if annual:
				av[name]['dataFile'] 		=  WOAFolder+'/woa13_all_n00_01.nc'
			else:
				av[name]['dataFile'] 		=  WOAFolder+'/nitrate_monthly_1deg.nc'
				
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['DIN',], 'convert': ukp.NoChange,'units':'mmol N/m^3'}
			av[name]['datadetails']  	= {'name': name, 'vars':['n_an',], 'convert': ukp.NoChange,'units':'mmol N/m^3'}
	
			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList
		
			#av[name]['layers'] 		= ['Surface','300m',]#'1000m',]#'Surface - 300m',]'100m',
			#av[name]['regions'] 		= regionList#['Global',]#'NorthAtlanticOcean','SouthAtlanticOcean',]#'NorthAtlantic']
			av[name]['metrics']		= metricList #['mean','median', ]

			av[name]['datasource'] 		= 'WOA'
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3
		
		if 'Si' in analysisKeys:
			name = 'Silicate'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', MEDUSAFolder_pref, annual)				
			if annual:	
				av[name]['dataFile'] 		= WOAFolder+'woa13_all_i00_01.nc'
			else:
				av[name]['dataFile'] 		= WOAFolder+'wsilicate_monthly_1deg.nc'
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['SIL',],  'convert': ukp.NoChange,'units':'mmol Si/m^3'}
			av[name]['datadetails']  	= {'name': name, 'vars':['i_an',], 'convert': ukp.NoChange,'units':'mmol Si/m^3'}
		
			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'WOA'
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3		
	
		if 'O2' in analysisKeys:
			name = 'Oxygen'
			if annual:
				av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', MEDUSAFolder_pref, annual)		
				av[name]['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
				
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['OXY',], 'convert': ukp.NoChange,'units':'mmol O2/m^3'}	
			av[name]['datadetails']  	= {'name': name, 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol O2/m^3'}

			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'WOA'
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3


		if 'DIC' in analysisKeys:
	
			def convertkgToM3(nc,keys):
				return nc.variables[keys[0]][:]* 1.027
				
			name = 'DIC'
		
			av[name]['modelFiles'] 		= listModelDataFiles(jobID, 'ptrc_T', MEDUSAFolder_pref, annual)				
			av[name]['dataFile'] 		= GLODAPv2Dir+ 'GLODAPv2.tco2.historic.nc'
				
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= glodapv2Coords
	
			av[name]['modeldetails'] 	= {'name': 'DIC', 'vars':['DIC',],  'convert': ukp.NoChange,'units':'mmol C/m^3'}
			av[name]['datadetails']  	= {'name': 'DIC', 'vars':['tco2',], 'convert': ukp.convertkgToM3,'units':'mmol C/m^3'}
	
			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'GLODAP'
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn		
			av[name]['Dimensions']		= 3
		
		if 'Alk' in analysisKeys:
			def convertmeqm3TOumolkg(nc,keys):
				return nc.variables[keys[0]][:]* 1.027
		
			name = 'Alkalinity'
			if annual:		
				av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', MEDUSAFolder_pref, annual)		
				av[name]['dataFile'] 	=  GlodapDir+'Alk.nc'
			else:
				print "Alkalinity data not available for monthly Analysis"
				assert 0
				
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= glodapCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['ALK',], 'convert': ukp.NoChange,'units':'meq/m^3',}
			av[name]['datadetails']  	= {'name': name, 'vars':['Alk',], 'convert': convertmeqm3TOumolkg,'units':'meq/m^3',}
	
		#	av[name]['layers'] 		=  ['Surface','100m','300m','1000m',]
		#	av[name]['regions'] 		= regionList
			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList		
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'GLODAP'
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3		
		
						
		if 'T' in analysisKeys:
			name = 'Temperature'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', MEDUSAFolder_pref, annual)										
			if annual:		
				#av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'woa13_decav_t00_01v2.nc'
			else:
				#av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'temperature_monthly_1deg.nc'
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['votemper',], 'convert': ukp.NoChange,'units':'degrees C'}
			av[name]['datadetails']  	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}
	
			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList	
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'WOA'
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3
					
		if 'S' in analysisKeys:
			name = 'Salinity'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', MEDUSAFolder_pref, annual)												
			if annual:
				#av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'woa13_decav_s00_01v2.nc'
			else:
				#av[name]['modelFiles']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'salinity_monthly_1deg.nc'
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['vosaline',], 'convert': ukp.NoChange,'units':'PSU'}	
			av[name]['datadetails']  	= {'name': name, 'vars':['s_an',], 'convert': ukp.NoChange,'units':'PSU'}

			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= regionList		
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'WOA'
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3
							
					
					

		for name in av.keys():
			print "------------------------------------------------------------------"	
			print "analysis-Timeseries.py:\tBeginning to call timeseriesAnalysis for ", name

			if len(av[name]['modelFiles']) == 0:
				print "analysis-Timeseries.py:\tWARNING:\tmodel files are not found:",av[name]['modelFiles'], jobID
				if strictFileCheck: assert 0		

			modelfilesexists = [os.path.exists(f) for f in av[name]['modelFiles']]
			if False in modelfilesexists:
				print "analysis-Timeseries.py:\tWARNING:\tnot model files do not all exist:",av[name]['modelFiles'] 
				if strictFileCheck: assert 0
			
			
			if av[name]['dataFile']!='':
			   if not os.path.exists(av[name]['dataFile']):
				print "analysis-Timeseries.py:\tWARNING:\tdata file is not found:",av[name]['dataFile']
				if strictFileCheck: assert 0
			
			#####
			# time series and traffic lights.
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
				grid		= av[name]['modelgrid'],
				gridFile	= av[name]['gridFile'],
				clean 		= False,
				noNewFiles	= True,
			)
			dataD[(jobID,name )] = tsa.dataD
			modeldataD[(jobID,name )] = tsa.modeldataD
	
	#####
	# Data now loaded, making plots next:
	
	#print modeldataD.keys()
	for k in modeldataD.keys():
		print "Model Data D:",k
	
	for name in av.keys():
		timesD  = {}
		arrD	= {}
		
		for jobID in jobs:
			if name in ['Iron','Nitrate','Silicate','Oxygen','Temperature','Salinity', 'Alkalinity','DIC',]:
				mdata = modeldataD[(jobID,name )][('Global', 'Surface', 'mean')]
				title = ' '.join(['Global', 'Surface', 'Mean',  getLongName(name)])
			else:
				mdata = modeldataD[(jobID,name )][('regionless', 'layerless', 'metricless')]
				title = getLongName(name)
			timesD[jobID] 	= sorted(mdata.keys())
			arrD[jobID]	= [mdata[t] for t in timesD[jobID]]
		
		#####
		# To account for changing units.
		if name in ['TotalAirSeaFluxCO2', 'TotalAirSeaFlux']: 
		    for j in arrD.keys():
			if j in ['u-ad980','u-af123','u-af725',
				 'u-ae742','u-af139','u-af578', 'u-af728']:
				arrD[j] = np.ma.array(arrD[j]) * 5.09369e-7 
		
		for ts in ['Together','Separate']:
		    for ls in ['Both','movingaverage',]:#'','Both',]:			
			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename=ukp.folder('images/TimeseriesCompare/')+name+'_'+ts+'_'+ls+'.png',
				units = '',
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
			)
	

timeseries_compare()

