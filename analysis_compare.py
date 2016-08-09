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


def timeseries_compare():
	### strategy here is a simple wrapper.
	#Total Global Primary production, Drake passage transport and Air Sea CO2 flux.
	jobs = ['u-af123', 'u-af139','u-af420','u-af421']
	annual = True
	
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
	#			analysisKeys.append('TotalIceArea')		# work in progress	
	#			analysisKeys.append('NorthernTotalIceArea')	# work in progress	
	#			analysisKeys.append('SouthernTotalIceArea')	# work in progress	


	analysisKeys.append('DrakePassageTransport')	# DrakePassageTransport		
	#analysisKeys.append('TotalAirSeaFlux')          # work in progress              
	#analysisKeys.append('IntPP_OSU')                # OSU Integrated primpary production    


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
		

		for name in av.keys():
			print "------------------------------------------------------------------"	
			print "analysis-Timeseries.py:\tBeginning to call timeseriesAnalysis for ", name

			if len(av[name]['modelFiles']) == 0:
				print "analysis-Timeseries.py:\tWARNING:\tmodel files are not found:",av[name]['modelFiles']
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
	print modeldataD.keys()
	for k in modeldataD.keys():
		print "Model Data D:",k,modeldataD[k]
		
		
	

timeseries_compare()

