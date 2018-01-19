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
"""
.. module:: analysis_compare
   :platform: Unix
   :synopsis: A script to produce an intercomparison of multiple runs the time series analyses.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

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
import os,sys,fnmatch
from getpass import getuser
from collections import defaultdict

#####	
# Load specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis
from timeseries import profileAnalysis
from timeseries import timeseriesPlots as tsp 
try:	from bgcvaltools.pftnames import getLongName
except:	from pftnames import getLongName
from bgcvaltools.mergeMonthlyFiles import mergeMonthlyFiles,meanDJF
from alwaysInclude import  alwaysInclude
from makeReport import comparehtml5Maker
import paths

from comparison.shifttimes import shifttimes
from config.configToDict import configToDict
from bgcvaltools.dataset import dataset

def listModelDataFiles(jobID, filekey, datafolder, annual,year=''):
	if year == '':
		if annual:
			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1y_*_"+filekey+".nc"))
		else:
			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1m_*_"+filekey+".nc"))
	else:
		if annual:
			print datafolder+jobID+"/"+jobID+"o_1y_*"+year+"????_"+filekey+".nc"
			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1y_*"+year+"????_"+filekey+".nc"))
		else:
			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1m_*"+year+"????_"+filekey+".nc"))
				
def timeseries_compare(colours,physics=True,bio=False,debug=False,year0=False,analysisname='',jobDescriptions={},lineThicknesses= defaultdict(lambda:1)):
	### strategy here is a simple wrapper.
	# It's a little cheat-y, as I'm copying straight from analysis_timeseries.py
	
	jobs = sorted(colours.keys())

	if analysisname=='':
		imageFolder = paths.imagedir+'/TimeseriesCompare/'
        	if len(jobs)==1:   imageFolder+= jobs[0]
	       	elif len(jobs)==2: imageFolder+= jobs[0]+'_and_'+jobs[1]
        	else: 		   imageFolder+= str(len(jobs))+'Runs_'+jobs[0]
	else:
		imageFolder = paths.imagedir+'/TimeseriesCompare/'+analysisname


	annual = True
	strictFileCheck = False

	analysisKeys = []

	if physics:
		analysisKeys.append('DrakePassageTransport')	# DrakePassageTransport		
	        analysisKeys.append('TotalIceArea')             # TotalIceArea  
	        analysisKeys.append('NorthernTotalIceArea')     # North TotalIceArea
	        analysisKeys.append('SouthernTotalIceArea')     # South TotalIceArea
	        analysisKeys.append('TotalIceExtent')           # work in progress      
        	analysisKeys.append('NorthernTotalIceExtent')   # work in progress      
	        analysisKeys.append('SouthernTotalIceExtent')   # work in progress      
	        analysisKeys.append('WeddelIceExent')   # work in progress      	        

        	analysisKeys.append('AMOC_26N')
	        analysisKeys.append('AMOC_32S')
                analysisKeys.append('ADRC_26N')                # AMOC 26N                        	        
        	analysisKeys.append('Temperature')             # WOA Temperature
	        analysisKeys.append('Salinity')                # WOA Salinity   
                analysisKeys.append('MLD')                      # MLD
        	analysisKeys.append('MaxMonthlyMLD')               # MLD Monthly max           
	        analysisKeys.append('MinMonthlyMLD')               # MLD Monthly min        

               	analysisKeys.append('ZonalCurrent')             # Zonal Veloctity
               	analysisKeys.append('MeridionalCurrent')        # Meridional Veloctity
               	analysisKeys.append('VerticalCurrent')          # Vertical Veloctity  
               	analysisKeys.append('GlobalMeanTemperature')    # Global Mean Temperature  
               	analysisKeys.append('VolumeMeanTemperature')    # Global Mean Temperature                 	             	 	       
               	analysisKeys.append('GlobalMeanSalinity')       # Global Mean Salinity
               	analysisKeys.append('IcelessMeanSST')    	# Global Mean Surface Temperature with no ice               		
               	
		analysisKeys.append('sowaflup')			# Net Upward Water Flux 
		analysisKeys.append('sohefldo')			# Net downward Water Flux 			
		analysisKeys.append('sofmflup')			# Water flux due to freezing/melting
		analysisKeys.append('sosfldow')			# Downward salt flux
		analysisKeys.append('soicecov')			# Ice fraction			               
                analysisKeys.append('sossheig')                 # SSH

		
	if bio:
		analysisKeys.append('TotalAirSeaFlux')          # work in progress             
                analysisKeys.append('NoCaspianAirSeaFluxCO2')   # work in progress                      
 
                analysisKeys.append('AirSeaFlux')               # work in progress              
		analysisKeys.append('IntPP_OSU')                # OSU Integrated primpary production    
		analysisKeys.append('GlobalExportRatio')

		analysisKeys.append('N')                        # WOA Nitrate
		analysisKeys.append('Si')                       # WOA Siliate
		analysisKeys.append('O2')                       # WOA Oxygen
		analysisKeys.append('Iron')
		analysisKeys.append('Alk')	
		analysisKeys.append('DIC')

		analysisKeys.append('CHD')
		analysisKeys.append('CHN')
		analysisKeys.append('CHL')		
		analysisKeys.append('DiaFrac')
                analysisKeys.append('DMS')
		analysisKeys.append('Dust')                     # Dust
        	analysisKeys.append('TotalDust')                # Total Dust

#                analysisKeys.append('DMS_ARAN')

                analysisKeys.append('DTC')
							
      	  	analysisKeys.append('TotalOMZVolume')           # Total Oxygen Minimum zone Volume
       	 	analysisKeys.append('OMZThickness')             # Oxygen Minimum Zone Thickness
        	analysisKeys.append('OMZMeanDepth')             # Oxygen Minimum Zone mean depth      
        
        if debug:
        	####
        	# Supercedes other flags.
		analysisKeys = []
#		analysisKeys.append('ERSST')
                analysisKeys.append('DrakePassageTransport')    # DrakePassageTransport         
                analysisKeys.append('AMOC_26N')
#                analysisKeys.append('NoCaspianAirSeaFluxCO2')   # work in progress                      
#               	analysisKeys.append('VolumeMeanTemperature')    # Global Mean Temperature                 	             	 	       
#		analysisKeys.append('CHD')
#		analysisKeys.append('CHN')
#		analysisKeys.append('DiaFrac')
#                analysisKeys.append('AMOC_26N')
#                analysisKeys.append('MLD')
#                analysisKeys.append('DMS')
#                analysisKeys.append('N')                        # WOA Nitrate
#                analysisKeys.append('Si')                       # WOA Siliate

#                analysisKeys.append('TotalAirSeaFlux')          # work in progress              
#                analysisKeys.append('AirSeaFlux')          # work in progress              

                #analysisKeys.append('ADRC_26N')                # AMOC 26N                        
#                analysisKeys.append('VerticalCurrent')          # Vertical Veloctity           
#                analysisKeys.append('sossheig')                 # SSH
                analysisKeys.append('NorthernTotalIceArea')     # North TotalIceArea
                analysisKeys.append('SouthernTotalIceArea')     # South TotalIceArea
                analysisKeys.append('TotalIceArea')     	#  TotalIceArea
                analysisKeys.append('NorthernTotalIceExtent')   # work in progress      
                analysisKeys.append('SouthernTotalIceExtent')   # work in progress      
                analysisKeys.append('WeddelIceExent')   # work in progress       
                analysisKeys.append('TotalIceExtent')           # work in progress      

	
                analysisKeys.append('GlobalMeanTemperature')
#                analysisKeys.append('GlobalMeanSalinity')

#               	analysisKeys.append('quickSST')    		# Area Weighted Mean Surface Temperature
#       	  	analysisKeys.append('TotalOMZVolume')           # Total Oxygen Minimum zone Volume
#       	 	analysisKeys.append('OMZThickness')             # Oxygen Minimum Zone Thickness
#        	analysisKeys.append('OMZMeanDepth')             # Oxygen Minimum Zone mean depth    
#		analysisKeys.append('O2')                       # WOA Oxygen        	
#       	if bio ==False:return
#       	if physics == True:return  
                               	
	layerList 	= ['Surface',]
	metricList 	= ['mean',]
  	regionList	= ['Global',]

	vmtregionList = ['Global', 'ignoreInlandSeas','Equator10','SouthernOcean','ArcticOcean',  'Remainder','NorthernSubpolarAtlantic','NorthernSubpolarPacific','WeddelSea',]
	level3 = ['DMS',]	
	
	#####
	# paths:
	orcaGridfn 	= paths.orcaGridfn #'/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'
	if annual:	WOAFolder = paths.WOAFolder_annual
	else:		WOAFolder = paths.WOAFolder	
  	




	
	

	#####
	# Coordinate dictionairy
	# These are python dictionairies, one for each data source and model.
	# This is because each data provider seems to use a different set of standard names for dimensions and time.
	# The 'tdict' field is short for "time-dictionary". 
	#	This is a dictionary who's indices are the values on the netcdf time dimension.
	#	The tdict indices point to a month number in python numbering (ie January = 0)
	# 	An example would be, if a netcdf uses the middle day of the month as it's time value:
	#		tdict = {15:0, 45:1 ...}	


	medusaCoords 	= {'t':'index_t', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	medusaUCoords 	= {'t':'index_t', 'z':'depthu', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	medusaVCoords 	= {'t':'index_t', 'z':'depthv', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	medusaWCoords 	= {'t':'index_t', 'z':'depthw', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	
	icCoords 	= {'t':'time_counter', 'z':'nav_lev', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.	
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	takahashiCoords	= {'t':'index_t', 'z':'index_z','lat': 'LAT', 'lon': 'LON', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}	
	glodapCoords	= {'t':'index_t', 'z':'depth',  'lat': 'latitude', 'lon': 'longitude', 'cal': 'standard','tdict':[] }
	glodapv2Coords	= {'t':'time',    'z':'Pressure','lat':'lat',      'lon':'lon',        'cal': '',        'tdict':{0:0,} }
	mldCoords	= {'t':'index_t', 'z':'index_z','lat':'lat','lon': 'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	cciCoords	= {'t':'index_t', 'z':'index_z','lat': 'lat',      'lon': 'lon', 'cal': 'standard','tdict':['ZeroToZero'] }
	dmsCoords	= {'t':'time',    'z':'depth',  'lat':'Latitude',  'lon': 'Longitude','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	godasCoords 	= {'t':'index_t',    'z':'level',  'lat': 'lat',      'lon': 'lon', 'cal': 'standard','tdict':['ZeroToZero'] }




	dataD = {}		
	modeldataD = {}
		
	for jobID in jobs:
	
		#####
		# Location of images directory
		# the imagedir is where the analysis images will be saved.
		imagedir	 = ukp.folder(paths.imagedir+'/'+jobID+'/timeseries')
		shelvedir 	= ukp.folder(paths.shelvedir+"/timeseries/"+jobID)
		
		


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
								
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_U', paths.ModelFolder_pref, annual)
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


		icekeys = ['NorthernTotalIceArea','SouthernTotalIceArea','TotalIceArea','NorthernTotalIceExtent','SouthernTotalIceExtent','WeddelIceExent','TotalIceExtent']
		if len(set(icekeys).intersection(set(analysisKeys))):
		    for name in icekeys:
		    	if name not in analysisKeys:continue
		
			nc = Dataset(orcaGridfn,'r')
			area = nc.variables['e2t'][:] * nc.variables['e1t'][:]			
			tmask = nc.variables['tmask'][0,:,:]
			lat = nc.variables['nav_lat'][:,:]
                        lon = nc.variables['nav_lon'][:,:]
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
		
			def calcTotalIceExtent(nc,keys):	#Global
				return np.ma.masked_where((tmask==0)+(nc.variables[keys[0]][:].squeeze()<0.15),area).sum()/1E12
			
			def calcTotalIceExtentN(nc,keys): # North
				return np.ma.masked_where((tmask==0)+(nc.variables[keys[0]][:].squeeze()<0.15)+(lat<0.),area).sum()/1E12

			def calcTotalIceExtentS(nc,keys): # South
				return np.ma.masked_where((tmask==0)+(nc.variables[keys[0]][:].squeeze()<0.15)+(lat>0.),area).sum()/1E12

			weddelmask = (lat<-80.)+(lat>-65.)+(lon <-60.)+(lon > -20.)
			def calcTotalIceExtentWS(nc,keys): # South
				return np.ma.masked_where((tmask==0)+(nc.variables[keys[0]][:].squeeze()<0.15)+weddelmask,area).sum()/1E12
											
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)												
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
		    	if name in ['NorthernTotalIceExtent',]:
				av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov',], 'convert': calcTotalIceExtentN,'units':'1E6 km^2'}		    	
			#	av[name]['regions'] 		=  ['NorthHemisphere',]	
				
		    	if name in ['SouthernTotalIceExtent',]:
				av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov',], 'convert': calcTotalIceExtentS,'units':'1E6 km^2'}		    	
			#	av[name]['regions'] 		=  ['SouthHemisphere',]				

		    	if name in ['TotalIceExtent',]:
				av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov',], 'convert': calcTotalIceExtent,'units':'1E6 km^2'}		    	
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
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)								
				
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
				av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)							
				if noOSU:	av[name]['dataFile']            = ''
				else:		av[name]['dataFile'] 		= paths.OSUDir +"/standard_VGPM.SeaWIFS.global.average.nc"
			
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

			


			av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)								
			if annual:
				av[name]['dataFile'] 		=  paths.TakahashiFolder+'takahashi_2009_Anual_sumflux_2006c_noHead.nc'		
			else:	
				av[name]['dataFile'] 		=  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'				
				print "Air Sea Flux CO2 monthly not implemented"
				assert 0
				#av[name]['dataFile'] 		=  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
				
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

		if 'NoCaspianAirSeaFluxCO2' in analysisKeys:
			name = 'NoCaspianAirSeaFluxCO2'
			nc = Dataset(paths.orcaGridfn,'r')
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
				return arr



			av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)

			av[name]['modelcoords'] 	= medusaCoords

			av[name]['modeldetails']        = {'name': 'AirSeaFluxCO2', 'vars':['CO2FLUX',], 'convert': eOrcaTotal,'units':'Pg C/yr'}

			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['ignoreCaspian',]
			av[name]['metrics']		= ['sum',]
			av[name]['model']		= 'MEDUSA'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 2

			av[name]['datacoords'] 		= {'name': '', 'units':''}
			av[name]['datadetails'] 	=  {'name': '',	'units':''}
			av[name]['dataFile']		= ''
			av[name]['datasource']  	= ''

		if 'AirSeaFlux' in analysisKeys:

		        #nc = dataset(paths.orcaGridfn,'r')
		        #area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
		        #nc.close()

		        def eOrcaTotal(nc,keys):
		                factor =  12./1000.
		                arr = nc.variables['CO2FLUX'][:].squeeze()      # mmolC/m2/d
		                #if arr.ndim ==3:
		                #       for i in np.arange(arr.shape[0]):
		                #               arr[i] = arr[i]*area
		                #elif arr.ndim ==2: arr = arr*area
		                #else: assert 0
		                return arr * factor

		        def takaTotal(nc,keys):
		                arr = nc.variables['TFLUXSW06'][:].squeeze()    # 10^12 g Carbon year^-1
		                arr = -1.E12* arr / 365.                                #g Carbon/day
		                factor = -1.E12/(365. ) # convert to #/ 1.E12
		                area = nc.variables['AREA_MKM2'][:].squeeze() *1E12     # 10^6 km^2
		                fluxperarea = arr/area
		                #arr = arr*area #* 1.E24        # converts area into m^2
		                #print arr.sum(), arr.sum()*factor
		                return fluxperarea
		                # area 10^6 km^2
		                # flux:  10^15 g Carbon month^-1. (GT)/m2/month


		        name = 'AirSeaFluxCO2'

		        av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)
		        if annual:
		                av[name]['dataFile']            =  paths.TakahashiFolder+'takahashi_2009_Anual_sumflux_2006c_noHead.nc'
		        else:
		                av[name]['dataFile']            =  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
		                print "Air Sea Flux CO2 monthly not implemented"
		                assert 0
		                #av[name]['dataFile']           =  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'

		        av[name]['modelcoords']         = medusaCoords
		        av[name]['datacoords']          = takahashiCoords
		        av[name]['modeldetails']        = {'name': 'AirSeaFluxCO2', 'vars':['CO2FLUX',], 'convert': eOrcaTotal,'units':'g C/m2/day'}
		        av[name]['datadetails']         = {'name': 'AirSeaFluxCO2', 'vars':['TFLUXSW06','AREA_MKM2'], 'convert': takaTotal,'units':'g C/m2/day'}
		        av[name]['layers']              = ['Surface',]
		        av[name]['regions']             = regionList
		        av[name]['metrics']             = metricList
		        av[name]['datasource']          = ''
		        av[name]['model']               = 'MEDUSA'
		        av[name]['modelgrid']           = 'eORCA1'
		        av[name]['gridFile']            = paths.orcaGridfn
		        av[name]['Dimensions']          = 2

							
		if 'CHD' in analysisKeys or  'CHN' in analysisKeys:
		    for name in ['CHD','CHN',]:
			if name not in analysisKeys: continue
		
			av[name]['modelFiles']  	= listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)		
			av[name]['dataFile'] 		= ''
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= ''
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':[name,], 'convert': ukp.NoChange,'units':'mg C/m^3'}
			av[name]['datadetails']  	= {'name': '', 'units':''}
	
			av[name]['layers'] 		= ['Surface','100m',] 	# CCI is surface only, it's a satellite product.
			av[name]['regions'] 		= regionList 
			av[name]['metrics']		= metricList	#['mean','median', ]
		
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 2		

		if 'CHL' in analysisKeys:
			name = 'Chlorophyll'
			av[name]['modelFiles']  	= listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
			av[name]['dataFile'] 		= ''

			av[name]['modelcoords'] 	= medusaCoords
			av[name]['datacoords'] 		= maredatCoords

			av[name]['modeldetails'] 	= {'name': name, 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}
			av[name]['datadetails']  	= {'name': '', 'units':''}
		
			av[name]['layers'] 		= ['Surface','100m','200m',]
			av[name]['regions'] 		= regionList
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= ''
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 3
				
		if 'DiaFrac' in analysisKeys:
					
			name = 'DiaFrac'
			def caldiafrac(nc,keys):
				return 100.*nc.variables[keys[0]][:].squeeze()/(nc.variables[keys[0]][:].squeeze()+nc.variables[keys[1]][:].squeeze())
		
			av[name]['modelFiles']  	= listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)		
			av[name]['dataFile'] 		= ''
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= ''
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['CHD','CHN',], 'convert': caldiafrac,'units':'%'}
			av[name]['datadetails']  	= {'name': '', 'units':''}
	
			av[name]['layers'] 		= ['Surface','100m',] 	# CCI is surface only, it's a satellite product.
			av[name]['regions'] 		= regionList 
			av[name]['metrics']		= metricList	#['mean','median', ]
		
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 2	
        	if 'DTC' in analysisKeys:
	            for name in ['DTC',]:
	                if name not in analysisKeys: continue

                	av[name]['modelFiles']          = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
        	        av[name]['dataFile']            = ''

	                av[name]['modelcoords']         = medusaCoords
                	av[name]['datacoords']          = ''

        	        av[name]['modeldetails']        = {'name': name, 'vars':[name,], 'convert': ukp.mul1000,'units':'umol-C/m3'}
	                av[name]['datadetails']         = {'name': '', 'units':''}

                	av[name]['layers']              = ['3000m',]#'100m',]         # CCI is surface only, it's a satellite product.
        	        av[name]['regions']             = regionList
	                av[name]['metrics']             = metricList    #['mean','median', ]

        	        av[name]['datasource']          = ''
	                av[name]['model']               = 'MEDUSA'

                	av[name]['modelgrid']           = 'eORCA1'
        	        av[name]['gridFile']            = paths.orcaGridfn
	                av[name]['Dimensions']          = 3
		
		
		
		if  'Iron' in analysisKeys:
			name = 'Iron'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)								
				
			av[name]['dataFile'] 		= paths.icFold+"/UKESM_fields_1860_eORCA1_small.nc"
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
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)				
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
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)				
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
				av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)		
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

		if 'OMZMeanDepth' in analysisKeys:
			if annual:
				av['OMZMeanDepth']['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
				av['OMZMeanDepth']['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
			else:
				print "OMZ Thickness not implemented for monthly data"
				assert 0
			
			nc = Dataset(paths.orcaGridfn,'r')
			depths   	= np.abs(nc.variables['gdepw' ][:])
			tmask 		= nc.variables['tmask'][:]						
			nc.close()			

			omzthreshold = 20.
				
			def modelMeanOMZdepth(nc,keys):
				o2 = nc.variables[keys[0]][:].squeeze()
				meandepth = np.ma.masked_where((o2>omzthreshold)+o2.mask + (tmask==0),depths).mean(0)
				if meandepth.max() in [0.,0]: return np.array([0.,])						
				return np.ma.masked_where(meandepth==0., meandepth)

			def woaMeanOMZdepth(nc,keys):
				o2 = nc.variables[keys[0]][:].squeeze() *44.661
				pdepths = np.zeros_like(o2) 
				lons = nc.variables['lon'][:]
				lats = nc.variables['lat'][:]			
				wdepths = np.abs(nc.variables['depth'][:])
			
				for y,lat in enumerate(lats):
				    for x,lon in enumerate(lons):			  
					pdepths[:,y,x] = wdepths
				wmeanDepth = np.ma.masked_where((o2>omzthreshold)+o2.mask,pdepths).mean(0).data
				print "woaMeanOMZdepth",wmeanDepth.min(),wmeanDepth.mean(),wmeanDepth.max()
				#assert 0
			
				if wmeanDepth.max() in [0.,0]: return np.array([1000.,])
				return np.ma.masked_where(wmeanDepth==0., wmeanDepth)
				
			av['OMZMeanDepth']['modelcoords'] 	= medusaCoords 	
			av['OMZMeanDepth']['datacoords'] 	= woaCoords
	
			av['OMZMeanDepth']['modeldetails'] 	= {'name': 'OMZMeanDepth', 'vars':['OXY',],  'convert': modelMeanOMZdepth,'units':'m'}
			av['OMZMeanDepth']['datadetails']  	= {'name': 'OMZMeanDepth', 'vars':['o_an',], 'convert': woaMeanOMZdepth,'units':'m'}		
	
			av['OMZMeanDepth']['layers'] 		= ['layerless',] 
			av['OMZMeanDepth']['regions'] 		= regionList
			av['OMZMeanDepth']['metrics']		= metricList

			av['OMZMeanDepth']['datasource'] 	= 'WOA'
			av['OMZMeanDepth']['model']		= 'MEDUSA'

			av['OMZMeanDepth']['modelgrid']		= 'eORCA1'
			av['OMZMeanDepth']['gridFile']		= paths.orcaGridfn	
			av['OMZMeanDepth']['Dimensions']	= 2		
		
		
		
		if 'OMZThickness' in analysisKeys or 'OMZThickness50' in analysisKeys:
			if 'OMZThickness' in analysisKeys and 'OMZThickness50' in analysisKeys:
				print "Only run one of these at a time"
				assert 0
			
		
			if annual:
				av['OMZThickness']['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
				av['OMZThickness']['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
			else:
				print "OMZ Thickness not implemented for monthly data"
				assert 0
			
			nc = Dataset(paths.orcaGridfn,'r')
			thickness   	= nc.variables['e3t' ][:]
			tmask 		= nc.variables['tmask'][:]								
			nc.close()			

			if 'OMZThickness' in analysisKeys: 	omzthreshold = 20.
			if 'OMZThickness50' in analysisKeys: 	omzthreshold = 50.
				
			def modelOMZthickness(nc,keys):
				o2 = nc.variables[keys[0]][:].squeeze()
				totalthick = np.ma.masked_where((o2>omzthreshold)+o2.mask+ (tmask==0),thickness).sum(0).data
				if totalthick.max() in [0.,0]: return np.array([0.,])
			
				return np.ma.masked_where(totalthick==0., totalthick)
				#return np.ma.masked_where((arr>omzthreshold) + (arr <0.) + arr.mask,thickness).sum(0)

			
			def woaOMZthickness(nc,keys):
				o2 = nc.variables[keys[0]][:].squeeze() *44.661
				pthick = np.zeros_like(o2) 
				lons = nc.variables['lon'][:]
				lats = nc.variables['lat'][:]			
				zthick  = np.abs(nc.variables['depth_bnds'][:,0] - nc.variables['depth_bnds'][:,1])

				for y,lat in enumerate(lats):
				    for x,lon in enumerate(lons):			  
					pthick[:,y,x] = zthick
				totalthick = np.ma.masked_where((o2>omzthreshold)+o2.mask,pthick).sum(0).data
			
				if totalthick.max() in [0.,0]: return np.array([0.,])
				return np.ma.masked_where(totalthick==0., totalthick)

			av['OMZThickness']['modelcoords'] 	= medusaCoords 	
			av['OMZThickness']['datacoords'] 		= woaCoords
	
			av['OMZThickness']['modeldetails'] 	= {'name': 'OMZThickness', 'vars':['OXY',], 'convert': modelOMZthickness,'units':'m'}
			av['OMZThickness']['datadetails']  	= {'name': 'OMZThickness', 'vars':['o_an',], 'convert': woaOMZthickness,'units':'m'}
	
			av['OMZThickness']['layers'] 		= ['layerless',] 
			av['OMZThickness']['regions'] 		= regionList
			av['OMZThickness']['metrics']		= metricList

			av['OMZThickness']['datasource'] 		= 'WOA'
			av['OMZThickness']['model']		= 'MEDUSA'

			av['OMZThickness']['modelgrid']		= 'eORCA1'
			av['OMZThickness']['gridFile']		= paths.orcaGridfn	
			av['OMZThickness']['Dimensions']		= 2		
		
		
		

		if 'TotalOMZVolume' in analysisKeys or 'TotalOMZVolume50' in analysisKeys:
			if 'TotalOMZVolume' in analysisKeys and 'TotalOMZVolume50' in analysisKeys:
				print "Only run one of these at a time"
				assert 0
				
			if annual:
				av['TotalOMZVolume']['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
				av['TotalOMZVolume']['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
			else:
				print "OMZ volume not implemented for monthly data"
				assert 0
			
			nc = Dataset(paths.orcaGridfn,'r')
			try:	
				pvol   = nc.variables['pvol' ][:]
				tmask = nc.variables['tmask'][:]
			except:
				tmask = nc.variables['tmask'][:]			
				area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
				pvol = nc.variables['e3t'][:] *area			
				pvol = np.ma.masked_where(tmask==0,pvol)
			nc.close()			

			if 'TotalOMZVolume' in analysisKeys:	omzthreshold = 20.
			if 'TotalOMZVolume50' in analysisKeys:	omzthreshold = 50.	
			
			def modelTotalOMZvol(nc,keys):
				arr = np.ma.array(nc.variables[keys[0]][:].squeeze())
				return np.ma.masked_where((arr>omzthreshold) + pvol.mask + arr.mask,pvol).sum()
	
			
			def woaTotalOMZvol(nc,keys):
				arr = nc.variables[keys[0]][:].squeeze() *44.661
				#area = np.zeros_like(arr[0])
				pvol = np.zeros_like(arr) 
				#np.ma.masked_wjhere(arr.mask + (arr <0.)+(arr >1E10),np.zeros_like(arr))
				lons = nc.variables['lon'][:]
				lats = nc.variables['lat'][:]			
				#lonbnds = nc.variables['lon_bnds'][:]
				latbnds = nc.variables['lat_bnds'][:]
				zthick  = np.abs(nc.variables['depth_bnds'][:,0] - nc.variables['depth_bnds'][:,1])
			
				for y,lat in enumerate(lats):
					area = ukp.Area([latbnds[y,0],0.],[latbnds[y,1],1.])
					for z,thick in enumerate(zthick):
						pvol[z,y,:] = np.ones_like(lons)*area*thick
					
				return np.ma.masked_where(arr.mask + (arr >omzthreshold)+(arr <0.),pvol).sum()
				
			av['TotalOMZVolume']['modelcoords'] 	= medusaCoords 	
			av['TotalOMZVolume']['datacoords'] 		= woaCoords
	
			av['TotalOMZVolume']['modeldetails'] 	= {'name': 'TotalOMZVolume', 'vars':['OXY',], 'convert': modelTotalOMZvol,'units':'m^3'}
			av['TotalOMZVolume']['datadetails']  	= {'name': 'TotalOMZVolume', 'vars':['o_an',], 'convert': woaTotalOMZvol,'units':'m^3'}
	
			av['TotalOMZVolume']['layers'] 		= ['layerless',] 
			av['TotalOMZVolume']['regions'] 	= ['regionless',]
			av['TotalOMZVolume']['metrics']		= ['metricless', ]

			av['TotalOMZVolume']['datasource'] 		= 'WOA'
			av['TotalOMZVolume']['model']		= 'MEDUSA'

			av['TotalOMZVolume']['modelgrid']		= 'eORCA1'
			av['TotalOMZVolume']['gridFile']		= paths.orcaGridfn	
			av['TotalOMZVolume']['Dimensions']		= 1	
			


		if 'DIC' in analysisKeys:
	
			def convertkgToM3(nc,keys):
				return nc.variables[keys[0]][:]* 1.027
				
			name = 'DIC'
		
			av[name]['modelFiles'] 		= listModelDataFiles(jobID, 'ptrcG_T', paths.ModelFolder_pref, annual)				
			av[name]['dataFile'] 		= paths.GLODAPv2Dir+ 'GLODAPv2.tco2.historic.nc'
				
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
				av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)		
				av[name]['dataFile'] 	=  paths.GlodapDir+'Alk.nc'
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
		
						
		if 'Temperature' in analysisKeys:
			name = 'Temperature'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)										
			if annual:		
				#av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'woa13_decav_t00_01v2.nc'
			else:
				#av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
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

                if 'ERSST' in analysisKeys:
                        name = 'ERSST'
                        av[name]['modelFiles']  = "/group_workspaces/jasmin/esmeval/example_data/bgc/ERSST.v4/sst.mnmean.v4.nc"
		        ERSSTCoords     = {'t':'time',    'z':'',  'lat': 'lat',      'lon': 'lon', 'cal': 'standard','tdict':['ZeroToZero'] }

                        av[name]['modelcoords']         = ERSSTCoords
                        #av[name]['datacoords']          = woaCoords

                        av[name]['modeldetails']        = {'name': name, 'vars':['sst',], 'convert': ukp.NoChange,'units':'degrees C'}
                        av[name]['datadetails']         = {'name': name, 'vars':[], 'units':''}

                        av[name]['layers']              = ['layerless',]
                        av[name]['regions']             = regionList
                        av[name]['metrics']             = metricList

                        av[name]['datasource']          = ''
                        av[name]['model']               = 'ERSST'

                        av[name]['modelgrid']           = 'ERSST_2g'
                        av[name]['gridFile']            = '/group_workspaces/jasmin/esmeval/example_data/bgc/ERSST.v4/sst.mnmean.v4.nc'
                        av[name]['Dimensions']          = 2
		

		if 'GlobalMeanTemperature' in analysisKeys:
			name = 'GlobalMeanTemperature'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)										
			av[name]['dataFile'] 	= ''		
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords

			nc = Dataset(paths.orcaGridfn,'r')
			try:	
				pvol   = nc.variables['pvol' ][:]
				tmask = nc.variables['tmask'][:]
			except:
				tmask = nc.variables['tmask'][:]			
				area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
				pvol = nc.variables['e3t'][:] *area			
				pvol = np.ma.masked_where(tmask==0,pvol)
			nc.close()
		
			def sumMeanLandMask(nc,keys):
				temperature = np.ma.masked_where(tmask==0,nc.variables[keys[0]][:].squeeze())
				return (temperature*pvol).sum()/(pvol.sum())
			
			av[name]['modeldetails'] 	= {'name': name, 'vars':['votemper',], 'convert': sumMeanLandMask,'units':'degrees C'}
			av[name]['datadetails']  	= {'name': '', 'units':''}		
	
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]
			av[name]['metrics']		= ['metricless',]

			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 1

		if 'VolumeMeanTemperature' in analysisKeys:
			name = 'VolumeMeanTemperature'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
			av[name]['dataFile'] 	= ''
			av[name]['modelcoords'] = medusaCoords
			av[name]['datacoords'] 	= woaCoords

		        def sumMeanLandMask(nc,keys):
		        	assert 0

			av[name]['modeldetails'] 	= {'name': name, 'vars':['votemper',], 'convert': sumMeanLandMask,'units':'degrees C'}
			av[name]['datadetails']  	= {'name': '', 'units':''}
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= vmtregionList
			av[name]['metrics']		= ['wcvweighted',]
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 2


		if 'GlobalMeanSalinity' in analysisKeys:
			name = 'GlobalMeanSalinity'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)										
			av[name]['dataFile'] 	= ''		
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords

			nc = Dataset(paths.orcaGridfn,'r')
			try:	
				pvol   = nc.variables['pvol' ][:]
				tmask = nc.variables['tmask'][:]
			except:
				tmask = nc.variables['tmask'][:]			
				area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
				pvol = nc.variables['e3t'][:] *area			
				pvol = np.ma.masked_where(tmask==0,pvol)
			nc.close()
		
			def sumMeanLandMask(nc,keys):
				temperature = np.ma.masked_where(tmask==0,nc.variables[keys[0]][:].squeeze())
				return (temperature*pvol).sum()/(pvol.sum())
			
			av[name]['modeldetails'] 	= {'name': name, 'vars':['vosaline',], 'convert': sumMeanLandMask,'units':'degrees C'}
			av[name]['datadetails']  	= {'name': '', 'units':''}		
	
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]
			av[name]['metrics']		= ['metricless',]

			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 1
			
		if 'quickSST' in analysisKeys:
			name = 'quickSST'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)										
		

			nc = Dataset(paths.orcaGridfn,'r')
			ssttmask = nc.variables['tmask'][0]			
			area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
			area = np.ma.masked_where(ssttmask==0,area)
			nc.close()
		
			def meanLandMask(nc,keys):
				#### works like no change, but applies a mask.
				#print "meanLandMask:",ssttmask.shape,nc.variables[keys[0]][0,0].shape
				temperature = np.ma.masked_where(ssttmask==0,nc.variables[keys[0]][0,0].squeeze())
				print "meanLandMask:",nc.variables['time_counter'][:],temperature.mean(),(temperature*area).sum()/(area.sum())
				return (temperature*area).sum()/(area.sum())
			
					
			if annual:		
				#av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
				av[name]['dataFile'] 		= ''#WOAFolder+'woa13_decav_t00_01v2.nc'
			else:
				#av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
				av[name]['dataFile'] 		= ''#WOAFolder+'temperature_monthly_1deg.nc'
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['votemper',], 'convert': meanLandMask,'units':'degrees C'}
			av[name]['datadetails']  	= {'name': '', 'units':''}
	
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]	
			av[name]['metrics']		= ['metricless',]

			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 1

		if 'IcelessMeanSST' in analysisKeys:
			name = 'IcelessMeanSST'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)										
			av[name]['dataFile'] 	= ''		
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords

			nc = Dataset(paths.orcaGridfn,'r')
			tmask = nc.variables['tmask'][:]			
			area_full = nc.variables['e2t'][:] * nc.variables['e1t'][:]
			nc.close()
		
			def calcIcelessMeanSST(nc,keys):
				#### works like no change, but applies a mask.
				icecov = nc.variables['soicecov'][:].squeeze()
				sst = nc.variables['votemper'][:,0,].squeeze()
				sst = np.ma.masked_where((tmask[0]==0)+(icecov>0.15)+sst.mask,sst)
				area=  np.ma.masked_where(sst.mask,area_full)
				val = (sst*area).sum()/(area.sum())
				print "calcIcelessMeanSST", sst.shape,area.shape, val
				return val
		
			
			av[name]['modeldetails'] 	= {'name': name, 'vars':['soicecov','votemper',], 'convert': calcIcelessMeanSST,'units':'degrees C'}
			av[name]['datadetails']  	= {'name': '', 'units':''}		
			#av[name]['datadetails']  	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}
	
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= ['regionless',]
			av[name]['metrics']		= ['metricless',]

			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 1
				
					
		if 'Salinity' in analysisKeys:
			name = 'Salinity'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)												
			if annual:
				#av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'woa13_decav_s00_01v2.nc'
			else:
				#av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
				av[name]['dataFile'] 		= WOAFolder+'salinity_monthly_1deg.nc'
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= woaCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['vosaline',], 'convert': ukp.NoChange,'units':'PSU'}	
			av[name]['datadetails']  	= {'name': name, 'vars':['s_an',], 'convert': ukp.NoChange,'units':'PSU'}

			av[name]['layers'] 		=  layerList
			av[name]['regions'] 		= vmtregionList
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'WOA'
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= orcaGridfn
			av[name]['Dimensions']		= 3

		if 'MLD' in analysisKeys:

		        def mldapplymask(nc,keys):
		                mld = np.ma.array(nc.variables[keys[0]][:])
		                return np.ma.masked_where((nc.variables[keys[1]][:]==0.)+mld.mask+(mld==1.E9),mld)


		        name = 'MLD'
		        av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
		        av[name]['dataFile']            = paths.MLDFolder+"mld_DT02_c1m_reg2.0-annual.nc"       #mld_DT02_c1m_reg2.0.nc"

		        av[name]['modelcoords']         = medusaCoords
		        av[name]['datacoords']          = mldCoords

		        av[name]['modeldetails']        = {'name': 'mld', 'vars':['somxl010',],   'convert': ukp.NoChange,'units':'m'}
		        av[name]['datadetails']         = {'name': 'mld', 'vars':['mld','mask',], 'convert': mldapplymask,'units':'m'}

		        av[name]['layers']              = ['layerless',]#'Surface - 1000m','Surface - 300m',]#'depthint']
		        mldregions =regionList
			mldregions.extend(['WeddelSea',])# 'LabradorSea', 'NorwegianSea'])				        
#		        mldregions.extend(['NordicSea', 'LabradorSea', 'NorwegianSea'])
		        av[name]['regions']             = mldregions
		        av[name]['metrics']             = metricList

		        av[name]['datasource']          = 'IFREMER'
		        av[name]['model']               = 'NEMO'

		        av[name]['modelgrid']           = 'eORCA1'
		        av[name]['gridFile']            = paths.orcaGridfn
		        av[name]['Dimensions']          = 2

	        if 'MaxMonthlyMLD' in analysisKeys or 'MinMonthlyMLD' in analysisKeys:

        	        #/group_workspaces/jasmin2/ukesm/BGC_data/u-ad371/monthlyMLD/MetOffice_data_licence.325210916
                	monthlyFiles = glob(paths.ModelFolder_pref+'/'+jobID+'/monthlyMLD/'+jobID+'o_1m_*_grid_T.nc')
	                if len(monthlyFiles):
        	                maxmldfiles = mergeMonthlyFiles(monthlyFiles,outfolder='',cal=medusaCoords['cal'])

                	        for name in ['MaxMonthlyMLD','MinMonthlyMLD']:
                        	        if name not in analysisKeys:continue

	                                def mldapplymask(nc,keys):
        	                                mld = np.ma.array(nc.variables[keys[0]][:]).max(0)
                	                        mld = np.ma.masked_where((nc.variables[keys[1]][:]==0.)+mld.mask+(mld==1.E9),mld)
                        	                return mld

                                	def mldmonthlymask(nc,keys):
                                        	mld = np.ma.array(np.ma.abs(nc.variables[keys[0]][:])).max(0)
	                                        mld = np.ma.masked_where(mld.mask+(mld.data>1.E10),mld)
        	                                return mld

                	                def mldapplymask_min(nc,keys):
                        	                mld = np.ma.array(nc.variables[keys[0]][:]).min(0)
                                	        mld = np.ma.masked_where((nc.variables[keys[1]][:]==0.)+mld.mask+(mld==1.E9),mld)
                                        	return mld

	                                def mldmonthlymask_min(nc,keys):
        	                                mld = np.ma.array(np.ma.abs(nc.variables[keys[0]][:])).min(0)
                	                        mld = np.ma.masked_where(mld.mask+(mld.data>1.E10),mld)
                        	                return mld


	                                av[name]['modelFiles']          = maxmldfiles #listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
        	                        av[name]['dataFile']            = paths.MLDFolder+"mld_DT02_c1m_reg2.0.nc"      #mld_DT02_c1m_reg2.0.nc"

                	                av[name]['modelcoords']         = medusaCoords
                        	        av[name]['datacoords']          = mldCoords

        	                        if name =='MaxMonthlyMLD':
	                                        av[name]['modeldetails']        = {'name': 'mld', 'vars':['somxl010',],   'convert': mldmonthlymask,'units':'m'}
                	                        av[name]['datadetails']         = {'name': 'mld', 'vars':['mld','mask',], 'convert': mldapplymask,'units':'m'}
                        	        if name =='MinMonthlyMLD':
                                	        av[name]['modeldetails']        = {'name': 'mld', 'vars':['somxl010',],   'convert': mldmonthlymask_min,'units':'m'}
                                        	av[name]['datadetails']         = {'name': 'mld', 'vars':['mld','mask',], 'convert': mldapplymask_min,'units':'m'}

	                                av[name]['layers']              = ['layerless',]
                                	av[name]['regions']             = regionList
	                                av[name]['metrics']             = metricList
        	                        av[name]['datasource']          = 'IFREMER'
                	                av[name]['model']               = 'NEMO'
                        	        av[name]['modelgrid']           = 'eORCA1'
                                	av[name]['gridFile']            = paths.orcaGridfn
	                                av[name]['Dimensions']          = 2



		if 'ZonalCurrent' in analysisKeys:
			name = 'ZonalCurrent'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_U', paths.ModelFolder_pref, annual)												
			if annual:
				av[name]['dataFile'] 		= paths.GODASFolder+'ucur.clim.nc'
			
			av[name]['modelcoords'] 	= medusaUCoords 	
			av[name]['datacoords'] 		= godasCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['vozocrtx',], 'convert': ukp.mul1000,'units':'mm/s'}	
			av[name]['datadetails']  	= {'name': name, 'vars':['ucur',], 'convert': ukp.NoChange,'units':'mm/s'}

			av[name]['layers'] 		= layerList
			av[name]['regions'] 		= regionList		
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'GODAS'
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 3
		

		if 'MeridionalCurrent' in analysisKeys:
			name = 'MeridionalCurrent'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_V', paths.ModelFolder_pref, annual)												
			if annual:
				av[name]['dataFile'] 		= paths.GODASFolder+'vcur.clim.nc'
			
			av[name]['modelcoords'] 	= medusaVCoords 	
			av[name]['datacoords'] 		= godasCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['vomecrty',], 'convert': ukp.mul1000,'units':'mm/s'}	
			av[name]['datadetails']  	= {'name': name, 'vars':['vcur',], 'convert': ukp.NoChange,'units':'mm/s'}

			av[name]['layers'] 		= layerList
			av[name]['regions'] 		= regionList		
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'GODAS'
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 3
			
		if 'VerticalCurrent' in analysisKeys:
			name = 'VerticalCurrent'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_W', paths.ModelFolder_pref, annual)												
			if annual:
				av[name]['dataFile'] 		= paths.GODASFolder+'dzdt.clim.nc'
			
			av[name]['modelcoords'] 	= medusaWCoords 	
			av[name]['datacoords'] 		= godasCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['vovecrtz',], 'convert': ukp.mul1000000,'units':'um/s'}
			av[name]['datadetails']  	= {'name': name, 'vars':['dzdt',], 'convert': ukp.NoChange,'units':'um/s'}

			av[name]['layers'] 		= layerList
			av[name]['regions'] 		= regionList		
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'GODAS'
			av[name]['model']		= 'NEMO'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 3
			


							


	
		if 'AMOC_26N' in analysisKeys or 'AMOC_32S' in analysisKeys or 'ADRC_26N' in analysisKeys:
			# Note that this will only work with the eORCAgrid.
			latslice26N = slice(227,228)
			latslice32S = slice(137,138)
			e3v,e1v,tmask,alttmask = {},{},{},{}
		    	for name in ['AMOC_26N','AMOC_32S','ADRC_26N']:
		    		if name not in analysisKeys:continue

				####
				if name in ['AMOC_26N','ADRC_26N']: 	latslice = latslice26N
				if name == 'AMOC_32S': 	latslice = latslice32S

				# Load grid data
				nc = Dataset(paths.orcaGridfn,'r')
				e3v[name] = nc.variables['e3v'][:,latslice,:]	# z level height 3D
				e1v[name] = nc.variables['e1v'][latslice,:]	#
				tmask[name] = nc.variables['tmask'][:,latslice,:]
				nc.close()

				# load basin mask
				nc = Dataset('data/basinlandmask_eORCA1.nc','r')
				alttmask[name] = nc.variables['tmaskatl'][latslice,:]	# 2D Atlantic mask
				nc.close()

			def calc_amoc32S(nc,keys):
				name = 'AMOC_32S'
				zv = np.ma.array(nc.variables['vomecrty'][...,latslice32S,:]) # m/s
				atlmoc = np.array(np.zeros_like(zv[0,:,:,0]))
				e2vshape = e3v[name].shape
				for la in range(e2vshape[1]):		#ji, y
	 			  for lo in range(e2vshape[2]):	#jj , x,
	 			    if int(alttmask[name][la,lo]) == 0: continue
				    for z in range(e2vshape[0]): 	# jk
	 			    	if int(tmask[name][z,la,lo]) == 0: 	   continue
	 			    	if np.ma.is_masked(zv[0,z,la,lo]): continue
	 			    	atlmoc[z,la] = atlmoc[z,la] - e1v[name][la,lo]*e3v[name][z,la,lo]*zv[0,z,la,lo]/1.E06

	 			####
	 			# Cumulative sum from the bottom up.
	 			for z in range(73,1,-1):
	 				atlmoc[z,:] = atlmoc[z+1,:] + atlmoc[z,:]
				return np.ma.max(atlmoc)

			def amoc26N_array(nc,keys,amocname='AMOC_26N'):
				zv = np.ma.array(nc.variables['vomecrty'][...,latslice26N,:]) # m/s
				atlmoc = np.array(np.zeros_like(zv[0,:,:,0]))
				e2vshape = e3v[amocname].shape
				for la in range(e2vshape[1]):		#ji, y
	 			  for lo in range(e2vshape[2]):		#jj , x,
	 			    if int(alttmask[amocname][la,lo]) == 0: continue
				    for z in range(e2vshape[0]): 	# jk
	 			    	if int(tmask[amocname][z,la,lo]) == 0: 	   continue
	 			    	if np.ma.is_masked(zv[0,z,la,lo]): continue
	 			    	atlmoc[z,la] = atlmoc[z,la] - e1v[amocname][la,lo]*e3v[amocname][z,la,lo]*zv[0,z,la,lo]/1.E06

	 			####
	 			# Cumulative sum from the bottom up.
	 			for z in range(73,1,-1):
	 				atlmoc[z,:] = atlmoc[z+1,:] + atlmoc[z,:]
				#return np.ma.max(atlmoc)
				return atlmoc

			def calc_amoc26N(nc,keys):
				return np.ma.max(amoc26N_array(nc,keys,amocname='AMOC_26N'))
			
			def calc_min_amoc26N(nc,keys):
				return np.ma.min(amoc26N_array(nc,keys,amocname='ADRC_26N'))

			
		    	for name in ['AMOC_26N','AMOC_32S','ADRC_26N']:
		    		if name not in analysisKeys:continue

				av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_V', paths.ModelFolder_pref, annual)
				av[name]['dataFile'] 	= ''

				av[name]['modelcoords'] = medusaCoords
				av[name]['datacoords'] 	= medusaCoords

				if name == 'AMOC_26N': 	av[name]['modeldetails']= {'name': name, 'vars':['vomecrty',], 'convert': calc_amoc26N,'units':'Sv'}
				if name == 'ADRC_26N': 	av[name]['modeldetails']= {'name': name, 'vars':['vomecrty',], 'convert': calc_min_amoc26N,'units':'Sv'}			
				if name == 'AMOC_32S': 	av[name]['modeldetails']= {'name': name, 'vars':['vomecrty',], 'convert': calc_amoc32S,'units':'Sv'}

				av[name]['datadetails']  	= {'name':'','units':'',}
				av[name]['layers'] 		=  ['layerless',]
				av[name]['regions'] 		= ['regionless',]
				av[name]['metrics']		= ['metricless',]
				av[name]['datasource'] 		= ''
				av[name]['model']		= 'NEMO'
				av[name]['modelgrid']		= 'eORCA1'
				av[name]['gridFile']		= paths.orcaGridfn
				av[name]['Dimensions']		= 1	
					
		if 'DMS' in analysisKeys:
			name = 'DMS'
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)#[::30]
			if annual:		
				av[name]['dataFile'] 		= paths.DMSDir+'DMSclim_mean.nc'
			else:	av[name]['dataFile'] 		= ''
			
			av[name]['modelcoords'] 	= medusaCoords 	
			av[name]['datacoords'] 		= dmsCoords
	
			av[name]['modeldetails'] 	= {'name': name, 'vars':['DMS_ARAN',], 'convert': ukp.mul1000000,'units':'umol/m3'}
			av[name]['datadetails']  	= {'name': name, 'vars':['DMS',], 'convert': ukp.NoChange,'units':'umol/m3'}
	
			av[name]['layers'] 		= ['layerless',]
			av[name]['regions'] 		= regionList	
			av[name]['metrics']		= metricList

			av[name]['datasource'] 		= 'Lana'
			av[name]['model']		= 'MEDUSA'

			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 2				

		if 'TotalDust' in analysisKeys:
		        name = 'TotalDust'
		        av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)[:]
		        av[name]['dataFile']    = paths.Dustdir+'mahowald.orca100_annual.nc'

		        av[name]['modelcoords']         = medusaCoords
		        av[name]['datacoords']          = medusaCoords

		        nc = Dataset(paths.orcaGridfn,'r')
		        masked_area = nc.variables['e2t'][:] * nc.variables['e1t'][:]*nc.variables['tmask'][0]
		        nc.close()

		        def datadustsum(nc,keys):
		                #factors are:
		                # 0.035: iron as a fraction of total dust
		                # 1e6: convert from kmol -> mmol
		                # 1e-12: convert from mol to Gmol
		                # 0.00532: solubility factor for iron
		                # 55.845: atmoic mass of iron (g>mol conversion)
		                # (24.*60.*60.*365.25): per second to per year

		                dust = nc.variables[keys[0]][:]
		                dust[:,:,234:296,295:348] = 0.
		                dust[:,:,234:248,285:295] = 0.
		                dust[:,:,228:256,290:304] = 0.
		                return (masked_area*dust).sum() *0.035 * 1.e6*1.e-12 *0.00532*(24.*60.*60. *365.25)/  55.845

		        def modeldustsum(nc,keys):
		                dust = nc.variables[keys[0]][:]
		                dust[:,234:296,295:348] = 0.
		                dust[:,234:248,285:295] = 0.
		                dust[:,228:256,290:304] = 0.
		                return (masked_area*dust).sum() *1.E-12 *365.25

		        av[name]['modeldetails']        = {'name': name, 'vars':['AEOLIAN',], 'convert': modeldustsum,'units':'Gmol Fe/yr'}
		        av[name]['datadetails']         = {'name': name, 'vars':['dust_ann',], 'convert': datadustsum,'units':'Gmol Fe/yr'}

		        av[name]['layers']              = ['layerless',]
		        av[name]['regions']             = ['regionless',]
		        av[name]['metrics']             = ['metricless',]

		        av[name]['datasource']          = 'Mahowald'
		        av[name]['model']               = 'MEDUSA'

		        av[name]['modelgrid']           = 'eORCA1'
		        av[name]['gridFile']            = paths.orcaGridfn
		        av[name]['Dimensions']          = 1
		        
		if 'Dust' in analysisKeys:
		        name = 'Dust'
		        av[name]['modelFiles']  = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref, annual)[:]
		        av[name]['dataFile']            = paths.Dustdir+'mahowald.orca100_annual.nc'

		        av[name]['modelcoords']         = medusaCoords
		        av[name]['datacoords']          = medusaCoords
		        av[name]['modeldetails']        = {'name': name, 'vars':['AEOLIAN',], 'convert': ukp.NoChange,'units':'mmol Fe/m2/d'}

		        def mahodatadust(nc,keys):
		                #factors are:
		                # 0.035: iron as a fraction of total dust
		                # 1e6: convert from kmol -> mmol
		                # 0.00532: solubility factor or iron
		                # 55.845: atmoic mass of iron (g>mol conversion)
		                # (24.*60.*60.): per second to per day
		                dust = nc.variables[keys[0]][:]
		                dust[0,0,194:256,295:348] = 0.
		                dust[0,0,194:208,285:295] = 0.
		                dust[0,0,188:216,290:304] = 0.
		                return dust *0.035 * 1.e6 *0.00532*(24.*60.*60.) / 55.845

		        av[name]['datadetails']         = {'name': name, 'vars':['dust_ann',], 'convert': mahodatadust ,'units':'mmol Fe/m2/d'}

		        av[name]['layers']              = ['layerless',]
		        av[name]['regions']             = regionList
		        av[name]['metrics']             = metricList

		        av[name]['datasource']          = 'Mahowald'
		        av[name]['model']               = 'MEDUSA'

		        av[name]['modelgrid']           = 'eORCA1'
		        av[name]['gridFile']            = paths.orcaGridfn
		        av[name]['Dimensions']          = 2


	
		#####
		# North Atlantic Salinity
		#sowaflup = "Net Upward Water Flux" ;
		#sohefldo = "Net Downward Heat Flux" ;
		#sofmflup = "Water flux due to freezing/melting" ;
		#sosfldow = "Downward salt flux" ;
			
		naskeys = ['sowaflup','sohefldo','sofmflup','sosfldow','soicecov','sossheig']
		if len(set(naskeys).intersection(set(analysisKeys))):
		    for name in naskeys:
		    	if name not in analysisKeys:continue

			#nc = Dataset(paths.orcaGridfn,'r')
			#area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
			#tmask = nc.variables['tmask'][0,:,:]
			#lat = nc.variables['nav_lat'][:,:]
			#nc.close()

			nas_files = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
			try:	nc = Dataset(nas_files[0],'r')
			except:
				print "nc does not exist:", name
				continue
			
			if name not in nc.variables.keys():
				print "analysis_timeseries.py:\tWARNING: ",name ,"is not in the model file."
				continue
			av[name]['modelFiles']  	= nas_files
			av[name]['dataFile'] 		= ''

			av[name]['modelcoords'] 	= medusaCoords
			av[name]['datacoords'] 		= medusaCoords		
		       
			nasUnits = {	'sowaflup':"kg/m2/s",
					'sohefldo':"W/m2",
					'sofmflup':"kg/m2/s",
					'sosfldow':"PSU/m2/s",
					'soicecov':'',				
					'sossheig':'m',
				   }
		
			av[name]['modeldetails'] 	= {'name': name[:], 'vars':[name[:],], 'convert': ukp.NoChange,'units':nasUnits[name][:]}

			av[name]['regions'] 		=  ['NordicSea', 'LabradorSea', 'NorwegianSea','Global',]

			av[name]['datadetails']  	= {'name':'','units':'',}
			av[name]['layers'] 		=  ['layerless',]
			av[name]['metrics']		= metricList
			av[name]['datasource'] 		= ''
			av[name]['model']		= 'NEMO'
			av[name]['modelgrid']		= 'eORCA1'
			av[name]['gridFile']		= paths.orcaGridfn
			av[name]['Dimensions']		= 2
		

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
			#dataD[(jobID,name )] = tsa.dataD
			modeldataD[(jobID,name )] = tsa.modeldataD
	
	#####
	# Data now loaded, making plots next:
	for k in modeldataD.keys():
		print "Model Data D:",k

	#####
	# 2D fields for North Atlantic analysis 
	doNorthAtlanticSalt = False
	if doNorthAtlanticSalt:
		nasregionList	= [
				'NordicSea', 'LabradorSea', 'NorwegianSea'
				]
					
		for name in ['sowaflup','sohefldo','sofmflup','sosfldow','soicecov','MLD','sossheig',]:
 		  continue
		  if name not in av.keys():continue
		  for region in nasregionList:
		    for layer in ['layerless',]:
			timesD  = {}
			arrD	= {}
			
			if name in ['',]: metric = 'metricless'
			else:	metric = 'mean'
			
			for jobID in jobs:
                                print name, region,layer, jobID#, mdata, title

				try:
					mdata = modeldataD[(jobID,name )][(region, layer,metric)]
					title = ' '.join([getLongName(n) for n in [region, layer, metric, name]])
				except: continue
				print name, region,layer, jobID, mdata, title
				
				timesD[jobID] 	= sorted(mdata.keys())
				arrD[jobID]	= [mdata[t] for t in timesD[jobID]]
				if jobID == 'u-aj588':
					arrD[jobID]     = np.ma.masked_where(arrD[jobID]==0.,arrD[jobID])

			
			if len(arrD.keys()) ==0:continue	
			
                        units = av[name]['modeldetails']['units']
			
			for ts in ['Together',]:
			    for ls in ['DataOnly',]:
				tsp.multitimeseries(
					timesD, 		# model times (in floats)
					arrD,			# model time series
					data 	= -999,		# in situ data distribution
					title 	= title,
					filename=ukp.folder(imageFolder+'/NAS')+'_'.join(['NAS',name,region,layer,ts,ls+'.png']),
					units = units,
					plotStyle 	= ts,
					lineStyle	= ls,
					colours		= colours,
	                                thicknesses     = lineThicknesses,
				)
	
	####
	for name in ['Temperature','Salinity',]:
		  if name not in av.keys():continue
		  for region in vmtregionList:
		    for layer in ['Surface','500m','1000m']:
			timesD  = {}
			arrD	= {}
		
			for jobID in jobs:
				try:	mdata = modeldataD[(jobID,name )][(region, layer, 'mean')]
				except: continue
				title = ' '.join([region, layer, 'Mean',  getLongName(name)])
				
				#timesD[jobID] 	= sorted(mdata.keys())
				#arrD[jobID]	= [mdata[t] for t in timesD[jobID]]
                                times,datas = shifttimes(mdata, jobID,year0=year0)
                                timesD[jobID]   = times #mdata.keys())
                                arrD[jobID]     = datas #t] for t in timesD[jobID]]

				
                                if jobID == 'u-aj588':
                                        arrD[jobID]     = np.ma.masked_where(arrD[jobID]==0.,arrD[jobID])

	                        #times,datas = shifttimes(mdata, jobID,year0=year0)

			if len(arrD.keys()) ==0:continue
                        units = av[name]['modeldetails']['units']

			for ts in ['Together',]:
			    for ls in ['DataOnly',]:
				tsp.multitimeseries(
					timesD, 		# model times (in floats)
					arrD,			# model time series
					data 	= -999,		# in situ data distribution
					title 	= title,
					filename=ukp.folder(imageFolder)+'_'.join([name,region,layer,ts,ls+'.png']),
					units = units,
					plotStyle 	= ts,
					lineStyle	= ls,
					colours		= colours,
	                                thicknesses     = lineThicknesses,
				)
					
#		assert 0	
		
	####
	# Standard surface:
	for name in av.keys():
		timesD  = {}
		arrD	= {}
		units = av[name]['modeldetails']['units']
		title= ''	
		for jobID in jobs:
			if name in ['Iron','Nitrate','Silicate',
					'Oxygen','Temperature','Salinity','O2',
					 'Alkalinity','DIC',
					 'CHD','CHN','DiaFrac','CHL','Chlorophyll',
	 			 	 'ZonalCurrent','MeridionalCurrent','VerticalCurrent']:
				mdata = modeldataD[(jobID,name )][('Global', 'Surface', 'mean')]
				title = ' '.join(['Global', 'Surface', 'Mean',  getLongName(name)])
			elif name in [  'OMZThickness', 'OMZMeanDepth', 'DMS','MLD','DMS_ARAN','Dust','MaxMonthlyMLD','MinMonthlyMLD',]:
				try:
					mdata = modeldataD[(jobID,name )][('Global', 'layerless', 'mean')]
					title = ' '.join(['Global', getLongName(name)])	
				except: continue
                        elif name in ['DTC',]:
                                mdata = modeldataD[(jobID,name )][('Global', '3000m', 'mean')]
                                title = ' '.join(['Global', '3000m', 'Mean',  getLongName(name)])
                        elif name in ['AirSeaFlux','AirSeaFluxCO2','Alk','Alkalinity',]:
	                        try:
					mdata = modeldataD[(jobID,name )][('ignoreInlandSeas', 'Surface', 'mean')]
		                        title = ' '.join(['ignoreInlandSeas', 'Surface', 'Mean',  getLongName(name)])
				except: continue
                        elif name in ['NoCaspianAirSeaFluxCO2',]:
                                try:
                                        mdata = modeldataD[(jobID,name )][('ignoreCaspian', 'layerless', 'sum')]
                                        title = ' '.join(['Total',  getLongName(name), '(No Caspian)'])
                                except: continue
                        elif name in ['VolumeMeanTemperature',]:
                                try:
                                        mdata = modeldataD[(jobID,name )][('Global', 'layerless', 'wcvweighted')]
                                        title = ' '.join(['Global',  getLongName(name), ])
                                except: continue
						
			elif name in [ 'sowaflup','sohefldo','sofmflup','sosfldow','sossheig', 'soicecov',]:
				
                                try:
					mdata = modeldataD[(jobID,name )][('Global', 'layerless', 'mean')]
				except: continue
                                title = ' '.join(['Global', getLongName(name)])
				
				#####
				# Special hack for these guys.
				#nasregionList	= ['NordicSea', 'LabradorSea', 'NorwegianSea'	]			
				#mdata = modeldataD[(jobID,name )][('regionless', 'layerless', 'mean')]
				#title = getLongName(name)				
			else:
				print name, jobID, name
				mdata = modeldataD[(jobID,name )][('regionless', 'layerless', 'metricless')]
				title = getLongName(name)

			times,datas = shifttimes(mdata, jobID,year0=year0)
			timesD[jobID] 	= times
			arrD[jobID]	= datas

		#####
		# To account for changing units.
		if name in ['TotalAirSeaFluxCO2', 'TotalAirSeaFlux']: 
		    for j in arrD.keys():
			if j in ['u-ad980','u-af123','u-af725',
				 'u-ae742','u-af139','u-af578', 'u-af728']:
				arrD[j] = np.ma.array(arrD[j]) * 5.09369e-7 
               
		if name in ['DMS',]:
                    for j in arrD.keys():
                        if j in ['u-ag914',]:
                                for i,(t,d) in enumerate(zip(timesD[jobID], arrD[jobID])):
                                        if float(t) < 1600.:continue
                                        arrD[j][i] =d/1000000.
		
		for ts in ['Together',]:#'Separate']:
		    for ls in ['DataOnly',]:#'','Both',]:
                        if ls=='' and name not in level3: continue

			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename= ukp.folder(imageFolder)+name+'_'+ts+'_'+ls+'.png',
				units 	= units,
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
                                thicknesses     = lineThicknesses,
			)

                        if name in ['NoCaspianAirSeaFluxCO2','AMOC_26N','ADRC_26N','AMOC_32S']:
			    targetdict = {'NoCaspianAirSeaFluxCO2':0.,
					'AMOC_26N':17.,
					'ADRC_26N':-999,
                                        'AMOC_32S':-999,
					}
                            tsp.multitimeseries(
                                timesD,                 # model times (in floats)
                                arrD,                   # model time series
                                data    = targetdict[name] ,         # in situ gata distribution
				dataname= 'Target',
                                title   = title,
                                filename= ukp.folder(imageFolder)+name+'_'+ts+'_'+'movingav30years'+'.png',
                                units   = units,
                                plotStyle       = ts,
                                lineStyle       = 'movingav30years',
                                colours         = colours,
                                thicknesses     = lineThicknesses,
                            )

	

	
			
							
	####
	# Oxygen at Depth:
  	regionList	= [
  			'Global', 'ignoreInlandSeas',
	  		'SouthernOcean','ArcticOcean',
			'Equator10', 'Remainder',
  			'NorthernSubpolarAtlantic','NorthernSubpolarPacific',
  			]	
	for name in ['Oxygen','O2',]:
	  if name not in av.keys():continue
	  for region in regionList:
	    for layer in ['Surface','500m','1000m']:
		timesD  = {}
		arrD	= {}
		
		for jobID in jobs:
			try:	mdata = modeldataD[(jobID,name )][(region, layer, 'mean')]
			except: continue
			title = ' '.join([region, layer, 'Mean',  getLongName(name)])
			times,datas = shifttimes(mdata, jobID,year0=year0)
			timesD[jobID] 	= times
			arrD[jobID]	= datas
			
			
                units = av[name]['modeldetails']['units']
		
		for ts in ['Together',]:#'Separate']:
		    for ls in ['DataOnly',]:#'','Both',]:						
			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename=ukp.folder(imageFolder+'/Oxygen')+'_'.join([name,region,layer,ts,ls+'.png']),
				units = units,
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
				thicknesses	= lineThicknesses,
			)
	for name in ['VolumeMeanTemperature',]:
	  if name not in av.keys():continue
	  for region in vmtregionList:
	    for layer in ['layerless',]:
		timesD  = {}
		arrD	= {}
		
		for jobID in jobs:
			try:	mdata = modeldataD[(jobID,name )][(region, layer, 'wcvweighted')]
			except: continue
			title = ' '.join([region, layer, getLongName(name)])
			times,datas = shifttimes(mdata, jobID,year0=year0)
			timesD[jobID] 	= times
			arrD[jobID]	= datas
			
                units = av[name]['modeldetails']['units']
		for ts in ['Together',]:#'Separate']:
		    for ls in ['DataOnly',]:#'','Both',]:						
			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename=ukp.folder(imageFolder)+'_'.join([name,region,layer,ts,ls+'.png']),
				units = units,
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
				thicknesses	= lineThicknesses,
			)
        	
        	
	for name in ['DiaFrac','CHD','CHN','CHL','N','Si','Iron','Alk','DIC','Chlorophyll','DMS','Nitrate','Salinity','Silicate','MLD','mld','MaxMonthlyMLD','MinMonthlyMLD','Dust',]:
	  if name not in av.keys():continue
	  
	  for region in regionList:
	    for layer in ['Surface','100m','200m','layerless',]:
	    
		timesD  = {}
		arrD	= {}
		
		for jobID in jobs:
			try:mdata = modeldataD[(jobID,name )][(region, layer, 'mean')]
			except: continue
			title = ' '.join([region, layer, 'Mean',  getLongName(name)])
			times,datas = shifttimes(mdata, jobID,year0=year0)
			timesD[jobID] 	= times
			arrD[jobID]	= datas
						
                units = av[name]['modeldetails']['units']
		
		for ts in ['Together',]:#'Separate']:
		    for ls in ['DataOnly',]:#'','Both',]:			
			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename=ukp.folder(imageFolder+'/BGC')+'_'.join([name,region,layer,ts,ls+'.png']),
				units = units,
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
                                thicknesses     = lineThicknesses,

			)
	for name in ['MLD','mld','MaxMonthlyMLD','MinMonthlyMLD',]:
	  if name not in av.keys():continue
	  for region in vmtregionList:
	    for layer in ['layerless',]:
	    
		timesD  = {}
		arrD	= {}
		
		for jobID in jobs:
			try:mdata = modeldataD[(jobID,name )][(region, layer, 'mean')]
			except: continue
			title = ' '.join([region, layer, 'Mean',  getLongName(name)])
			times,datas = shifttimes(mdata, jobID,year0=year0)
			timesD[jobID] 	= times
			arrD[jobID]	= datas
						
                units = av[name]['modeldetails']['units']
		
		for ts in ['Together',]:#'Separate']:
		    for ls in ['DataOnly',]:#'','Both',]:			
			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename=ukp.folder(imageFolder+'/MLD')+'_'.join([name,region,layer,ts,ls+'.png']),
				units = units,
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
                                thicknesses     = lineThicknesses,
			)
			

	for name in ['DMS',]:
	  if name not in av.keys():continue
	  for region in regionList:
	    for layer in ['Surface','100m','200m',]:
	    
		timesD  = {}
		arrD	= {}
		
		for jobID in jobs:
			try:mdata = modeldataD[(jobID,name )][(region, layer, 'mean')]
			except: continue
			title = ' '.join([region, layer, 'Mean',  getLongName(name)])
	
			timesD[jobID] 	= sorted(mdata.keys())
			arrD[jobID]	= [mdata[t] for t in timesD[jobID]]
                units = av[name]['modeldetails']['units']
		for ts in ['Together',]:#'Separate']:
		    for ls in ['DataOnly',]:#'','Both',]:			
			tsp.multitimeseries(
				timesD, 		# model times (in floats)
				arrD,			# model time series
				data 	= -999,		# in situ data distribution
				title 	= title,
				filename=ukp.folder(imageFolder+'/DMS')+'_'.join([name,region,layer,ts,ls+'.png']),
				units = units,
				plotStyle 	= ts,
				lineStyle	= ls,
				colours		= colours,
                                thicknesses     = lineThicknesses,
			)
	try:
		AllImages = glob(imageFolder, recursive=True)
	except:
		AllImages = []
		for root, dirnames, filenames in os.walk(imageFolder):
    			for filename in fnmatch.filter(filenames, '*.png'):
     			   AllImages.append(os.path.join(root, filename))		
	comparehtml5Maker(
		jobIDs = jobs,
		reportdir = ukp.folder('CompareReports/'+analysisname),
		files = AllImages,
		clean = False,
		doZip= False,
		jobDescriptions = jobDescriptions,
		jobColours = colours,
		)
	



def flatten(lats,lons,dataA,dataB):
	m =  np.ma.array(dataA).mask
	m += np.ma.array(dataB).mask	
	m += np.ma.masked_invalid(dataA/dataB).mask
	
	return  np.ma.masked_where(m, lats).compressed(),\
		np.ma.masked_where(m, lons).compressed(),\
		np.ma.masked_where(m, dataA).compressed(),\
		np.ma.masked_where(m, dataB).compressed()
	
		
		
def CompareTwoRuns(jobIDA,jobIDB,physics=True,bio=False,yearA='',yearB='',debug=True):
	#
	#spatial maps of specific fields

	filetype = []
	if physics:	
		filetype.extend(['grid_T','grid_U','grid_V','grid_W',])
	if bio:	
		filetype.extend(['diad_T','ptrc_T',])		
	filesA = {}
	filesB = {}

        imageFolder = 'images/TimeseriesCompare/'
        imageFolder+= jobIDA+'-'+yearA+'_and_'+jobIDB+'-'+yearB
        
        			
	for ft in filetype:
	        filesA[ft] = listModelDataFiles(jobIDA, ft, paths.ModelFolder_pref, True,year=yearA)[0]
	        filesB[ft] = listModelDataFiles(jobIDB, ft, paths.ModelFolder_pref, True,year=yearB)[0]
	        print "files A:",filesA[ft]
	        print "files B:",filesB[ft]	        
	      
		ncA = Dataset(filesA[ft], 'r')
		ncB = Dataset(filesB[ft], 'r')		
		keys = ukp.intersection(ncA.variables.keys(),ncB.variables.keys())

		lats = ncA.variables['nav_lat'][:]		
		lons = ncA.variables['nav_lon'][:]
				
		for key in keys:
			dataA=0.
			dataB=0.
			if key in alwaysInclude: continue
			if key in ['bounds_lon', 'bounds_lat']:continue
			if ncA.variables[key].ndim==4:
				dataA = ncA.variables[key][0,0,:,:]
				dataB = ncB.variables[key][0,0,:,:]			

			elif ncA.variables[key].ndim==3:
				dataA = ncA.variables[key][0,:,:]
				dataB = ncB.variables[key][0,:,:]	
			elif ncA.variables[key].ndim==1:
				dataA = ncA.variables[key][:,:]
				dataB = ncB.variables[key][:,:]	
			else:
				print "can't plot:",key, ncA.variables[key].ndim
				continue	
			try:	title = ncA.variables[key].long_name
			except:	title = getLongName(key) 
	
			dmin = min([dataA.min(),dataB.min()])
			dmax = max([dataA.max(),dataB.max()])

			print key, lats.shape,lons.shape,dataA.shape,dataB.shape	
			la,lo,data,datb = flatten(lats,lons,dataA,dataB)
			print key, la.shape,lo.shape,data.shape,datb.shape
			
			if 0 in [len(la),len(lo),len(data),len(datb)]:continue
			filename = ukp.folder(imageFolder+'/'+ft)+ft+'-'+key+'.png' 
			ukp.robinPlotQuad(
				lo, la,
				data,
				datb,
				filename,
				titles=[jobIDA+' '+yearA,jobIDB+' '+yearB],
				title=title,
				lon0=0.,
				marble=False,
				drawCbar=True,
				cbarlabel='',
				doLog=False,
				vmin=dmin,vmax=dmax,
				maptype = 'Cartopy',
				scatter=False)	

		
def main():		
        standards = configToDict('config/jobIDcolours.ini')
	thicknesses = defaultdict(lambda: 0.75)
	thicknesses['u-ar783'] = 2.2
        thicknesses['u-at793'] = 2.2
        thicknesses['u-at760'] = 2.2
        thicknesses['u-at628'] = 2.2
        thicknesses['u-at629'] = 2.2


        jobDescriptions = configToDict('config/jobIDdescriptions.ini')   

	try:
		args = argv[1:]
		jobIDs = []
		suite = 'all'
		for job in args:
	                if job == 'debug': 	suite = 'debug'
        	        elif job == 'physics' :	suite='physics'
                	elif job in ['bgc','bio']: suite='bgc'
	                else:  jobIDs.append(job)
	
	except:	
		jobsIDs = []
		suite	= ''

	if len(jobIDs): 
		print "analysis_compare.py:", jobIDs, suite
		if suite=='all':
			phys= 1
			bio = 1
			debug = 0
                if suite=='physics':
                        phys= 1
                        bio = 0
                        debug = 0
                if suite=='bio':
                        phys= 0
                        bio = 1
                        debug = 0
                if suite=='debug':
                        phys= 0
                        bio = 0
                        debug = 1
		try:    
			colours = {i:standards[i] for i in jobIDs}
		except:
			colours = {}
			randomcolours = ['red','blue','purple','green','gold','sienna', 'orange', 'black']
			for i,job in enumerate(jobIDs):
				colours[job] =randomcolours[i]
		name = '_'.join(jobIDs) 
	        timeseries_compare(colours, physics=phys,bio=bio,year0=False,debug=debug,analysisname=name,jobDescriptions=jobDescriptions)
		print "Successful command line comparison"
		exit
	else:
#		jobs = ['u-an869','u-ao586']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='2500-3000',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='Re-couplingTestReduced_2')

                #jobs = ['u-ar783','u-ar538','u-am927i','u-am927ii','u-am927iii','u-aq853iii','u-ar766']
                jobs = ['u-at629','u-at793','u-at628','u-at760','u-am927i','u-am927ii','u-am927iii','u-aq853iii','u-ar766','u-as462','u-as858',]
                colours = {i:standards[i] for i in jobs}
        	thicknesses2 = defaultdict(lambda: 0.75)
        	thicknesses2['u-ar766'] = 1.7
                thicknesses2['u-at629'] = 2.2
                thicknesses2['u-at793'] = 2.2
                thicknesses2['u-at628'] = 2.2
                thicknesses2['u-at760'] = 2.2

                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='4800-5100',
                        jobDescriptions=jobDescriptions,
                        analysisname='UKESM_0.9.3_from4945',
                        lineThicknesses= thicknesses2)

                jobs = ['u-au213','u-au214','u-as462','u-ar977']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                       i:standards[i] for i in jobs},
                       physics=1,
                       bio=1,
                       debug=0,
                       year0='Drift',
                       jobDescriptions=jobDescriptions,
                       analysisname='DecemberRestartTest',
                       lineThicknesses= thicknesses
                )

#               jobs = ['u-au213','u-au214',]
#               colours = {i:standards[i] for i in jobs}
#               timeseries_compare({
#                       i:standards[i] for i in jobs},
#                       physics=1,
#                       bio=1,
#                       debug=1,
#                       year0='Drift3',
#                       jobDescriptions=jobDescriptions,
#                       analysisname='UKESM1_candidates_ice',
#                       lineThicknesses= thicknesses)


                jobs = ['u-ar977','u-as462','u-as858','u-at629','u-at793', 'u-at628',  'u-at760',]
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=1,
                        year0='Drift3',
                        jobDescriptions=jobDescriptions,
                        analysisname='UKESM1_candidates_ice',
                        lineThicknesses= thicknesses)
			
                jobs = ['u-am064','u-am927','u-aq853','u-ar766',]
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='piControl',
                        jobDescriptions=jobDescriptions,
                        analysisname='piControlSpinup',
                        lineThicknesses= thicknesses)

		jobs = ['u-ar977','u-as462','u-as858','u-at629','u-at793', 'u-at628',  'u-at760',]
		colours = {i:standards[i] for i in jobs}
		timeseries_compare({
		        i:standards[i] for i in jobs},
		        physics=1,
		        bio=1,
		        debug=0,
		        year0='Drift3',
		        jobDescriptions=jobDescriptions,
        		analysisname='UKESM1_candidates',
                        lineThicknesses= thicknesses)

                jobs = ['u-at643', 'u-at629', 'u-at646','u-ar977','u-as858',]#'u-ar462']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='Drift2',
                        jobDescriptions=jobDescriptions,
                        analysisname='Set3_UM10.7_STRATROP_OO_0.9.3',
			)

#                jobs = ['u-as051', 'u-as412', 'u-as558']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False, #'Drift2',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM0.9.2-CN_historical',
#                        lineThicknesses= thicknesses)
		#assert 0 
                jobs = ['u-at643', 'u-at629', 'u-at646','u-ar977']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='Drift2',
                        jobDescriptions=jobDescriptions,
                        analysisname='Set3_UM10.7_STRATROP_plusOceanOnly',
                        lineThicknesses= thicknesses)


#                jobs = ['u-at572','u-at482',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='Set1_UM10.7_dmsmin1.25',
#                        lineThicknesses= thicknesses)

#                jobs = ['u-at572','u-at574',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='Set2_UM10.7_CN',
#                        lineThicknesses= thicknesses)

#                jobs = ['u-at643', 'u-at629', 'u-at646',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='Set3_UM10.7_STRATROP',
#                        lineThicknesses= thicknesses)
                        
                jobs = ['u-ar977','u-as462','u-as643',]
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='Drift',
                        jobDescriptions=jobDescriptions,
                        analysisname='DriftTest',
                        lineThicknesses= thicknesses)

                jobs = ['u-ar783','u-ar855','u-ar799','u-ar835']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0=True,
                        jobDescriptions=jobDescriptions,
                        analysisname='StratTropTests',
                        lineThicknesses= thicknesses)


#                jobs = ['u-an869','u-ao586']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='OOvFC1',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='OceanOnlyVsFullCoupled_1',
#                        lineThicknesses= thicknesses)

                jobs = ['u-ar538','u-ar783']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='OOvFC2',
                        jobDescriptions=jobDescriptions,
                        analysisname='OceanOnlyVsFullCoupled_2',
                        lineThicknesses= thicknesses)

                jobs = ['u-an869','u-ao586','u-ak900','u-ar538','u-ar783']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='FullSpinUp',
                        jobDescriptions=jobDescriptions,
                        analysisname='OriginalOceanOnlySpinUp_debug',
                        lineThicknesses= thicknesses)

#                jobs = ['u-an869','u-ao586','u-ak900','u-ar538']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='FullSpinUp',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='OriginalOceanOnlySpinUp')

                jobs = ['u-ar783','u-ar538','u-am927i','u-am927ii','u-am927iii','u-aq853iii',]
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='4945-5110',
                        jobDescriptions=jobDescriptions,
                        analysisname='UKESM_0.9.2_from4945',
                        lineThicknesses= thicknesses)

                jobs = ['u-ar977','u-ar538','u-ar975',]
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='from4950',
                        jobDescriptions=jobDescriptions,
                        analysisname='UKESM_NOC_runs',)



#                jobs = ['u-ar538','u-ar748','u-ar773','u-ar761','u-ar769', 'u-ar976','u-ar977']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='4830-5000',#'UKESM_0.9.1',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM_0.9.1_NewSuiteDebugging')


                jobs = ['u-ar538','u-ar976','u-ar977']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0=False,#'UKESM_0.9.1',
                        jobDescriptions=jobDescriptions,
                        analysisname='UKESM_0.9.1_longerTermSuff')

                #jobs = ['u-an869','u-ar783','u-ar538','u-am064','u-am927','u-aq853']
                #colours = {i:standards[i] for i in jobs}
                #timeseries_compare({
                #        i:standards[i] for i in jobs},
                #        physics=1,
                #        bio=1,
                #        debug=0,
                #        year0=False, #UKESM_0.9.1',
                #        jobDescriptions=jobDescriptions,
                #        analysisname='UKESM_0.9.2')
		
#                jobs = ['u-ar951','u-ar379',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='True',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='IceburgVsNoIceburgFix')


#                jobs = ['u-aq853','u-ar549',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='True',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='CoupledRunwithOOIC_vs_CoupledRun')

#                jobs = ['u-ar538','u-ar748','u-ar773','u-ar761','u-ar769', 'u-ar976','u-ar977']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='4830-5000',#'UKESM_0.9.1',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM_0.9.1_NewSuiteDebugging')
#



                #jobs = ['u-an869','u-ar537','u-ar538','u-ar541', 'u-ar615','u-ar616' ]
                #colours = {i:standards[i] for i in jobs}
                #timeseries_compare({
                #        i:standards[i] for i in jobs},
                #        physics=1,
                #        bio=1,
                #        debug=0,
                #        year0='ignoreStart',
                #        jobDescriptions=jobDescriptions,
                #        analysisname='u-an869_debugging')		

                jobs = ['u-an869','u-ao586','u-ak900','u-ar538']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0='FullSpinUp',
                        jobDescriptions=jobDescriptions,
                        analysisname='OriginalOceanOnlySpinUp')

#                jobs = ['u-ao912','u-ap865','u-aq142','u-aq710']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM0.8-PrescribedVeg')

                jobs = ['u-am064','u-am515','u-aq378','u-am927','u-aq853']
                colours = {i:standards[i] for i in jobs}
                timeseries_compare({
                        i:standards[i] for i in jobs},
                        physics=1,
                        bio=1,
                        debug=0,
                        year0=False,
                        jobDescriptions=jobDescriptions,
                        analysisname='CMIP6_spin-up_runs_GC3.1_UM10.7')

#                jobs = ['u-ao365','u-am004','u-ao837','u-aq226']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='True',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM_0.9')
#		#assert 0

#                jobs = ['u-ap951','u-ao586',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='True',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='u-ap951_vs_u-ao586')

#                jobs = ['u-an869','u-ao586','u-ap951',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM_spinUp_Updated')

#                jobs = ['u-ap951','u-an869',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='u-ap951_vs_u-an869')

		

#                jobs = ['u-ao912','u-ap450',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM0.8-PMOA')

#                jobs = ['u-an869','u-ap951',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM0.8-MEDUSA')


#                jobs = ['u-ao912','u-ao913','u-ao914','u-ao884','u-ao886',]#'u-ao837',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='LandChlTests')

#                jobs = ['u-ao365','u-am004','u-ao837',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='UKESM0.8',
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM_0.8')

#                jobs = ['u-an869','u-ao586','u-ak900',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0='FullSpinUp', 
#                        jobDescriptions=jobDescriptions,
#                        analysisname='OriginalOceanOnlySpinUp')

#                jobs = ['u-an869','u-an908','u-an911','u-an631','u-an629','u-an619','u-an989',]
#                timeseries_compare(
#                        {i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#			analysisname='AllOceanOnlySpinUps')

#                jobs = [ 'u-an989','u-an908','u-an619']
#                timeseries_compare(
#                	{i:standards[i] for i in jobs},
#                	physics=1,
#                	bio=1,
#                	debug=0,
#                	year0=False,
#                	jobDescriptions=jobDescriptions,
#                	analysisname='IC_Observations')
                
                
#                jobs = ['u-an631','u-an869','u-an629','u-an911']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                	i:standards[i] for i in jobs},
#                	physics=1,
#                	bio=1,
#                	debug=0,                	
#                	year0=False,
#                	jobDescriptions=jobDescriptions,
#                	analysisname='IC_u-ak900')
                                

#                jobs = ['u-am927','u-an619','u-an629','u-an631',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='atmos_u-am927')

#                jobs = ['u-am515','u-an908','u-an869','u-an911','u-an989']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='atmos_u-am515')

#                jobs = ['u-ao365','u-ao404',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare({
#                        i:standards[i] for i in jobs},
#                        physics=1,
#                        bio=1,
#                        debug=0,
#                        year0=False,
#                        jobDescriptions=jobDescriptions,
#                        analysisname='UKESM_0.7')
						

#                jobs = ['u-am001','u-am004',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare(colours, physics=1,bio=0,year0=False,debug=1,analysisname='UKESM0.6_vs_GC3.1')
		
#                jobs = ['u-ak900','u-an631','u-an629','u-an619',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare(colours, physics=1,bio=1,year0='2000-2600normu-ak900',debug=0,jobDescriptions=jobDescriptions, analysisname='OceanOnlySpinUp_u-ak900')


#                jobs = ['u-ai611','u-aj391','u-al901','u-am064','u-am927','u-am515',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare(colours, physics=1,bio=0,year0=False,debug=0,jobDescriptions=jobDescriptions, analysisname='TillsCoupledRuns_physics')
                
#		jobs = ['u-am927','u-am515',]
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare(colours, physics=0,bio=1,year0=True,debug=0,jobDescriptions=jobDescriptions, analysisname='TillsCoupledRuns_BGC')


	
#                jobs = ['u-ai567','u-aj588','u-am696','u-am792','u-am892','u-ak900']
#                colours = {i:standards[i] for i in jobs}
#                timeseries_compare(colours, physics=1,bio=1,year0='First100Years',debug=0,jobDescriptions=jobDescriptions, analysisname='CirculationTests-june17')



		





	

if __name__=="__main__": main()


