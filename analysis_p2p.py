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

#Standard Python modules:
from sys import argv,exit
from os.path import exists
from calendar import month_name
from socket import gethostname
from glob import glob
from netCDF4 import Dataset
import numpy as np


#Specific local code:
import UKESMpython as ukp
from p2p import makePatternStatsPlots, testsuite_p2p
from p2p.patternAnalyses import InterAnnualPatterns,BGCvsPhysics
from pftnames import months
from p2p.shelveToDictionary import shelveToDictionary
#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	

# 	from 




def analysis_jasmin(
		models	= ['NEMO','MEDUSA',],
		jobID 	= 'xkrus',
		years 	= ['2077'], #'2075','2076',
		modelGrid = 'ORCA1',
		annual 	= False,
		noPlots = False,
		analysisSuite='bio',
		):
	
	# DMS model:
	#model= 'MEDUSA'
	#jobID = 'xkrum'
	#year = 'clim'		
	#modelGrid = 'ORCA1'
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobID+'-' + year+'/'

	# ORCA1:
	#models= ['MEDUSA','NEMO']
	#jobID = 'xkrus'
	#years = ['2077-annual'] #'2075','2076',
	#modelGrid = 'ORCA1'
	#MEDUSAFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	#NEMOFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	#annual = True


	
	# ORCA025:
	#model= 'MEDUSA'
	#jobID = 'xjwki'
	#year = '1979'		
	#modelGrid = 'ORCA025'	
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/"+jobID+"_postProc/"+year+"/"
	
	
	
	#####
	# Which analysis to run
	
	if analysisSuite.lower() in ['default',]:
		doCHL 		= True
		doCHL_CCI 	= 0				
		doDiatoms	= 0#True
		doMicrozoo	= 0#True
		doMesozoo	= 0#True
		doN		= True
		doSi		= True	
		doFe		= True		
		doPCO2		= 0#True
		doAlk		= True			# Glodap Alkalinity
		doDIC		= True			# Globap tCO2
		doAirSeaFlux	= True			# work in progress
		doIntPP_iMarNet	= True			# Integrated primpary production from iMarNEt
		doIntPP_OSU	= True			# OSU Integrated primpary production	
		doO2		= True	
		doSal		= True
		doTemp		= True
		doMLD		= True
	elif analysisSuite.lower() in ['debug',]:
		doCHL 		= 0#True
		doCHL_CCI 	= 0				
		doDiatoms	= 0#True
		doMicrozoo	= 0#True
		doMesozoo	= 0#True		
		doN		= 0#True
		doSi		= 0#True	
		doFe		= 0#True		
		doPCO2		= 0#True
		doO2		= 0#True	
		doSal		= 0#True
		doTemp		= 0#True
		doMLD		= 0#True
		doAlk		= 0#True			# Glodap Alkalinity
		doDIC		= 0#True			# Globap tCO2
		doAirSeaFlux	= 0#True		# work in progress
		doIntPP_iMarNet	= 0#True		# Integrated primpary production from iMarNEt
		doIntPP_OSU	= True		# OSU Integrated primpary production	
	elif analysisSuite.lower() in ['bio',]:
		doCHL 		= True
		doCHL_CCI 	= 0				
		doDiatoms	= True
		doMicrozoo	= True
		doMesozoo	= True		
		doN		= 0#True
		doSi		= 0#True	
		doFe		= 0#True		
		doPCO2		= 0#True
		doO2		= 0#True	
		doSal		= 0#True
		doTemp		= 0#True
		doMLD		= 0#True
		doAlk		= True			# Glodap Alkalinity
		doDIC		= True			# Globap tCO2
		doAirSeaFlux	= True			# work in progress
		doIntPP_iMarNet	= True			# Integrated primpary production from iMarNEt
		doIntPP_OSU	= True			# OSU Integrated primpary production	
	elif analysisSuite.lower() in ['annual',]:
		#####
		# Data sets with annual coverage available (instead o monthly)
		doCHL 		= 0
		doCHL_CCI 	= True		
		doDiatoms	= 0
		doMicrozoo	= 0
		doMesozoo	= 0		
		doN		= True
		doSi		= True	
		doFe		= 0#True		
		doPCO2		= 0#True
		doO2		= True	
		doSal		= True
		doTemp		= True
		doMLD		= True
		doAlk		= True			# Glodap Alkalinity
		doDIC		= True			# Globap tCO2
		doAirSeaFlux	= True			# work in progress
		doIntPP_iMarNet	= True			# Integrated primpary production from iMarNEt
		doIntPP_OSU	= True			# OSU Integrated primpary production	
	else:		
		doCHL 		= 0#True
		doCHL_CCI 	= 0				
		doDiatoms	= 0#True
		doMicrozoo	= 0#True
		doMesozoo	= 0#True
		doN		= True
		doSi		= 0#True	
		doFe		= 0#True
		doPCO2		= 0#True
		doO2		= 0#True	
		doSal		= 0#True
		doTemp		= 0#True
		doMLD		= 0#True
		doAlk		= True			# Glodap Alkalinity
		doDIC		= True			# Globap tCO2
		doAirSeaFlux	= True			# work in progress
		doIntPP_iMarNet	= True			# Integrated primpary production from iMarNEt
		doIntPP_OSU	= True			# OSU Integrated primpary production	
	
	
	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)	
	#if annual:	
	#	print "WARNING: annual data not yet tested."
	#	assert False
	#if annual:depthLevels 	= ['',]
	#else:		depthLevels 	= ['Transect','Surface','100m','200m','500m',]


	#####
	# Location of data files.
	if gethostname().find('pmpc')>-1:	
		print "analysis-JASMIN.py:\tBeing run at PML on ",gethostname()
		
		if annual:	
			#####
			# No need to stitch together multiple months into one file:
			WOAFolder 		= "/data/euryale7/scratch/ledm/WOA/annual/"
			MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
			NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"			
		else:		
			WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
			MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"_postProc/"
			NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"_postProc/"			
		#MAREDATFolder 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"
		#GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
		#TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
		#MLDFolder	= "/data/euryale7/scratch/ledm/IFREMER-MLD/"
		
		ObsFolder = "/data/euryale7/backup/ledm/Observations/"
		MAREDATFolder 	= ObsFolder+"/MAREDAT/MAREDAT/"
		GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
		MLDFolder	= ObsFolder+"/IFREMER-MLD/"
		iMarNetFolder	= ObsFolder+"/LestersReportData/"
		GlodapDir	= ObsFolder+"/GLODAP/"
		GLODAPv2Dir	= ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
		OSUDir		= ObsFolder+"OSU/"
		CCIDir		= ObsFolder+"CCI/"
				
		if jobID in ["xkrus",]:
			# Old school ORCA1 grid
			orcaGridfn 	='data/mesh_mask_ORCA1_75.nc'
		else:
			# New eORCA1 grid		
			orcaGridfn 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
				
		workDir		= "/data/euryale7/scratch/ledm/ukesm_postProcessed/"
		imgDir		= ukp.folder('images')
		
	if gethostname().find('ceda.ac.uk')>-1:
		print "analysis-JASMIN.py:\tBeing run at CEDA on ",gethostname()
			
		ObsFolder 	= "/group_workspaces/jasmin/esmeval/example_data/bgc/"
		esmvalFolder 	= "/group_workspaces/jasmin/esmeval/data/"		
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
		
	
		# Directory for output files:
		workDir 	= ukp.folder(esmvalFolder+"ukesm_postProcessed/")
		imgDir		= ukp.folder('images')		

		if jobID in ["xkrus",]:
			# Old school ORCA1 grid
			orcaGridfn 	='/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_ORCA1_75.nc'
		else:
			# New eORCA1 grid		
			orcaGridfn 	= '/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'
			
						
	#####
	# Set which spatial and temporal limitations to plot.
	transects 	= ['AtlanticTransect', 'PacificTransect',]
	justAll		= ['All',]				# All is not a slice, it has no cut on location, time, or depth.
	AllStandard	= ['All','Standard','ignoreInlandSeas']	
	HighLatWinter	= ['All','HighLatWinter',]
	tsRegions	= ['Global','Equator10', 'Remainder','ArcticOcean','NorthernSubpolarAtlantic','NorthernSubpolarPacific','ignoreInlandSeas','SouthernOcean',]
					
	medusaCoords 	= {'t':'index_t', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}		
	cciCoords	= {'t':'index_t', 'z':'index_z','lat': 'lat',      'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	glodapCoords	= {'t':'index_t', 'z':'depth',  'lat': 'latitude', 'lon': 'longitude', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero'] }	
	osuCoords	= {'t':'index_t', 'z':'index_z','lat': 'latitude', 'lon': 'longitude', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero'] }		
	glodapv2Coords	= {'t':'index_t', 'z':'Pressure','lat':'lat',      'lon': 'lon',        'cal': '',        'tdict':{0:0,} }	
	shelvesAV = []	
	
	for year in years:		
		#####
		# Location of model files.
		if annual:
			MEDUSAFolder 	= MEDUSAFolder_pref+jobID+"/"
			NEMOFolder 	= NEMOFolder_pref+jobID+"/"			
		else:			
			MEDUSAFolder 	= MEDUSAFolder_pref+year+'/'
			NEMOFolder  	= NEMOFolder_pref+year+'/'		
		

		#####
		# AutoVivification is a form of nested dictionary.
		# We use AutoVivification here to determine which files to analyse and which fields in those files.
		# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
		av = ukp.AutoVivification()
		if doCHL:
			name = 'Chlorophyll_pig'		
			av[name]['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
			if modelGrid == 'ORCA1':	av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_CHL.nc"
			if modelGrid == 'ORCA025':	av[name]['MEDUSA']['File']	= MEDUSAFolder+"xjwki_1979_CH.nc"
			
			av[name]['Data']['coords'] 		= maredatCoords
			av[name]['MEDUSA']['coords']		= medusaCoords
			
			av[name]['MEDUSA']['details']		= {'name': 'CHL', 'vars':['CHL',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			
			av[name]['Data']['details']		= {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}			

			av[name]['Data']['source'] 		= 'MAREDAT'
			av[name]['MEDUSA']['source']		= 'MEDUSA'

			av[name]['depthLevels'] 		= ['',]
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['plottingSlices'] 		= tsRegions

		if doCHL_CCI:						
			name = 'Chlorophyll_cci'
			if annual:
				av[name]['Data']['File'] 	= CCIDir+"ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-annual-fv2.0.nc"	
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]					
			else:
				av[name]['Data']['File'] 	= CCIDir+'ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-all-fv2.0.nc'
				av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_CHL.nc"
			
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['',]
			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= cciCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'CCI'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']		= {'name': name, 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}		
			av[name]['Data']['details']		= {'name': name, 'vars':['chlor_a',], 'convert':  ukp.NoChange,'units':'mg C/m^3'}			
		
		if doDiatoms:
			name = 'Diatoms'
			if annual: 
				print "No diatoms iron file",
				assert 0		
			av[name]['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
			av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_PHD.nc"
			
			av[name]['depthLevels'] 		= ['',]	
			av[name]['MEDUSA']['grid']		= modelGrid						
			av[name]['plottingSlices']		= AllStandard

			av[name]['Data']['coords'] 	= maredatCoords
			av[name]['MEDUSA']['coords']	= medusaCoords

			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['PHD',],     'convert': ukp.N2Biomass,'units': 'mg C/m^3'}						
			av[name]['Data']['details']	= {'name': name, 'vars':['BIOMASS',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			

			av[name]['Data']['source'] 	= 'MAREDAT'
			av[name]['MEDUSA']['source']	= 'MEDUSA'

			
		if doMicrozoo:
			name = 'Microzoo'		
			if annual: 
				print "No microzoo iron file",
				assert 0		
			av[name]['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
			av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_ZMI.nc"	
			
			av[name]['MEDUSA']['grid']	= modelGrid		
			av[name]['depthLevels'] 		= ['',]	
			av[name]['plottingSlices'] 	= AllStandard

			av[name]['Data']['coords'] 	= maredatCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['ZMI',],     'convert': ukp.N2Biomass,'units': 'mg C/m^3'}			
			av[name]['Data']['details']	= {'name': name, 'vars':['BIOMASS',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			

			av[name]['Data']['source'] 	= 'MAREDAT'
			av[name]['MEDUSA']['source']	= 'MEDUSA'
			
			
		if doMesozoo:
			name = 'Mesozoo'
			if annual: 
				print "No mesozoo iron file",
				assert 0
						
			av[name]['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
			av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_ZME.nc"	

			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['',]
			av[name]['plottingSlices'] 	= AllStandard

			av[name]['Data']['coords'] 	= maredatCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['ZME',],     'convert': ukp.N2Biomass,'units': 'mg C/m^3'}			
			av[name]['Data']['details']	= {'name': name, 'vars':['BIOMASS',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			

			av[name]['Data']['source'] 	= 'MAREDAT'
			av[name]['MEDUSA']['source']	= 'MEDUSA'
			
			
		if doN:
			name = 'Nitrate'		
			if annual:	
				av[name]['Data']['File'] 	= WOAFolder+'woa13_all_n00_01.nc'
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]	
			else:		
				av[name]['Data']['File'] 	= WOAFolder+'nitrate_monthly_1deg.nc'	
				if modelGrid == 'ORCA1':	av[name]['MEDUSA']['File'] = MEDUSAFolder+jobID+'_' + year+"_DIN.nc"	
				if modelGrid == 'ORCA025':	av[name]['MEDUSA']['File'] = MEDUSAFolder+jobID+'_'+ year+"_DIN.nc"							
			
			av[name]['MEDUSA']['grid']	= modelGrid		
			av[name]['depthLevels'] 	= ['Surface','Transect','PTransect','SOTransect',]
			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['DIN',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['n_an',], 'convert': ukp.NoChange,}	# no units?

						
		if doSi:
			name = 'Silicate'
			if annual:
				av[name]['Data']['File'] 	= WOAFolder+'woa13_all_i00_01.nc'
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]					
			else:
				av[name]['Data']['File'] 	= WOAFolder+'silicate_monthly_1deg.nc'	
				av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_SIL.nc"
			
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']
			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['SIL',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['i_an',], 'convert': ukp.NoChange,}	# no units?
			
						
		if doFe:
			name = 'Iron'
			if annual: 
				print "No annual iron file",
				assert 0
				
			av[name]['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
			av[name]['MEDUSA']['File'] 		= MEDUSAFolder+jobID+'_' + year+"_FER.nc"	
			
			av[name]['depthLevels'] 		= ['',]
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['plottingSlices']		= justAll
			
			av[name]['Data']['coords'] 	= {'t': 'MONTH','z':'DEPTH','lat':'Latitude','lon':'Longitude','cal':'standard','tdict': ukp.tdicts['OneToZero']}
			av[name]['MEDUSA']['coords']	= medusaCoords
	
			
			av[name]['Data']['source'] 	= 'GEOTRACES'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['FER',], 'convert': ukp.mul1000,'units':'umol F/m^3'}			
			av[name]['Data']['details']	= {'name': name, 'vars':['Fe_D_CONC_BOTTLE',], 'convert': ukp.NoChange,}	# no units?

			
		if doO2:
			name = 'Oxygen'		
			if annual:
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]					
				av[name]['Data']['File'] 	=  WOAFolder+'woa13_all_o00_01.nc'
			else:	
				av[name]['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
				av[name]['MEDUSA']['File']	= MEDUSAFolder+jobID+"_"+year+"_OXY.nc"
			
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']
			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['OXY',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']		= {'name': name, 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol/m^3'}


		if doAlk:
			name = 'Alkalinity'	
			def convertmeqm3TOumolkg(nc,keys):
				return nc.variables[keys[0]][:]* 1.027
		
			if annual:		
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]
				av[name]['Data']['File'] 		=  GlodapDir+'Alk.nc'
			else:
				print "Alkalinity data not available for monthly Analysis"
				assert 0
				
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']
				

			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= glodapCoords
			av[name]['MEDUSA']['coords']	= medusaCoords				

			av[name]['Data']['source'] 	= 'GLODAP'
			av[name]['MEDUSA']['source']	= 'MEDUSA'

			av[name]['MEDUSA']['details'] 	= {'name': name, 'vars':['ALK',], 'convert': ukp.NoChange,'units':'meq/m^3',}
			av[name]['Data']['details']  	= {'name': name, 'vars':['Alk',], 'convert': convertmeqm3TOumolkg,'units':'meq/m^3',}		
						
		if doDIC:
			name = 'DIC'
		
			if annual:		
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]
				av[name]['Data']['File'] 	=  GLODAPv2Dir+'GLODAPv2.tco2.historic.nc'
			else:
				print "DIC data not available for monthly Analysis"
				assert 0
				
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']
				

			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= glodapv2Coords
			av[name]['MEDUSA']['coords']	= medusaCoords				

			av[name]['Data']['source'] 	= 'GLODAPv2'
			av[name]['MEDUSA']['source']	= 'MEDUSA'

			av[name]['MEDUSA']['details'] 	= {'name': name, 'vars':['DIC',],  'convert': ukp.NoChange,'units':'mmol C/m^3'}
			av[name]['Data']['details']  	= {'name': name, 'vars':['tco2',], 'convert': ukp.NoChange,'units':'micro-mol kg-1'}
	
		

		if doIntPP_OSU:
			name = 'IntegratedPrimaryProduction_OSU'
			
			
			#####
			# Files:
			if annual:
				av[name]['MEDUSA']['File']  	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_diad_T.nc"))[0]				
				av[name]['Data']['File']  		= OSUDir +"/standard_VGPM.SeaWIFS.global.average.nc"
			else:
				print "IntegratedPrimaryProduction (OSU) data not available for monthly Analysis"
				assert 0
							
			#####
			# Calculating depth in PP in medusa			
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
				elif arr.ndim==1:
					index_x = nc.variables['index_x'][:]
					index_y = nc.variables['index_y'][:]
					for i,a in enumerate(arr):
						arr[i] = a * area[index_y[i],index_x[i]]
				else: assert 0
				return arr


			#####
			# converting data to same units.						
			nc = Dataset(av[name]['Data']['File'] ,'r')
			lats = nc.variables['latitude'][:]
			osuareas = np.zeros((1080, 2160))
			osuarea = (111100. / 6.)**2. # area of a pixel at equator. in m2
			for a in np.arange(1080):osuareas[a] = np.ones((2160,))*osuarea*np.cos(np.deg2rad(lats[a]))
		
			def osuconvert(nc,keys):
				arr = nc.variables[keys[0]][:] 
				tlen = 1 # arr.shape[0]
				arr  = arr/tlen * 365.	/ 1000. /     1E15
				if arr.ndim ==3:
					for i in np.arange(arr.shape[0]):
						arr[i] = arr[i]*osuareas
				elif arr.ndim ==2: arr = arr*osuareas
				elif arr.ndim ==1:
					index_x = nc.variables['index_x'][:]
					index_y = nc.variables['index_y'][:]
					for i,a in enumerate(arr):
						print i,a,[index_y[i],index_x[i]]
						arr[i] = a * osuareas[index_y[i],index_x[i]]				
				else: 
					assert 0
				return arr
						
			av[name]['MEDUSA']['coords'] 	= medusaCoords 	
			av[name]['Data']['coords']	= osuCoords

			av[name]['MEDUSA']['grid']	= modelGrid		
			av[name]['depthLevels'] 	= ['',]
			
			if annual:	av[name]['plottingSlices'] 	= tsRegions
			else:		av[name]['plottingSlices'] 	= HighLatWinter
						
			av[name]['Data']['source'] 	= 'OSU'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
			
			av[name]['MEDUSA']['details'] 	= {'name': name, 'vars':['PRN' ,'PRD'], 'convert': medusadepthInt,'units':'gC/yr'}
			av[name]['Data']['details']  	= {'name': name, 'vars':['NPP',], 'convert': osuconvert,'units':'gC/yr'}
					
						
		

			
					

		if doAirSeaFlux:
			assert 0
			
			
	
						
					
#		#if doPCO2:
#			av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
#			av['pCO2']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"
#			av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
#			av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]	
#			av['pCO2']['depthLevels'] 	= ['',]
#			av['pCO2']['MEDUSA']['grid']	= modelGrid				
#			#av['pCO2']['plottingSlices'] 	= []
#
#	#	if doIntPP:
#	#		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
#	#		av['intpp']['Data']['Vars'] 	= ['PPint',]
#	#		#av['intpp']['ERSEM']['Vars'] 	= ['netPP',] This field is net, not integrated.
#	#		#av['intpp']['ERSEM']['File']	= ERSEMFolder+'_ERSEMMisc.nc' This file will need to be remade.
#	#		av['intpp']['ERSEM']['grid']	= modelGrid
#	#		av['intpp']['depthLevels'] 		= ['',]
		
		
		if doSal:
			name = 'Salinity'
			if annual:
				av[name]['NEMO']['File'] 	= sorted(glob(NEMOFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_grid_T.nc"))[0]
				av[name]['Data']['File'] 	= WOAFolder+'woa13_decav_s00_01v2.nc'	
			else:	
				av[name]['Data']['File'] 	= WOAFolder+'salinity_monthly_1deg.nc'	
				av[name]['NEMO']['File'] 	= NEMOFolder+jobID+"_"+year+'_SAL.nc'	

			av[name]['NEMO']['grid'] 		= modelGrid
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']	 
			av[name]['plottingSlices'] 	= tsRegions
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['NEMO']['coords']	= medusaCoords
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['NEMO']['source']	= 'NEMO'			
	
			av[name]['NEMO']['details']	= {'name': name, 'vars':['vosaline',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['s_an',], 'convert': ukp.NoChange,}	# no units?
			
					
		if doTemp:
			name = 'Temperature'
			if annual:
				av[name]['NEMO']['File'] 	= sorted(glob(NEMOFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_grid_T.nc"))[0]				
				av[name]['Data']['File'] 	= WOAFolder+'woa13_decav_t00_01v2.nc'	
			else:	
				av[name]['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
				av[name]['NEMO']['File'] 	= NEMOFolder+jobID+"_"+year+'_TEMP.nc'	

			av[name]['NEMO']['grid'] 	= modelGrid	
			av[name]['depthLevels'] 	= ['Surface','Transect','PTransect']	
			av[name]['plottingSlices'] 	= tsRegions

			av[name]['Data']['coords'] 	= woaCoords
			av[name]['NEMO']['coords']	= medusaCoords
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['NEMO']['source']	= 'NEMO'			
	
			av[name]['NEMO']['details']	= {'name': name, 'vars':['votemper',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,}	# no units?
			
						   
		if doMLD:
			name = 'MLD'		
			if annual:	
				av[name]['NEMO']['File'] 	= sorted(glob(NEMOFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_grid_T.nc"))[0]							
				av[name]['Data']['File'] 		= MLDFolder+"mld_DT02_c1m_reg2.0-annual.nc"
			else:	
				av[name]['Data']['File'] 		= MLDFolder+"mld_DT02_c1m_reg2.0.nc"
				av[name]['NEMO']['File'] 		= NEMOFolder+jobID+"_"+year+'_MLD.nc'	
					
			av[name]['NEMO']['grid'] 		= modelGrid
			av[name]['depthLevels'] 		= ['',]
			av[name]['plottingSlices'] 		= tsRegions

			av[name]['Data']['coords'] 	= {'t':'index_t', 'z':'index_z','lat':'lat','lon':'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
			av[name]['NEMO']['coords']	= medusaCoords
			av[name]['Data']['source'] 	= 'IFREMER'
			av[name]['NEMO']['source']	= 'NEMO'			
	
			av[name]['NEMO']['details']	= {'name': name, 'vars':['somxl010',], 'convert': ukp.NoChange,'units':'m'}			
			av[name]['Data']['details']	= {'name': name, 'vars':['mld','mask',], 'convert': ukp.applymask,'units':'m'}	# no units?
			
		
		for model in models:
			workingDir 	= ukp.folder(workDir+model+'-'+jobID+'-'+year)
			imageFolder 	= ukp.folder(imgDir+'/'+jobID)

			shelvesAV.extend(
		    		testsuite_p2p(
					model 		= model,
					jobID 		= jobID,
					year  		= year,
					av 		= av,
					plottingSlices	= [],	# set this so that testsuite_p2p reads the slice list from the av.
					workingDir 	= workingDir,
					imageFolder	= imageFolder,
					noPlots		= noPlots,	# turns off plot making to save space and compute time.
					gridFile	= orcaGridfn,	# enforces custom gridfile.
					annual		= annual,
			 	)
			)
	
		
	#BGCvsPhysics(shelvesAV, jobID, modelGrid )
	#if len(years)>1: InterAnnualPatterns(shelvesAV, jobID, years,modelGrid)	# plots interannual comparison and time series.
#	def outPutForJASMIN(shelvesAV):
#	outdict = shelveToDictionary(shelvesAV)
	
		
if __name__=="__main__":
	#analysis_jasmin()
	try: 	jobID = argv[1]
	except:	jobID =	'u-ab749'
	try:	year = argv[2]
	except:	year = '2007'
	analysis_jasmin(models	= ['NEMO','MEDUSA',],
		jobID 	= jobID,
		years 	= [year,], #'2075','2076',
		modelGrid = 'eORCA1',
		annual 	= True,
		noPlots = False,
		analysisSuite='debug',)
	
			

