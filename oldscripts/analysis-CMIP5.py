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

#Specific local code:
from UKESMpython import populateSlicesList, AutoVivification, folder, reducesShelves,listShelvesContents,mnStr
from testsuite_p2p import testsuite_p2p
from p2p import makePatternStatsPlots
from p2p.patternAnalyses import InterAnnualPatterns,BGCvsPhysics,modelIntercomparisonAnnual
from bgcvaltools.pftnames import months


#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	

# 	from 

def analysis():
	# DMS model:
	#model= 'MEDUSA'
	#jobID = 'xkrum'
	#year = 'clim'		
	#modelGrid = 'ORCA1'
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobID+'-' + year+'/'

	# ORCA1:
	CMIP5models = [ 'HadGEM2-ES','ERSEM', 'MEDUSA',
		'GFDL-ESM2G', 'GFDL-ESM2M', 'CESM1-BGC', 'NorESM1-ME',
		'MPI-ESM-MR', 'IPSL-CM5A-MR', 'MRI-ESM1', # 'HadGEM2-ES',
		'MPI-ESM-LR', 'HadGEM2-CC',  'IPSL-CM5A-LR',
		'IPSL-CM5B-LR',  'CMCC-CESM',  'CNRM-CM5', 'BNU-ESM',	# These 4 didn't make it into the A Cabre paper.
		]
		
					
	#models= ['MEDUSA','NEMO']
	jobID = 'CMIP5'
	years = ['2005-annual'] #'2075','2076',
	modelGrid = 'Flat1deg'
	CMIP5Folder= "/data/euryale7/scratch/ledm/CMIP5_postProcessed/"
	#NEMOFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	annual = True
	

		
	#####
	# Which analysis to run
	doCHL 		= 0#True
	doMAREDAT 	= 0#True
	doN		= 0#True
	doSi		= 0#True	
	doFe		= 0#True		
	doPCO2		= 0#True
	doIntPP		= 0#True
	doO2		= 0#True
		
	doU		= 0#True
	doV		= True	
	doSal		= 0#True
	doTemp		= 0#True/data/euryale7/scratch/ledm/CMIP5_postProcessed/O2/cmip5_HadGEM2-ES_2005.nc
	doMLD		= 0#True
	
	
	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)	
	if annual:	depthLevels 	= ['',] #'2000m','Surface',
	else:		depthLevels 	= ['Transect','Surface','100m','200m','500m',]

	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	if annual:	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/annual/"
	else:		WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GODASFolder	= "/data/euryale7/scratch/ledm/GODAS_postProcessed/"
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
	LesterFolder 	= "/data/euryale7/scratch/ledm/LestersReportData/"			
	
	shelvesAV = []
	
	for year in years:
	    for model in CMIP5models[:]:		
		#####
		# Location of model files.
		#modelFolder 	= CMIP5Folder_pref+year+'/'
		#NEMOFolder  	= NEMOFolder_pref+year+'/'			
		
		#####
		# AutoVivification is a form of nested dictionary.
		# We use AutoVivification here to determine which files to analyse and which fields in those files.
		# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
		av = AutoVivification()


		if doO2:
			if annual:	av['oxygen']['Data']['File'] 	=  WOAFolder+'woa13_all_o00_01.nc'
			#if annual:	av['oxygen']['Data']['File'] 	=  "/data/euryale7/scratch/ledm/WOA/annual/tmp/woa13_all_o00_01_lowhalf.nc"
			else:		av['oxygen']['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
			av['oxygen']['Data']['Vars'] 	= ['o_an',] 
			av['oxygen']['CMIP5_'+model]['Vars'] 	= ['o2',]
			av['oxygen']['CMIP5_'+model]['File']	= CMIP5Folder+"O2/"+year+"/cmip5_"+model+"_2005.nc"
			av['oxygen']['CMIP5_'+model]['grid']	= modelGrid
			av['oxygen']['depthLevels'] 		= depthLevels
			gridFile = av['oxygen']['CMIP5_'+model]['File']
		
		if doU:
			if annual:	av['U']['Data']['File'] 	=  GODASFolder+'ucur.2005.nc'
			else:	assert False
			av['U']['Data']['Vars'] 	= ['uo',] 
			av['U']['CMIP5_'+model]['Vars'] = ['uo',]
			av['U']['CMIP5_'+model]['File']	= CMIP5Folder+"U/uo_"+model+"_2005.nc"
			av['U']['CMIP5_'+model]['grid']	= modelGrid
			av['U']['depthLevels'] 		= depthLevels
			gridFile = av['U']['CMIP5_'+model]['File']

		if doV:
			if annual:	av['V']['Data']['File'] 	=  GODASFolder+'vcur.2005.nc'
			else:	assert False
			av['V']['Data']['Vars'] 	= ['vo',] 
			av['V']['CMIP5_'+model]['Vars'] = ['vo',]
			av['V']['CMIP5_'+model]['File']	= CMIP5Folder+"V/vo_"+model+"_2005.nc"
			av['V']['CMIP5_'+model]['grid']	= modelGrid
			av['V']['depthLevels'] 		= depthLevels
			gridFile = av['V']['CMIP5_'+model]['File']
			
						
				
		#####
		# Set which spatial and temporal limitations to plot.
		if annual:
			plottingSlices = populateSlicesList(
				 plotDefaults		=True,	
				 plotMonths		=0,#True,
				 plotdepthRanges	=True,
				 plotpercentiles	=0,#True	
				 plotLatRegions		=0,# True
				 plotQualityCuts	=0,#True	
				 plotSeas		=0,#True		 
				 plotOceans		=True,
				 plotHemispheres	=0,#True,
				 plotSeasons		=0,# True
				 plotOceanSeasons	=0,# True		 		 
				 plotOceanMonths   	=0,#True	
				 plotHemispheresMonths  =0,#True,
				 plotTransects		= True
			)
				
		else:
			plottingSlices = populateSlicesList(
				 plotDefaults		=True,	
				 plotMonths		=True,
				 plotdepthRanges	=0,
				 plotpercentiles	=0,#True	
				 plotLatRegions		=0,# True
				 plotQualityCuts	=0,#True	
				 plotSeas		=0,#True		 
				 plotOceans		=True,
				 plotHemispheres	=0,# True
				 plotSeasons		=0,# True
				 plotOceanSeasons	=0,# True		 		 
				 plotOceanMonths   	=0,#True	
				 plotHemispheresMonths  =True,
				 plotTransects		= True				 
			)

		#for model in models:
#		workingDir = folder("/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/"+'CMIP5_'+model+'-'+jobID+'-'+year+'lowhalf')
#		imageFolder 	= folder('images/'+'CMIP5/CMIP5_'+model+'-'+jobID+'lowhalf')
		workingDir = folder("/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/"+'CMIP5_'+model+'-'+jobID+'-'+year)
		imageFolder 	= folder('images/'+'CMIP5/CMIP5_'+model+'-'+jobID)
		
		shelvesAV.extend(
			testsuite_p2p(
				model 		= 'CMIP5_'+model,
				jobID 		= jobID,
				year  		= year,
				av 		= av,
				plottingSlices	= plottingSlices,
				workingDir 	= workingDir,
				imageFolder	= imageFolder,
				gridFile	= gridFile
		 	)
		)


	#BGCvsPhysics(shelvesAV, jobID, modelGrid )
	modelIntercomparisonAnnual(shelvesAV,folder('images/CMIP5/CMIP5_intercomparison'))
		
	if len(years)>1: InterAnnualPatterns(shelvesAV, jobID, years,modelGrid)	# plots interannual comparison and time series.
	
		
if __name__=="__main__":
	analysis()	
	
	
			

