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

#Specific local code:
from UKESMpython import populateSlicesList, AutoVivification, folder# ,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict, slicesDict,reducesShelves

#from pftnames import MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names,getmt

from testsuite_p2p import testsuite_p2p
#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	
#	This code is the 

# 	from 
def analysis():
	model= 'NEMO'
	jobID = 'xhonc'
	year = '1997'		# arbitrary year chosen for 
			

	#####
	# Which analysis to run
	# Physics:
	doLight		= 0#True	# no file yet
	doSal		= True
	doTemp		= True
	doMLD		= True	

	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)
	depthLevels 	= ['Surface','100m','200m','500m',]#

	#####
	# Location of data files.
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	

	#####
	# Location of model files.	
	NEMOFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobID +'/'+year +'/'+jobID +'_'+year
	
	

	#####
	# AutoVivification is a form of nested dictionary.
	# We use AutoVivification here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	if doSal:
		av['salinity']['Data']['File'] 		= WOAFolder+'salinity_monthly_1deg.nc'	
		av['salinity']['NEMO']['File'] 		= NEMOFolder+'_NEMO.nc'	
		av['salinity']['Data']['Vars'] 		= ['s_an',]
		av['salinity']['NEMO']['Vars'] 		= ['vosaline',]
		av['salinity']['depthLevels'] 		= depthLevels	 
		av['salinity']['NEMO']['grid'] 		= 'ORCA1'
		
	if doTemp:
		av['temperature']['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
		av['temperature']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'	
		av['temperature']['Data']['Vars'] 	= ['t_an',]	
		av['temperature']['NEMO']['Vars'] 	= ['votemper',]
		av['temperature']['depthLevels'] 	= depthLevels	
		av['temperature']['NEMO']['grid'] 	= 'ORCA1'	
					   
	if doMLD:	
		av['mld']['Data']['File'] 		= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['mld']['NEMO']['File'] 		= NEMOFolder+'_NEMO.nc'	
		av['mld']['NEMO']['File'] 		= NEMOFolder+'_NEMO.nc'			
		av['mld']['Data']['Vars'] 		= ['mld','mask',]
		av['mld']['NEMO']['Vars'] 		= ['somxl010',]	
		av['mld']['depthLevels'] 		= ['',]
		av['mld']['NEMO']['grid'] 		= 'ORCA1'

				
		#av['mld_DR003']['Data']['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['mld_DR003']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'			
		#av['mld_DR003']['Data']['Vars'] 	= ['mld','mask',]
		#av['mld_DR003']['NEMO']['Vars'] 	= ['somxl010',]	
		#av['mld_DR003']['depthLevels'] 		= ['',]

		#av['mld_DReqDTm02']['Data']['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['mld_DReqDTm02']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'			
		#av['mld_DReqDTm02']['Data']['Vars'] 	= ['mld','mask',]
		#av['mld_DReqDTm02']['NEMO']['Vars'] 	= ['somxl010',]	
		#av['mld_DReqDTm02']['depthLevels'] 		= ['',]

	if doLight:	
		# Light file ? 
		#av['irradiation']['Data']['File'] 		= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['irradiation']['NEMO']['File'] 		= NEMOFolder+jobID+'_'+ year+"_MEDUSA_Light.nc"			
		av['irradiation']['Data']['Vars'] 		= []
		av['irradiation']['NEMO']['Vars'] 		= ['MED_QSR',]	
		av['irradiation']['depthLevels'] 			= ['',]
		av['irradiation']['NEMO']['grid'] 		= 'ORCA1'

	#####
	# Set which spatial and temporal limitations to plot.
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
		)
		
	workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/"+model+'-'+jobID+'-'+year)
	imageFolder 	= folder('images/'+model+'-'+jobID)
	shelvesAV = testsuite_p2p(
		model = model,
		jobID = jobID,
		year  = year,
		av = av,
		plottingSlices= plottingSlices,
		workingDir = workingDir,
		imageFolder= imageFolder)
	
	
			
if __name__=="__main__":
	analysis()	
	
							
