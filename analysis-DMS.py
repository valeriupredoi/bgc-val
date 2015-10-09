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
from UKESMpython import populateSlicesList, AutoVivification,folder # ,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict, slicesDict,reducesShelves
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
def analysisDMS():
	model= 'MEDUSA'
	jobID = 'xkrum'
	year = '2526'		# arbitrary year chosen for 
	
	# DMS:
	doDMS_clim	= True
	doDMS_pixels	= True
	doDMS_pixels2	= 0#True
	

	

	#####
	# AutoVivification is a form of nested dictionary.
	# We use AutoVivification here to determine which files to analyse and which fields in those files.
	# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	
	if doDMS_clim:
		dmsd= {'dms_and':'anderson','dms_ara':'aranamit','dms_hal':'halloran','dms_sim':'simodach'} 
		for dms in ['dms_and','dms_ara','dms_hal','dms_sim',]:
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'		
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['Data']['Vars'] 		= ['lanaetal',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]
			av[dms]['depthLevels'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] 		= 'Flat1deg'

	if doDMS_pixels:
		dmsd= {'dms_p_and':'anderson','dms_p_ara':'aranamit','dms_p_hal':'halloran','dms_p_sim':'simodach'} 
		for dms in ['dms_p_sim','dms_p_and','dms_p_ara','dms_p_hal',]:
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/DMS_Lana2011nc/DMSpixels.nc'		
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['Data']['Vars'] 		= ['DMS',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]			
			av[dms]['depthLevels'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] 		= 'Flat1deg'

	if doDMS_pixels2:
		#dmsd= {'dms_p_and':'anderson','dms_p_ara':'aranamit','dms_p_hal':'halloran','dms_p_sim':'simodach'} 
		dmsd= { 'dms_p_and1':'anderson1','dms_p_ara1':'aranamit1','dms_p_hal1':'halloran1','dms_p_sim1':'simodach1',
			'dms_p_and2':'anderson2','dms_p_ara2':'aranamit2','dms_p_hal2':'halloran2','dms_p_sim2':'simodach2',} 
		for d in ['dms_p_ara','dms_p_hal', 'dms_p_sim','dms_p_and',]:#
		  for n in ['1','2']:
		  	dms=d+n
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/DMS_Lana2011nc/DMSpixels.nc'		
			#av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum_v2.nc'				
			av[dms]['Data']['Vars'] 		= ['DMS',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]			
			#av[dms]['MEDUSA']['Vars'] 		= dmsd[dms]
			av[dms]['depthLevels'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] = 'Flat1deg'
	
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
	analysisDMS()			 
