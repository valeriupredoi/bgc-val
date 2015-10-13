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
from UKESMpython import populateSlicesList, AutoVivification,folder,reducesShelves # ,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict, slicesDict,reducesShelves
from pftnames import MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names#,getmt
from p2p import makeTargets, makePatternStatsPlots
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
	
	grid = 'Flat1deg'
	

	#####
	# AutoVivification is a form of nested dictionary.
	# We use AutoVivification here to determine which files to analyse and which fields in those files.
	# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	shelvesAV = []
	
	if doDMS_clim:
		dmsd= {'dms_and':'anderson','dms_ara':'aranamit','dms_hal':'halloran','dms_sim':'simodach'} 
		for dms in ['dms_and','dms_ara','dms_hal','dms_sim',]:
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'		
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['Data']['Vars'] 		= ['lanaetal',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]
			av[dms]['depthLevels'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] 		= grid

	if doDMS_pixels:
		dmsd= {'dms_p_and':'anderson','dms_p_ara':'aranamit','dms_p_hal':'halloran','dms_p_sim':'simodach'} 
		for dms in ['dms_p_sim','dms_p_and','dms_p_ara','dms_p_hal',]:
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/DMS_Lana2011nc/DMSpixels.nc'		
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['Data']['Vars'] 		= ['DMS',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]			
			av[dms]['depthLevels'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] 		= grid

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
			av[dms]['MEDUSA']['grid'] 		= grid
	
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
	shelvesAV.extend(testsuite_p2p(
		model = model,
		jobID = jobID,
		year  = year,
		av = av,
		plottingSlices= plottingSlices,
		workingDir = workingDir,
		imageFolder= imageFolder))
	
	# Here are some fields for comparing fields in the same model
	dmsmetrics = ['dms_and','dms_ara','dms_hal','dms_sim']
	dmspmetrics = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']
	PatternTypes = {'DMS_p':dmspmetrics,
			'DMS_e':dmsmetrics,
			}
			
	Summary= {}	
	Summary['AllAll'] = []	
	Summary['AllStandard'] = []		
				
	Summary['DMSAll'] = []
	Summary['DMSStandard']=[]	
	Summary['DMS_p_All'] = []
	Summary['DMS_p_Standard']=[]
				
	MonthsPatterns = {p:{} for p in PatternTypes.keys()}
	OceansPatterns = {p:{} for p in PatternTypes.keys()}
	SHMonthsPatterns = {p:{} for p in PatternTypes.keys()}
	NHMonthsPatterns = {p:{} for p in PatternTypes.keys()}
	#####
	SouthHemispheresMonths = [(h,m) for h in ['SouthHemisphere',] for m in months] 
	NorthHemispheresMonths = [(h,m) for h in ['NorthHemisphere',] for m in months]
			  	
		
	Summary['AllAll'] = 	 reducesShelves(shelvesAV, models =[model,], sliceslist = ['All',])
	Summary['AllStandard'] = reducesShelves(shelvesAV, models =[model,], sliceslist = ['Standard',])
	
	for TypeListName,names in PatternTypes.items():
	    for name in names:
	    	print "loading shelve meta data:",TypeListName,':',name
		SHMonthsPatterns[TypeListName][name] = reducesShelves(shelvesAV,  models =[model,],depthLevels = ['Surface',], names = [name,], sliceslist = SouthHemispheresMonths)
		NHMonthsPatterns[TypeListName][name] = reducesShelves(shelvesAV,  models =[model,],depthLevels = ['Surface',], names = [name,], sliceslist = NorthHemispheresMonths)
		
		OceansPatterns[TypeListName][name] = reducesShelves(shelvesAV,  models =[model,],depthLevels = ['Surface',], names = [name,], sliceslist = Ocean_names)
		MonthsPatterns[TypeListName][name] = reducesShelves(shelvesAV,  models =[model,],depthLevels = ['Surface',], names = [name,], sliceslist = months)
				

	depthLevel = ''
	print 'OceansPatterns:', OceansPatterns
	print 'MonthsPatterns:', MonthsPatterns
	for k in Summary.keys():
		filename = folder(imageFolder+'/Targets/'+year+'/Summary')+model+'_'+year+'_'+k+'.png'
		
	  	makeTargets(Summary[k], 
			filename,
			legendKeys = ['name',],
			debug=True)
			
	for k in OceansPatterns.keys():
		if len(OceansPatterns[k]) ==0: continue			
		filenamebase = folder(imageFolder+'/Patterns/'+year+'/OceansPatterns/')+k+'_'+model+'-'+jobID+'_'+year+'_'+depthLevel
		print 'OceansPatterns:',k, OceansPatterns[k],filenamebase
		
		makePatternStatsPlots(	OceansPatterns[k], # {legend, shelves}
			k,	#xkeysname
			Ocean_names,		#xkeysLabels=
			filenamebase,	# filename base	
			grid	= grid,							
			)
	for k in MonthsPatterns.keys():	
		if len(MonthsPatterns[k]) ==0: continue							
		filenamebase = folder(imageFolder+'/Patterns/'+year+'/MonthsPatterns/')+k+'_'+model+'-'+jobID+'_'+year+'_'+depthLevel
		print 'MonthsPatterns:',k, MonthsPatterns[k],filenamebase
		makePatternStatsPlots(	
			MonthsPatterns[k], # {legend, shelves}
			k,	#xkeysname
			months,		#xkeysLabels=
			filenamebase,	# filename base				
			grid	= grid,
			)	

	for k in SHMonthsPatterns.keys():	
		if len(SHMonthsPatterns[k]) ==0: continue							
		filenamebase = folder(imageFolder+'/Patterns/'+year+'/SHMonthsPatterns/')+k+'_'+model+'-'+jobID+'_'+year+'_'+depthLevel
		print 'SHMonthsPatterns:',k, SHMonthsPatterns[k],filenamebase
		makePatternStatsPlots(	
			SHMonthsPatterns[k], # {legend, shelves}
			'South hemisphere '+k,	#xkeysname
			months,		#xkeysLabels=
			filenamebase,	# filename base				
			grid	= grid,
			)
	for k in NHMonthsPatterns.keys():	
		if len(NHMonthsPatterns[k]) ==0: continue							
		filenamebase = folder(imageFolder+'/Patterns/'+year+'/NHMonthsPatterns/')+k+'_'+model+'-'+jobID+'_'+year+'_'+depthLevel
		print 'NHMonthsPatterns:',k, NHMonthsPatterns[k],filenamebase
		makePatternStatsPlots(	
			NHMonthsPatterns[k], # {legend, shelves}
			'North hemisphere '+k,	#xkeysname
			months,		#xkeysLabels=
			filenamebase,	# filename base				
			grid	= grid,
			)			

	
	
			
if __name__=="__main__":
	analysisDMS()			 
