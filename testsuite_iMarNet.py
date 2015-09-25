#!/usr/bin/ipython 
#
# Copyright 2014, Plymouth Marine Laboratory
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
from sys import argv
from os.path import exists
from calendar import month_name

#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict, slicesDict
from p2p import matchDataAndModel,makePlots,makeTargets, csvFromShelves, makePatternStatsPlots
#from 
from pftnames import MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names,getmt

###	Potential problems?
###		Reliance on ORCA1 grid
from shelve import open as shOpen
import numpy as np
# Will not run this code in MO.
from ncdfView import ncdfView

def testsuite_iMarNet(	models=['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10'],
			year=1998,
			#ERSEMjobID='xhonc',
			plotallcuts = True,):

	
	#####
	# Can use command line arguments to choose a model.
	#if len(argv[1:]): models  = argv[1:]
	#else:	models = ['MEDUSA','ERSEM','NEMO']
    	
    	#####
    	# Which jobs to look at. 
	#ERSEMjobID = 'xhonp'
	jobIDs={}
	jobIDs['Diat-HadOCC'] 	= 'v1a'		
	jobIDs['ERSEM'] 	= 'xhonc'
	jobIDs['HadOCC'] 	= 'v1a'
	jobIDs['MEDUSA'] 	= 'v1'	
	jobIDs['PlankTOM6'] 	= 'xhyyp'
	jobIDs['PlankTOM10'] 	= 'xiokb'

	
	#####
	# Plot p2p for all regions/oceans, or just everything and "standard" cuts.
	#plotallcuts = True
	
	#####
	# Which Year to investigate for each model.
	# In an ideal world, they would all be the same, except that my current run is stuck in the queue.
	year = str(year)	
	years = {m:year for m in models}


	
	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	LesterFolder 	= "/data/euryale7/scratch/ledm/LestersReportData/"		
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
	#####
	# Location of model files.
	iMarNetFolder   = "/data/perseus2/scratch/ledm/iMarNet/OUTPUT_v2/"
	
	iMarNetFiles = {}
	iMarNetFiles['Diat-HadOCC']	= './Diat-HadOCC_v1a/diat-hadocc_'+year+'_2D_mon.nc3'
	iMarNetFiles['ERSEM']		= './ERSEM_v4/'+jobIDs['ERSEM']+'/ersem_'+jobIDs['ERSEM']+'_2D_'+year+'.nc'
	iMarNetFiles['HadOCC']		= './HadOCC_v1a/hadocc_'+year+'_2D_mon.nc3'
	iMarNetFiles['MEDUSA']		= './MEDUSA_v1/medusa_'+year+'_2D_mon.nc3'
	iMarNetFiles['PlankTOM6']	= './PlankTOM6_v1/TOM6_xhyyp_'+year+'_mon.nc'
	iMarNetFiles['PlankTOM10']	= './PlankTOM_v3/TOM10_xiokbo_'+year+'_mon.nc'
	

	
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	#ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobIDs['ERSEM']+'/'+years['ERSEM']+'/'+jobIDs['ERSEM']+'_'+years['ERSEM']
	#NEMOFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobIDs['NEMO'] +'/'+years['NEMO'] +'/'+jobIDs['NEMO'] +'_'+years['NEMO']
	

	regions = ['Surface',]#'100m','200m','500m',]
	
	#####
	# Which analysis to run
	#working - 2015-09-23:
	doCHL 		= 0#True
	doN		= 0#True
	doSi		= 0#True
	
	#in progress:
	doPCO2		= 0#True
	doIntPP		= True
	doO2		= 0#True		# not yet implemented.
	doDIC		= 0#True		# What database?
	doP		= 0#True			# Phosphate is not an iMarNet value in the monthly surface field.
	doFe		= 0#True 		# Iron is a challenge in IMarNet because we only keep surface values.

	#####
	# Set which spatial and temporal limitations to plot.
	#plotallcuts = True
	if plotallcuts:
		 plotDefaults		=True		 
		 plotMonths		=True
		 plotdepthRanges	=0
		 plotpercentiles	=0#True	
		 plotLatRegions		=0# True
		 plotQualityCuts	=0#True	
		 plotSeas		=0#True		 
		 plotOceans		=0#True
		 plotHemispheres	=0# True
		 plotSeasons		=0# True
		 plotOceanSeasons	=0# True		 		 
		 plotOceanMonths   	=0#True	
		 plotHemispheresMonths  =0#True			 
	else: 	
		 plotDefaults		=True		 	
		 plotMonths		=0#True
		 plotdepthRanges	=0#True	
		 plotpercentiles	=0#True	
		 plotLatRegions		=0#True
		 plotQualityCuts	=0#True
		 plotSeas		=0#True		 
		 plotOceans		=0#True	
		 plotHemispheres	=0		 
		 plotSeasons		=0# True
		 plotOceanSeasons	=0# True		 
		 plotOceanMonths   	=0	 	 	 
		 plotHemispheresMonths  =0			 

	if plotDefaults:	newSlices = ['All', 'Standard',]# Defaults
	else:			newSlices = []
	if plotMonths: 	 	newSlices.extend(slicesDict['Months'])
	if plotdepthRanges: 	newSlices.extend(slicesDict['depthRanges'])
	if plotpercentiles: 	newSlices.extend(slicesDict['percentiles'])
	if plotLatRegions:	newSlices.extend(slicesDict['latregions'])	
	if plotQualityCuts: 	newSlices.extend(slicesDict['QualityCuts'])		
	if plotSeas: 	 	newSlices.extend(slicesDict['Seas'])			
	if plotOceans: 	 	newSlices.extend(slicesDict['Oceans'])
	if plotHemispheres: 	newSlices.extend(slicesDict['Hemispheres'])	
	if plotSeasons: 	newSlices.extend(slicesDict['Seasons'])
	if plotOceanSeasons:	newSlices.extend(slicesDict['OceanSeasons'])		
	if plotOceanMonths: 	newSlices.extend(slicesDict['OceanMonths'])
	if plotHemispheresMonths: newSlices.extend(slicesDict['HemispheresMonths'])	


			
	#####
	# AutoVivification is a form of nested dictionary.
	# we use av here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	modelSkips = {	'Diat-HadOCC':	[],
			'ERSEM':	[],
			'HadOCC':	['si','pCO2'],
			'MEDUSA':	[],
			'PlankTOM6':	[],
			'PlankTOM10':	[],
		}
		
	av = AutoVivification()
	if doCHL:
		av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
		av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
		for m in models:
			if 'chl' in modelSkips[m]:continue
			av['chl'][m]['Vars'] 		= ['chl',]						
			av['chl'][m]['File']		= iMarNetFolder+iMarNetFiles[m]
			av['chl'][m]['grid']		= 'ORCA1'
		av['chl']['regions'] 			= regions
		
	if doN:
		av['nitrate']['Data']['File'] 		=  WOAFolder+'nitrate_monthly_1deg.nc'		
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 
		for m in models:
			if 'no3' in modelSkips[m]:continue		
			av['nitrate'][m]['Vars'] 	= ['no3',]						
			av['nitrate'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['nitrate'][m]['grid']	= 'ORCA1'
		av['nitrate']['regions'] 		= regions

	if doP:
		av['phosphate']['Data']['File'] 	=  WOAFolder+'phosphate_monthly_1deg.nc'		
		av['phosphate']['Data']['Vars'] 	= ['p_an',] 
		for m in models:
			if 'po4' in modelSkips[m]:continue		
			av['phosphate'][m]['Vars'] 	= ['po4',]						
			av['phosphate'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['phosphate'][m]['grid']	= 'ORCA1'
		av['nitrate']['regions'] 		= regions
				
	if doSi:
		av['silicate']['Data']['File'] 		=  WOAFolder+'silicate_monthly_1deg.nc'		
		av['silicate']['Data']['Vars'] 		= ['i_an',] 
		for m in models:
			if 'si' in modelSkips[m]:continue
			av['silicate'][m]['Vars'] 	= ['si',]
			av['silicate'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['silicate'][m]['grid']	= 'ORCA1'
		av['silicate']['regions'] 		= regions
		
	if doIntPP:
		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
		av['intpp']['Data']['Vars'] 	= ['PPint',]
		#av['intpp']['Data']['File'] 	=  MAREDATFolder+'PP100108.nc'
		#av['intpp']['Data']['Vars'] 	= ['PP',]
		for m in models:
			av['intpp'][m]['Vars'] 	= ['intpp',]
			av['intpp'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['intpp'][m]['grid']	= 'ORCA1'
		av['intpp']['regions'] 		= regions
		
	if doO2:
		av['oxygen']['Data']['File'] 	=  LesterFolder+"dissolved_oxygen_annual_1deg.nc"
		av['oxygen']['Data']['Vars'] 	= ['o_an',] 
		for m in models:
			av['oxygen'][m]['Vars'] 	= ['o2',]						
			av['oxygen'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['oxygen'][m]['grid']	= 'ORCA1'
		av['oxygen']['regions'] 		= ['Surface',]
					
	if doFe:
		av['iron']['Data']['File'] 	=  GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		av['iron']['Data']['Vars'] 	= ['Fe_D_CONC_BOTTLE',] 
		for m in models:
			av['iron'][m]['Vars'] 	= ['dfe',]						
			av['iron'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['iron'][m]['grid']	= 'ORCA1'
		av['iron']['regions'] 		= ['Surface',]
	
	if doPCO2:
		av['pCO2']['Data']['File'] 	=  TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'		
		av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		
		for m in models:
			if 'pCO2' in modelSkips[m]:continue		
			av['pCO2'][m]['Vars'] 	= ['spco2',]						
			av['pCO2'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
			av['pCO2'][m]['grid']	= 'ORCA1'
		av['pCO2']['regions'] 		=  ['',]
			
	
	#####
	# Run a quick test of the AV fields.	
	print "Performing test of AV dict:"
	for model in models:
		nctmp = ncdfView(iMarNetFolder+iMarNetFiles[m], Quiet=True)
		for name in sorted(av.keys()):	
		  for k in av[name][m]['Vars']: 
		  	if k in nctmp.variables.keys():continue
		  	print "WARNING:", k, 'not in ', model, av[name][m]['File']

	for name in sorted(av.keys()):
		try:	
			#print 'opening:', av[name]['Data']['File']
			nctmp = ncdfView(av[name]['Data']['File'], Quiet=True)
		except: 
			print 'Warning:', name ,"['Data']['File'] doesn't exist"
			continue

		for k in av[name]['Data']['Vars']: 
		  	if k in nctmp.variables.keys():continue
		  	print "WARNING:", k, 'not in ', name, 	av[name]['Data']['File']
		    
		    
	
	#####
	# Start analysis here:
	shelvesAV = AutoVivification()
	
	for model in models:
		for name in sorted(av.keys()):
		
		    for region in av[name]['regions']:
			#####
			# Do some checks to make sure that the files all exist:
			print model,name
			try: 
				if not isinstance(av[name][model],dict): continue
			except KeyError:
				print "No ",name, 'in ',model
				continue	
			region = str(region)
			
			
			#####
			# Location of image Output files
			imageFolder 	= folder('images/testsuite_iMarNet/'+model+'-'+jobIDs[model])
			workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/testsuite_iMarNet/"+model+'-'+jobIDs[model]+'-'+years[model])		

		
			try:
			    if not exists(av[name]['Data']['File']):
				print "testsuite_iMarNet.py:\tWARNING:\tFile does not exist", av[name]['Data']['File']
				continue
			except:
				print "testsuite_iMarNet.py:\tWARNING:\tFile does not exist\tav[",name,"][",model,'][File]'
				continue			    	
			try:
			    if not exists(av[name][model]['File']):
				print "testsuite_iMarNet.py:\tWARNING:\tFile does not exist", av[name][model]+'[File]'
				continue
			except:
				print "testsuite_iMarNet.py:\tWARNING:\tFile does not exist:\tav[",name,"][",model,'][File]'
				continue			
			print "\n\n\ntestsuite_iMarNet.py:\tINFO:\tRunning:",name
			
			
			#####
			# matchDataAndModel:
			# Match (real) Data and Model. 
			# Does not produce and plots.

			# Grid Testtry:	
			grid = av[name][model]['grid']
			#except: grid = 'ORCA1'
			if grid in ['', [], {}, None]	: 
				print "testsuite_p2p.py:\tERROR:\tgrid not found:\tav[",name,"][",model,'][grid]: ',grid
				assert False
			
			b = matchDataAndModel(av[name]['Data']['File'], 
								av[name][model]['File'],
								name,
								DataVars  	= av[name]['Data']['Vars'],
								ModelVars 	= av[name][model]['Vars'],
								model 		= model,
								jobID		= jobIDs[model],
								year		= years[model],
								workingDir 	= folder(workingDir+name+region),
								region 		= region,
								grid		= grid)
							
			#####
			# makePlots:
			# Make some plots of the point to point datasets.
			# MakePlot runs a series of analysis, comparing every pair in DataVars and ModelVars
			#	 under a range of different masks. For instance, only data from Antarctic Ocean, or only data from January.
			# The makePlot produces a shelve file in workingDir containing all results of the analysis.
			imageDir	= folder(imageFolder +'P2Pplots/'+years[model]+'/'+name+region)			
			m = makePlots(	b.MatchedDataFile, 
					b.MatchedModelFile, 
					name, 
					#model, 
					newSlices 	= newSlices,
					jobID 		= jobIDs[model],
					model 		= 'IMARNET_'+model,						
					region 		= region,
					year 		= years[model], 
					#plotallcuts	= plotallcuts, 
					shelveDir 	= folder(workingDir+name+region),
					imageDir	= imageDir,
					compareCoords	=True)

			shelvesAV[model][name.replace(region,'')][region] = m.shelvesAV
			csvFile = folder(workingDir+'/CSV')+'summary_file.csv'
			print "attempting csvFromShelves:",m.shelves, csvFile
			c = csvFromShelves.csvFromShelves(m.shelves, csvFile ,['check',])

			
			#####
			# makeTargets:
			# Make a target diagram of all matches for this particular dataset. 
			filename = folder(imageFolder+'/Targets/'+years[model]+'/AllSlices')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'.png'
			t = makeTargets(	m.shelves, 
						filename,
						#name=name,
						legendKeys = ['newSlice','ykey',],
						debug=True)
						#imageDir='', 
						
			
			#####
			# Ocean and month targets for this particular dataset.
			MonthShelves = []
			OceanShelves = []
			OceanMonthShelves = {o:[] for o in Ocean_names}
			OceanSeasonsShelves = {o:[] for o in Ocean_names}			
			HemisphereMonthShelves = {o:[] for o in Hemispheres}			


			for newSlice in m.shelvesAV.keys(): 
			   for xkey in m.shelvesAV[newSlice].keys():
				for ykey in m.shelvesAV[newSlice][xkey].keys():        	      
				    shelve = m.shelvesAV[newSlice][xkey][ykey]			  
				    if newSlice in month_name: 	MonthShelves.append(shelve)
				    if newSlice in Ocean_names:	OceanShelves.append(shelve)
				    
				    # This part is for Ocean + Month cuts.
				    if type(newSlice) in [type(['a','b',]),type(('a','b',))]:newSlice= newSlice[0]+newSlice[1]
				    if newSlice in OceanMonth_names:
				    	#print 'Prepping Targets',newSlice,xkey,ykey
				    	for mn in month_name:
				    	    t = newSlice.find(mn)
					    #print 'Prepping Targets',newSlice, mn,newSlice[:t],t
				    	    if t>0: 
				   	 	#print newSlice, mn,newSlice[:t],t
				    	    	OceanMonthShelves[newSlice[:t]].append(shelve)
				    if newSlice in OceanSeason_names:
				    	#print 'Prepping Targets',newSlice,xkey,ykey
				    	for mn in Seasons:
				    	    t = newSlice.find(mn)
					    #print 'Prepping Targets',newSlice, mn,newSlice[:t],t
				    	    if t>0: 
				   	 	#print newSlice, mn,newSlice[:t],t
				    	    	OceanSeasonsShelves[newSlice[:t]].append(shelve)				    	    	
				    if newSlice in HemispheresMonths:
				    	#print 'Prepping Targets',newSlice,xkey,ykey
				    	for mn in month_name:
				    	    t = newSlice.find(mn)
					    #print 'Prepping Targets',newSlice, mn,newSlice[:t],t
				    	    if t>0: 
				   	 	#print newSlice, mn,newSlice[:t],t
				    	    	HemisphereMonthShelves[newSlice[:t]].append(shelve)
			#continue	   
			#print OceanMonthShelves
			#assert False	    	
				    	
			if len(MonthShelves):	    
			  	filename = folder(imageFolder+'/Targets/'+years[model]+'/Months')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_Months.png'
				makeTargets(	MonthShelves, 
						filename,
						legendKeys = ['newSlice',],					
						)
			  	filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/'+name)+'Months-'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
				makePatternStatsPlots(	{name:MonthShelves,}, # {legend, shelves}
							'Months',	#xkeysname
							slicesDict['Months'],#xkeysLabels=
							filenamebase,	# filename base	
							grid	= grid,												
							)						
		
												
			if len(OceanShelves):	
				filename = folder(imageFolder+'/Targets/'+years[model]+'/Oceans')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_Oceans.png'
				makeTargets(	OceanShelves, 
						filename,
						legendKeys = ['newSlice',],					
						)
		  		filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/'+name)+'Oceans-'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
				makePatternStatsPlots(	{name:OceanShelves,}, # {legend, shelves}
							'Oceans',	#xkeysname
							slicesDict['Oceans'],#xkeysLabels=
							filenamebase,	# filename base	
							grid	= grid,
							)

			if len(OceanMonthShelves.keys()):										
				for o in OceanMonthShelves.keys():
				    if len(OceanMonthShelves[o]):
					filename = folder(imageFolder+'/Targets/'+years[model]+'/OceanMonths/'+o)+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_'+o+'Months.png'
					makeTargets(	OceanMonthShelves[o], 
							filename,
							legendKeys = ['newSlice',],					
							)
				counts = [len(sh) for i,sh in OceanMonthShelves.items()]
				if max(counts)>0:
				  	filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/'+name)+'OceanMonths-'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
					makePatternStatsPlots(	OceanMonthShelves, # {legend, shelves}
							'Ocean Months',	#xkeysname
							slicesDict['Months'],#xkeysLabels=
							filenamebase,	# filename base						
							grid	= grid,							
							)
			#assert False						
			for o in OceanSeasonsShelves.keys():
			    if len(OceanSeasonsShelves[o]):
				filename = folder(imageFolder+'/Targets/'+years[model]+'/OceanSeasons/'+o)+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_'+o+'Seasons.png'
				makeTargets(	OceanSeasonsShelves[o], 
						filename,
						legendKeys = ['newSlice',],					
						)
				counts = [len(sh) for i,sh in OceanSeasonsShelves.items()]
				if max(counts)>0:
				  	filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/'+name)+'OceanSeasons-'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
					makePatternStatsPlots(	OceanSeasonsShelves, # {legend, shelves}
							'Ocean Seasons',	#xkeysname
							slicesDict['Seasons'],	#xkeysLabels=
							filenamebase,	# filename base	
							grid	= grid,												
							)
													
			for o in HemisphereMonthShelves.keys():
			    if len(HemisphereMonthShelves[o]):
				filename = folder(imageFolder+'/Targets/'+years[model]+'/Hemispheres/'+o)+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_'+o+'Hemispheres.png'
				makeTargets(	HemisphereMonthShelves[o], 
						filename,
						legendKeys = ['newSlice',],					
						)
				counts = [len(sh) for i,sh in HemisphereMonthShelves.items()]
				if max(counts)>0:
				  	filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/'+name)+'HemisphereMonths-'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
					makePatternStatsPlots(	HemisphereMonthShelves, # {legend, shelves}
							'Hemisphere Months',	#xkeysname
							slicesDict['Months'],			#xkeysLabels=
							filenamebase,		# filename base	
							grid	= grid,												
							)
																			
			#print OceanMonthShelves
		#####				
		# Here are some fields for comparing fields in the same model
		Summary= {}		
		Summary['MaredatAll'] = []
		Summary['MaredatStandard'] = []			
		Summary['WOAAll'] = []
		Summary['WOAStandard'] = []	

		Summary['AllAll'] = []	
		Summary['AllStandard'] = []			
		Summary['SurfaceMetricsAll'] = []	
		Summary['SurfaceMetricsStandard'] = []			
		surfacemetrics = ['chl', 'pCO2', 'nitrate',]
							
		dmsmetrics = ['dms_and','dms_ara','dms_hal','dms_sim']
		dmspmetrics = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']	
		nitrates   = ['nitrate'+ s   for s in regions]	
		phosphates = ['phosphate'+ s for s in regions]	
		silicates  = ['silicate'+ s  for s in regions]	
						
		Summary['DMSAll'] = []
		Summary['DMSStandard']=[]	
		Summary['DMS_p_All'] = []
		Summary['DMS_p_Standard']=[]

		PatternTypes = ['DMS_p','DMS_e','Maredat','WOA','Nitrates','Phosphates','Silicates']
		MonthsPatterns = {p:{} for p in PatternTypes}
		OceansPatterns = {p:{} for p in PatternTypes}	
		
		#'DMS_p':{},'DMS_e':{},'Maredat':{},'WOA':{},''}
		#OceansPatterns = {'DMS_p':{},'DMS_e':{},'Maredat':{},'WOA':{},}		

		for name_init in shelvesAV[model].keys():
		    for r in av[name_init]['regions']:
		        name =name_init+r
			if name in MaredatTypes: 
				OceansPatterns['Maredat'][name] = []
			  	MonthsPatterns['Maredat'][name] = []
			if name in WOATypes: 
				OceansPatterns['WOA'][name] = []
			  	MonthsPatterns['WOA'][name] = []
			if name in dmsmetrics: 
				OceansPatterns['DMS_e'][name] = []
			  	MonthsPatterns['DMS_e'][name] = []
			if name in dmspmetrics: 
				OceansPatterns['DMS_p'][name] = []
			  	MonthsPatterns['DMS_p'][name] = []

			if name in nitrates: 
				OceansPatterns['Nitrates'][name] = []
			  	MonthsPatterns['Nitrates'][name] = []
			if name in phosphates: 
				OceansPatterns['Phosphates'][name] = []
			  	MonthsPatterns['Phosphates'][name] = []
			if name in silicates: 
				OceansPatterns['Silicates'][name] = []
			  	MonthsPatterns['Silicates'][name] = []				  	
			

		for name in shelvesAV[model].keys():
		  for region in shelvesAV[model][name].keys():
		    for newSlice in shelvesAV[model][name][region].keys(): 
		      for xkey in shelvesAV[model][name][region][newSlice].keys():
			for ykey in shelvesAV[model][name][region][newSlice][xkey].keys():        	      
			  	shelve = shelvesAV[model][name][region][newSlice][xkey][ykey]
			       	namer =name+region			  	
	       		  	if newSlice == 'All':		Summary['AllAll'].append(shelve)
	       		  	if newSlice == 'Standard':	Summary['AllStandard'].append(shelve)

				if namer in surfacemetrics:
				  	if newSlice == 'All':		Summary['SurfaceMetricsAll'].append(shelve)
				  	if newSlice == 'Standard':	Summary['SurfaceMetricsStandard'].append(shelve)
				  		       		  	
				if namer in MaredatTypes:
	        		  	if newSlice == 'All':		Summary['MaredatAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['MaredatStandard'].append(shelve)
					if newSlice in Ocean_names:	OceansPatterns['Maredat'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['Maredat'][namer].append(shelve)
											        		  	
	        		if namer in WOATypes:
	        		  	if newSlice == 'All':		Summary['WOAAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['WOAStandard'].append(shelve)	
					if newSlice in Ocean_names:	OceansPatterns['WOA'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['WOA'][namer].append(shelve)
	        		
	        		if namer in dmsmetrics:
	        		  	if newSlice == 'All':		Summary['DMSAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['DMSStandard'].append(shelve)
					if newSlice in Ocean_names:	OceansPatterns['DMS_e'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['DMS_e'][namer].append(shelve)
	        		  		        				
	        		if namer in dmspmetrics:
	        		  	if newSlice == 'All':  		Summary[ 'DMS_p_All'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary[ 'DMS_p_Standard'].append(shelve)
					if newSlice in Ocean_names:	OceansPatterns['DMS_p'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['DMS_p'][namer].append(shelve)
	        		for woa in ['silicate','nitrate','phosphate','salinity','temperature','iron',]:
	        		   for ns in ['All', 'Standard']:
	        		   	if ns == newSlice and woa == namer.lower():
	        		   		try: 	Summary[woa+ns].append(shelve)
	        		   		except:	Summary[woa+ns]= [shelve,]

				if namer in nitrates:
					print 'loop Found nitrates:',namer, name, region, newSlice, shelve				
					if newSlice in Ocean_names:	OceansPatterns['Nitrates'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['Nitrates'][namer].append(shelve)	        
																	
				if namer in phosphates:
					if newSlice in Ocean_names:	OceansPatterns['Phosphates'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['Phosphates'][namer].append(shelve)	
					
				if namer in silicates:
					if newSlice in Ocean_names:	OceansPatterns['Silicates'][namer].append(shelve)
					if newSlice in months:		MonthsPatterns['Silicates'][namer].append(shelve)	
					  	#Patterns['DMS_p'][name].append(shelve)
					  	#if newSlice not in patternShelves: patternShelves.append(newSlice)
	        print 'OceansPatterns:', OceansPatterns
	        print 'MonthsPatterns:', MonthsPatterns	        
		for k in Summary.keys():
			filename = folder(imageFolder+'/Targets/'+years[model]+'/Summary')+model+'_'+years[model]+'_'+k+'.png'
			
	  		makeTargets(Summary[k], 
					filename,
					legendKeys = ['name',],
					debug=True)
					
		for k in OceansPatterns.keys():
			print 'OceansPatterns:',k			
			if len(OceansPatterns[k]) ==0: continue				
			filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/OceansPatterns/')+k+'_'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+region
			print 'OceansPatterns:',k, OceansPatterns[k],filenamebase
			
			makePatternStatsPlots(	OceansPatterns[k], # {legend, shelves}
					k,	#xkeysname
					Ocean_names,			#xkeysLabels=
					filenamebase,		# filename base	
					grid	= grid,											
					)
		for k in MonthsPatterns.keys():	
			print 'MonthsPatterns:',k
			if len(MonthsPatterns[k]) ==0: continue										
			filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/MonthsPatterns/')+k+'_'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+region
			print 'MonthsPatterns:',k, MonthsPatterns[k],filenamebase
			makePatternStatsPlots(	
					MonthsPatterns[k], # {legend, shelves}
					k,	#xkeysname
					months,			#xkeysLabels=
					filenamebase,		# filename base						
					grid	= grid,
					)	
																						
		#if model=='ERSEM':
		AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+jobIDs[model]+'.yaml')
		#else:
		#	AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+'.yaml')				

			

	#####
	SouthHemispheresMonths = [(h,m) for h in ['SouthHemisphere',] for m in months] 	
	NorthHemispheresMonths = [(h,m) for h in ['NorthHemisphere',] for m in months] 	
	
	# Model	Intercomparison Pattern plots
	xaxes = {'Oceans': Ocean_names, 
		 'Months':months, 
		 'HemispheresMonths':HemispheresMonths,
		 'SouthHemispheresMonths':SouthHemispheresMonths,
		 'NorthHemispheresMonths':NorthHemispheresMonths,
		 }

	ModelIntercomparisons = {pt:{} for  pt in xaxes.keys()}
	for pt in xaxes.keys():
	    for name in av.keys():
	      for r in av[name]['regions']:
	        namer = name+r
		ModelIntercomparisons[pt][namer] = {}
		for model in av[name].keys():
			if model.lower() in ['regions','data',]:continue
			ModelIntercomparisons[pt][namer][model] = []

	for model in shelvesAV.keys():
	 for name in shelvesAV[model].keys():
	  for region in shelvesAV[model][name].keys():
	    for newSlice in shelvesAV[model][name][region].keys(): 
	      for xkey in shelvesAV[model][name][region][newSlice].keys():
		for ykey in shelvesAV[model][name][region][newSlice][xkey].keys():        	      
		  	shelve = shelvesAV[model][name][region][newSlice][xkey][ykey]
		       	namer =name+region	
		       	if newSlice in NorthHemispheresMonths:
		       		print 'NorthHemispheresMonths:', model,name,region,newSlice,xkey,ykey,shelve, [pt,namer,model]
		       	for pt,sliceList in xaxes.items():
		       	    if newSlice in sliceList:	ModelIntercomparisons[pt][namer][model].append(shelve) 
	       		print 'ModelIntercomparisons loop:', model,name,region,newSlice,xkey,ykey,shelve   


	for pt in ModelIntercomparisons.keys():
	    for namer in ModelIntercomparisons[pt].keys():
	    	print 'ModelIntercomparisons:',pt, namer, ModelIntercomparisons[pt][namer].keys()
	        count = sum([len(shelves) for shelves in ModelIntercomparisons[pt][namer].values()])
	    	if count ==0:continue
		filenamebase = folder('images/testsuite_iMarNet/ModelIntercomparison/Patterns/'+namer)+namer+'_'+pt+'_'+str(year)
		print 'ModelIntercomparisons:',k, ModelIntercomparisons[pt][namer],xaxes[pt],filenamebase
		
		if pt.find('Oceans')>-1:xaxisLabels = Ocean_names
		else:			xaxisLabels = months
		makePatternStatsPlots(	
					ModelIntercomparisons[pt][namer], # {legend, shelves}
					namer,				# xkeysname
					xaxisLabels,			# xkeysLabels
					filenamebase,			# filename base
					grid	= 'ORCA1',
					)
	print "Working dir:",workingDir
			
			

def makeCSVfile():
	# this was done for the workshop.
	shelves = [	'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/Diat-HadOCC-v1a-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
			'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/ERSEM-xhonp-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
#			'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/ERSEM-xhont-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
			'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/HadOCC-v1a-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
			'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/MEDUSA-v1-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
			'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/PlankTOM10-xiokb-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
			'/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/PlankTOM6-xhyyp-1998/chlSurface/chl_All_chlvsChlorophylla.shelve',
			]
	mchl = {}
	dchl = {}
	lat = {}
	lon = {}
	keys = {}
	xtime = {}
	ytime = {}
			
	lengths = []
	models = []
	for s in shelves:
		sh = shOpen(s)
		m = s.replace('/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/','').replace('/chlSurface/chl_All_chlvsChlorophylla.shelve','')
		models.append(m)
		keys[m]  = sh['xkey']
		mchl[m] = sh['datax']
		dchl[m] = sh['datay']
		xtime[m] = sh['x_time']
		ytime[m] = sh['y_time']				
		lat[m] = sh['x_lat']
		lon[m] = sh['x_lon']		
		if len(sh['datax']) not in lengths: lengths.append(len(sh['datax']))
	out = 'Latitude, Longitude, month, MAREDAT CHL,'
	units = 'degrees N, degrees E, [0=Jan 11=Dec.], mgChl/m3,'
	models = sorted(models)
	out  += ','.join([m for m in models]) + '\n'
	units += ','.join(['mgChl/m3' for m in models]) + '\n'
	
	out+=units
	if len(lengths) != 1:assert False
	
	for i in range(lengths[0]):
		m = 'Diat-HadOCC-v1a-1998'
		arrs = [lat[m][i], lon[m][i],xtime[m][i], dchl[m][i],]
		arrs.extend([mchl[m][i] for m in models])
		line = ','.join([str(round(a,3)) for a in arrs])
		out+=line+'\n'
	
	f = open('iMarNetChl.csv','w')
	f.write(out)
	f.close()
	
if __name__=="__main__":
	#makeCSVfile()
	
	#assert False
	# Can use command line arguments to choose a model.
	models 		= []
	years 		= []
	#ERSEMjobIDs 	= []
	
	#####
	# Determine command line arguments
	for a in argv[1:]:	
		try:	
			y = int(a)
			years.append(a)
			continue			
		except:pass
		
		if a in ['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10']:
			models.append(a)
			continue			
		if a[:4] in ['xhon','xjez']:
			#ERSEMjobIDs.append(a)
			if 'ERSEM' not in models:models.append('ERSEM')
			continue
			
		print "Command line argument not understood:",a
	

	#####
	#Set Defaults:
	if not len(years): 	years = ['1998',]
	if not len(models): 	models =['PlankTOM10','PlankTOM6','Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA',]
	#if not len(ERSEMjobIDs):ERSEMjobIDs = ['xhonp',]	

	print "#############################"
	print "__main__ arguments: "
	print "models:        ",models
	print "year:          ",years
	#print "ERSEM jobID:   ",ERSEMjobIDs
	print "#############################"

	
	for year in years:
		testsuite_iMarNet(models = models,	year=year, ) 
	
	print 'The end.'
	



















	
