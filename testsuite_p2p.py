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
from sys import argv,exit
from os.path import exists
from calendar import month_name

#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, AutoVivToYaml,YamlToDict, slicesDict,reducesShelves
from p2p import matchDataAndModel,makePlots,makeTargets, csvFromShelves, makePatternStatsPlots
#from 
from pftnames import MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names#,getmt

###	Potential problems?
###		Reliance on ORCA1 grid


def testsuite_p2p(	model='ERSEM',#'MEDUSA','ERSEM','NEMO'],
			year='1997',
			jobID='xhonc',
			av = {},
			plottingSlices = [],
			workingDir  = '',
			imageFolder = '',
			noPlots = False,
			gridFile=''
			):

	"""
	This analysis package performs the point to point analysis for a single model, for one year, for one job ID.
	
	Arguments:
	    model: Model name
	    
	    year: year (4 digit string), doesn't need to be exact. For climatologies, we typically use 2525 or something absurd.
	    
	    jobID: 5 letter jobID as used on monsoon. 
	    
	    av:
		The AutoVivification (av) is crucial. It controls the analysis. It locates the files. It identifies the field names in the netcdf.
	
		The av has a few very specific requiments in terms of structure and key words.
		Here is an example, of ERSEM chlorophyll:
			av['chl']['Data']['File'] 	= Observation_Filename_path_.netcdf	
			av['chl']['Data']['Vars'] 	= ['Chlorophylla',]		
			av['chl']['ERSEM']['File'] 	= Model_Filename_path_.netcdf
			av['chl']['ERSEM']['Vars'] 	= ['chl',]
			av['chl']['ERSEM']['grid']	= 'ORCA1'				
			av['chl']['depthlevels'] 	= ['',]
		where: 
			'File' is the file path
			'Vars' is the variable as it is call in the netcdf.
			'grid' is the model grid name. These grids are linked to a grid mesh file for calculating cell volume, surface area and land masks.
			'depthlevels' a list of depth levels. This is needed because some WOA files are huges and desktop computers may run the p2p analysis of that data.
				depthlevels options are ['', 'Surface','100m','200m','500m',]
				depthlevels = ['',] indicates look at no depth slicing. (Not recommended for big (>500,000 points) datasets! Don't say I didn't warn you.)

	    plottingSlices:
		plottingSlices is a list of regional, temporal, or statistical slices to be given to the analysis for plotting.
		ie:	plottingSlices = ['All', 				# plots everything,
					  'NorthAtlantic', 			# plots North Atlantic
					  'February',				# plots February
					  ('NorthAtlantic', 'February'),	# plots North Atlantic in February
					  ]
			plottingSlices can be made automatically with UKESMpthon.populateSlicesList()
			plottingSlices can also be added to the av:
				av[name]['plottingSlices'] = a list of slices
	    workingDir: 
	    	workingDir is a location for the working files that are produced during the analysis. 
	    	if no working directory is provided, the default is: ~/WorkingFiles/model-jobID-yyear
	    		
	    imageFolder: 
	    	imageFolder is a location for all the images  that are produced during the analysis. 
	    	if no working directory is provided, the default is: ~/images/model-jobID  	

	    noPlots:
	    	noPlots is a boolean value to turn off the production of images.
	    	This can streamline the analysis routine, if plots are not needed.
	    	
	Returns:
		shelvesAV:
		another AutoVivification with the following structure:
		shelvesAV[model][name][depthLevel][newSlice][xkey] = shelvePath
	
	testsuite_p2p is not the place for intercomparisons of models, years, or jobID. 
	This can be done after calling testsuite_p2p.py, and by using the 

	"""

	print "#############################"
	print "testsuite_p2p:  "
	print "models:        ",model
	print "year:          ",year
	print "jobID:         ",jobID
	print "av keys:	      ",sorted(av.keys())
	print "#############################"	
		
	if len( av.keys())==0:
		print "No autovivification nested dictionary given. - See testsuite_p2p documentation or a working example."
		exit(0)
	

	# Location of processing files
	if len( workingDir) == 0:
		workingDir = folder("WorkingFiles/"+model+'-'+jobID+'-'+year)
		print "No working directory provided, creating default:",workingDir
		
	# Location of image Output files
	if noPlots is False:
	    if len(imageFolder)==0:
		imageFolder 	= folder('images/'+model+'-'+jobID)
		print "No image directory provided, creating default:",imageFolder


	#####
	# Start analysis here:
	shelvesAV = []#AutoVivification()
			
	for name in sorted(av.keys()):
		#####
		# Start with some tests of the av.
		
		#####
		# Testing av for presence of model keyword
	    	print "testsuite_p2p: \t",model,jobID, year, name#, av[name][model]
	    	try: 
			if not isinstance(av[name][model],dict):
				print "testsuite_p2p: \tWARNING:", model ,' not in av', av[name].keys()
				continue
			if len(av[name][model].keys()) ==0:
				print "testsuite_p2p: \tWARNING:", model ,' not in av', av[name].keys()			
				continue
	    	except KeyError:
			print "testsuite_p2p: \tWARNING:\tNo ",name, 'in ',model
			continue
			
		#####
		# Testing av for presence of data keyword	
	    	try: 
			if not isinstance(av[name]['Data'],dict):
				print "testsuite_p2p: \tWARNING:", 'Data' ,' not in av', av[name].keys()
				continue
			if len(av[name]['Data'].keys()) ==0:
				print "testsuite_p2p: \tWARNING:", 'Data' ,' not in av', av[name].keys()	
				continue
	    	except KeyError:
			print "testsuite_p2p: \tWARNING:\tNo ",'Data', 'in ',model
			continue
			

		#####
		# Testing av for presence of data/obs files.
		try:
		    if not exists(av[name]['Data']['File']):
			print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name]['Data']['File']
			continue
		except:
			print "testsuite_p2p.py:\tWARNING:\tDict entry does not exist\tav[",name,"][",model,'][File]'
			continue			    	
		try:
		    if not exists(av[name][model]['File']):
			print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name][model]['File']
			continue
		except:
			print "testsuite_p2p.py:\tWARNING:\tDict entry does not exist:\tav[",name,"][",model,'][File] :',av[name][model]['File'] 
			continue			


	    	#####					
		# Testing av for presence of grid	
	    	grid = av[name][model]['grid']
	    	if grid in ['', [], {}, None]	: 
			print "testsuite_p2p.py:\tERROR:\tgrid not found:\tav[",name,"][",model,'][grid]: ',grid
			assert False

	    	#####					
		# Testing av for presence of depthLevels	
		if len(av[name]['depthLevels']) ==0:
				av[name]['depthLevels']=['',]
				print "testsuite_p2p: \tWARNING: no 'depthLevels' provided in av, using defaults: ['',]"

		#####
		# Made it though the initial tests. Time to start the analysis.
		print "\n\n\ntestsuite_p2p.py:\tINFO:\tMade it though initial tests. Running:",model,jobID, year, name, av[name]['depthLevels']
		for depthLevel in av[name]['depthLevels']:
			depthLevel = str(depthLevel)
			
			
			#####
			# matchDataAndModel:
			# Match observations and model. 
			# Does not produce and plots.
			b = matchDataAndModel(av[name]['Data']['File'], 
								av[name][model]['File'],
								name,
								DataVars  	= av[name]['Data']['Vars'],
								ModelVars 	= av[name][model]['Vars'],
								model 		= model,
								jobID		= jobID,
								year		= year,
								workingDir 	= folder(workingDir+name),
								depthLevel 	= depthLevel,
								grid		= grid,
								gridFile	= gridFile)
							
			#####
			# makePlots:
			# Make some plots of the point to point datasets.
			# MakePlot runs a series of analysis, comparing every pair in DataVars and ModelVars
			#	 under a range of different masks. For instance, only data from Antarctic Ocean, or only data from January.
			# The makePlot produces a shelve file in workingDir containing all results of the analysis.
			if len( plottingSlices) ==0:
				if len(av[name]['plottingSlices'])==0:
					plottingSlices = populateSlicesList()
					print "No plotting slices provided, using defaults",plottingSlices
				else:	plottingSlices = av[name]['plottingSlices']
					
					
			imageDir	= folder(imageFolder +'P2Pplots/'+year+'/'+name+depthLevel)	
			m = makePlots(	b.MatchedDataFile, 
					b.MatchedModelFile, 
					name, 
					newSlices 	= plottingSlices,
					jobID		= jobID,
					model 		= model,						
					depthLevel 	= depthLevel,
					year 		= year, 
					shelveDir 	= folder(workingDir+name+depthLevel),
					imageDir	= imageDir,
					compareCoords	=True,
					noPlots		= noPlots)

			#shelvesAV[model][name][depthLevel] = m.shelvesAV
			shelvesAV.extend(m.shelvesAV)
			
			#####
			# no plots doesn't produce any plots, but does produce the list of shelves which can be used in Taylor/Target/Pattern diagrams.			
			if noPlots: continue


			#csvFile = folder(workingDir+'/CSV')+'summary_file.csv'
			#print "attempting csvFromShelves:",m.shelves, csvFile
			#c = csvFromShelves.csvFromShelves(m.shelves, csvFile ,['check',])

			
			#####
			# makeTargets:
			# Make a target diagram of all matches for this particular dataset. 
			#filename = folder(imageFolder+'/Targets/'+year+'/AllSlices')+model+'-'+jobID+'_'+year+'_'+name+depthLevel+'.png'
			#t = makeTargets(	m.shelves, 
			#			filename,
			#			#name=name,
			#			legendKeys = ['newSlice','ykey',],
			#			debug=True)
			#			#imageDir='', 
			
			#####
			# Produce a set of pattern and a target plots for each of the groups here.
			groups = {'Oceans':[],'Months':[],'Seasons':[],'NorthHemisphereMonths':[],'SouthHemisphereMonths':[],'depthRanges':[]}
			for g in groups:
			    	groups[g] = reducesShelves(shelvesAV,  models =[model,],depthLevels = [depthLevel,], names = [name,], sliceslist =slicesDict[g])
				print g, groups[g]
				
				if len(groups[g])==0:continue 
				 
				#####
				# makeTargets:
				# Make a target diagram of the shelves of this group. 
			  	filename = folder(imageFolder+'/Targets/'+year+'/'+name+depthLevel+'/'+g)+model+'-'+jobID+'_'+year+'_'+name+depthLevel+'_'+g+'.png'
				makeTargets(	groups[g], 
						filename,
						legendKeys = ['newSlice',],					
						)
				#####
				# makePattern plots:
				# Make a pattern  diagram of all matches for this particular dataset. 
				xkeys=''
				for o in ['Oceans','Months','depthRanges']:
					if g.find(o)>=0:  xkeys=o
				if xkeys=='':
					print "Could no find x axis keys!",g,'in',['Oceans','Months']
					
			  	filenamebase = folder(imageFolder+'/Patterns/'+year+'/'+name+depthLevel+'/'+g)+'Months-'+model+'-'+jobID+'_'+year+'_'+name+depthLevel
				makePatternStatsPlots(	{name :groups[g],}, # {legend, shelves}
							name+' '+g,	#xkeysname
							slicesDict[xkeys],		#xkeysLabels=
							filenamebase,	# filename base	
							grid	= grid,												
							)
			#####
			# After finding all the shelves, we can plot them on the same axis.				
		  	filenamebase = folder(imageFolder+'/Patterns/'+year+'/'+name+depthLevel+'/ANSH')+'ANSH-Months-'+model+'-'+jobID+'_'+year+'_'+name+depthLevel
		  	
			makePatternStatsPlots(	{'North Hemisphere' :groups['NorthHemisphereMonths'],
						 'South Hemisphere' :groups['SouthHemisphereMonths'],
						 'Global' :	     groups['Months'], }, # {legend, shelves}
						name+' Months',	#xkeysname
						slicesDict['Months'],#xkeysLabels=
						filenamebase,	# filename base	
						grid	= grid,												
						)

		if noPlots: continue
		#####
		# And now by depth levels:
		
		groups = ['Oceans','Months','Seasons','depthRanges']	#'NorthHemisphereMonths':[],'SouthHemisphereMonths':[]}		
		for g in groups:
			if len(av[name]['depthLevels'])<=1: continue	
			outShelves = {}
			for dl in av[name]['depthLevels']:
				outShelves[dl] = reducesShelves(shelvesAV,  models =[model,],depthLevels = [dl,], names = [name,], sliceslist =slicesDict[g])	
		  	filenamebase = folder(imageFolder+'/Patterns/'+year+'/'+name+'AllDepths/')+'AllDepths-'+g+'-'+model+'-'+jobID+'_'+year+'_'+name
			makePatternStatsPlots(	outShelves, 
						name+' '+g,
						slicesDict[g],		
						filenamebase,
						grid	= grid,												
						)
									
	#assert False
			#####
			# EVERYTHING IS GREAT UNTIL THIS POINT, THEN IT ALL GOES TITS UP .
						

	return shelvesAV
																			
	#####			
	# Here are some fields for comparing fields in the same model
	surfacemetrics = ['chl', 'pCO2', 'nitrate',]
	dmsmetrics = ['dms_and','dms_ara','dms_hal','dms_sim']
	dmspmetrics = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']
	physics = ['salinity', 'temperature','irradiation',]
	WOA_bgc 	= ['silicate','nitrate','phosphate','oxygen']	
	WOA_phys 	= ['salinity','temperature',]
	PatternTypes = {'DMS_p':dmspmetrics,
			'DMS_e':dmsmetrics,
			'Maredat':MaredatTypes,
			#'WOA':WOATypes,
			'WOA BGC':WOA_bgc,
			'WOA physics':WOA_phys,					
			#'Nitrates':nitrates,
			#'Phosphates':phosphates,
			#'Silicates':silicates,
			'SurfaceMetrics':surfacemetrics,
			}
			
	Summary= {}	
	Summary['MaredatAll'] = []
	Summary['MaredatStandard'] = []		
	Summary['WOAAll'] = []
	Summary['WOAStandard'] = []	
	Summary['AllAll'] = []	
	Summary['AllStandard'] = []		
	Summary['SurfaceMetricsAll'] = []	
	Summary['SurfaceMetricsStandard'] = []	
				
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

	print "#####"
	print "Finished testsuite_p2p:",model, year,jobID
	print "working directory",workingDir
	print "image directory:", imageFolder
	
	return shelvesAV		
			


	
	
	
if __name__=="__main__":
	sys.exit(0)
	
	# Can use command line arguments to choose a model.
	models 		= []
	years 		= []
	jobIDs 		= []
	
	#####
	# Determine command line arguments
	if len(argv[1:]) <3:
		print "Not enough arguments supplied. Please supply (in any order): \n\t1. Model name \n\t2. Job ID \n\t3. Year. "
		exit(0)
		
	for a in argv[1:]:	
		try:	
			y = int(a)
			years.append(a)
			continue			
		except:pass
		if a == 'clim':
			years.append(a)
			continue
		
		if str(a).upper() in ['MEDUSA','ERSEM','NEMO']:
			models.append(str(a).upper())
			continue			
		if a[:4] in ['xhon','xjez']:
			jobIDs.append(a)
			#if 'ERSEM' not in models:models.append('ERSEM')
			continue

		if a[:4] in ['xkru',]:
			jobIDs.append(a)
			if 'MEDUSA' not in models:models.append('MEDUSA')
			continue
						
		print "Command line argument not understood:",a
	

	#####
	#Set Defaults:
	if not len(years): 	years = ['1998',]
	if not len(models): 	models = ['MEDUSA','ERSEM','NEMO']
	if not len(jobIDs):	jobIDs = ['xhonp',]	

	print "#############################"
	print "__main__ arguments: "
	print "models:        ",models
	print "year:          ",years
	print "jobID:         ",jobIDs
	print "#############################"

	
	for year in years:

		testsuite_p2p(models = models,	year=year,jobID=jobIDs[0] ) 
		if len(jobIDs)==1:continue
		for e in jobIDs[1:]:
			testsuite_p2p(models = ['ERSEM',],year=year,jobIDs=e ) 
	
	print 'The end.'
	



















	
