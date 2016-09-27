#!/usr/bin/ipython 
#
# Copyright 2014, Plymouth Marine Laboratory
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
from sys import argv
import os
from calendar import month_name
import getpass

#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict
from p2p import matchDataAndModel,makePlots,makeTargets, csvFromShelves

from bgcvaltools.pftnames import MaredatTypes,WOATypes,Ocean_names,getmt

###	Potential problems?
###		Reliance on ORCA1 grid


#####
# This is a copy of the testsuite_jamin.py script, but tailored to take an AutoAsses dictionary. 

def testsuite_AutoAssess(run,):

	#####
	# Do some tests to determine which model we're lookinat at here, MEDUSA or NEMO.
	models = [run['ocean_model'],]
	if 'UM' in models:
		print "testsuite_AutoAssess:\tWARNING:\tThis code is not made for the UM model"
		return
	
	if 'NEMO' not in models and 'MEDUSA' not in models:
		print "testsuite_AutoAssess:\tWARNING:\tThis code is not made for the ",models," model, only for", ['NEMO','MEDUSA']
		return		
	
	
	#####
	# For each model, load the appropriate year and job ID.
	for m in models:
		jobIDs  = {m:run['runid']}
		try:	years 	= {m:str(run['end_year'].year)}
		except:	years 	= {m:str(run['end_year'])}	# This needs some work, at the moment it only runs the final year.
		


	
	#####
	# Name outgoing summary file, and ensure that the outgoing folder exists.
	csvFile 	= run['summary_file']	
	csvFileFold 	= folder(os.path.dirname(os.path.realpath((csvFile))))
	
	
	#####
	# Grid used to produce this run. Options are ORCA1 or ORCA025
	grid = run['grid']
	gridAvailable = ['ORCA1', 'ORCA025']
	if grid not in gridAvailable:
		print "testsuite_AutoAssess:\tWARNING:\tThis grid is not appropriate, options are ", gridAvailable
		assert False
		
	#####
	# Which regions to look at. More regions are
	regions = ['Surface',]
			
	#####
	# Plot p2p for all regions/oceans, or just everything and "standard" cuts.
	plotallcuts = False
	


	#####
	# Location of data files is now given in the run dictionary.
	esmvalFolder 	= folder(run['data_root'])		
	MAREDATFolder 	= folder(esmvalFolder+"/MAREDAT/")
	WOAFolder 	= folder(esmvalFolder+"/WOA/")
	GEOTRACESFolder = folder(esmvalFolder+"/GEOTRACES/GEOTRACES_PostProccessed/")
	TakahashiFolder = folder(esmvalFolder+"/Takahashi2009_pCO2/")
	iFERMERDFolder  = folder(esmvalFolder+"/IFREMER-MLD/")
	
	#####
	# Location of model files.	
	MEDUSAFolder	= run['ss_annual']
	#NEMOFolder	= folder(esmvalFolder+"NEMO/"+ jobIDs['NEMO'] +'/'+years['NEMO'] +'/'+jobIDs['NEMO'] +'_'+years['NEMO'])
	NEMOFolder	= run['ss_annual']
	
	# Directory for output files:
	postprocDir 	= folder(run['LOCALDATA']+"/ukesm_postProcessed_"+getpass.getuser())
	imageDir	= folder('images')


				
	#####
	# regions for large datasets. (WOA)
	# Moved into function call
	# regions = ['Surface',]#'100m','200m','500m',]
	
	#####
	# Which analyses to run
	if 'NEMO' in models:
		doSalTemp	= True
		doMLD		= True	
	else:
		doSalTemp	= False
		doMLD		= False
		
	if 'MEDUSA' in models:
		doCHL 		= True
		doMAREDAT 	= True
		doN		= True
		doPSF		= True	
		doPCO2		= True
	else:
		doCHL 		= False
		doMAREDAT 	= False
		doN		= False
		doPSF		= False	
		doPCO2		= False		
	
	
	#####
	# AutoVivification is a form of nested dictionary.
	# we use av here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	if doCHL:
		av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
		av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
		av['chl']['MEDUSA']['Vars'] 		= ['CHL',]	
		av['chl']['regions'] 			= ['',]
		
	if doMAREDAT:
		av['diatoms']['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
		av['diatoms']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
		av['diatoms']['MEDUSA']['Vars'] 	= ['PHD',]	
		av['diatoms']['regions'] 		= ['',]	
	

		
		av['microzoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
		av['microzoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['microzoo']['MEDUSA']['Vars'] 	= ['ZMI',]	
		av['microzoo']['regions'] 		= ['',]	
	
		av['mesozoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
		av['mesozoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['mesozoo']['MEDUSA']['Vars'] 	= ['ZME',]	
		av['mesozoo']['regions'] 		= ['',]

	if doN:
		av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
		av['nitrate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
		av['nitrate']['MEDUSA']['Vars'] 	= ['DIN',]									
		av['nitrate']['regions'] 		= regions
	if doPSF:
		av['silicate']['Data']['File'] 		= WOAFolder+'silicate_monthly_1deg.nc'	
		av['silicate']['Data']['Vars'] 		= ['i_an',] 		#l+'_mn',
		av['silicate']['MEDUSA']['Vars'] 	= ['SIL',]									
		av['silicate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
		av['silicate']['regions'] 		= regions			
					
		av['iron']['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		av['iron']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['iron']['Data']['Vars'] 		= ['Fe_D_CONC_BOTTLE',]
		av['iron']['MEDUSA']['Vars'] 		= ['FER',]	
		av['iron']['regions'] 			= ['',]
		
	if doSalTemp:
		av['salinity']['Data']['File'] 		= WOAFolder+'salinity_monthly_1deg.nc'	
		av['salinity']['NEMO']['File'] 		= NEMOFolder+jobIDs['NEMO'] +'_'+years['NEMO']+'_NEMO.nc'
		av['salinity']['Data']['Vars'] 		= ['s_an',]
		av['salinity']['NEMO']['Vars'] 		= ['vosaline',]
		av['salinity']['regions'] 		= regions	 

		av['temperature']['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
		av['temperature']['NEMO']['File'] 	= NEMOFolder+jobIDs['NEMO'] +'_'+years['NEMO']+'_NEMO.nc'	
		av['temperature']['Data']['Vars'] 	= ['t_an',]	
		av['temperature']['NEMO']['Vars'] 	= ['votemper',]
		av['temperature']['regions'] 		= regions	
						
				   
	if doMLD:	
		av['mld']['Data']['File'] 		= iFERMERDFolder+ "mld_DT02_c1m_reg2.0.nc"
		av['mld']['NEMO']['File'] 		= NEMOFolder+jobIDs['NEMO'] +'_'+years['NEMO']+'_NEMO.nc'		
		av['mld']['Data']['Vars'] 		= ['mld','mask',]
		av['mld']['NEMO']['Vars'] 		= ['somxl010',]	
		av['mld']['regions'] 			= ['',]


		
	if doPCO2:
		av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
		av['pCO2']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
		av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
		av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]	
		av['pCO2']['regions'] 		= ['',]
	

	
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
			imageFolder 	= folder(imageDir+model+'-'+jobIDs[model])
			postprocFolder  = folder(postprocDir + model+'-'+jobIDs[model]+'-'+years[model])

		
			try:
			    if not os.path.exists(av[name]['Data']['File']):
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name]['Data']['File']
				continue
			except:
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist\tav[",name,"][",model,'][File]'
				continue			    	
			try:
			    if not os.path.exists(av[name][model]['File']):
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name][model]+'[File]'
				continue
			except:
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist:\tav[",name,"][",model,'][File]'
				continue			
			print "\n\n\ntestsuite_p2p.py:\tINFO:\tRunning:",name
			
			
			#####
			# matchDataAndModel:
			# Match (real) Data and Model. 
			# Does not produce and plots.
			b = matchDataAndModel(av[name]['Data']['File'], 
								av[name][model]['File'],
								name,
								DataVars  	= av[name]['Data']['Vars'],
								ModelVars 	= av[name][model]['Vars'],
								model 		= model,
								jobID		= jobIDs[model],
								year		= years[model],
								workingDir 	= folder(postprocFolder+name),
								region 		= region,
								grid 		= grid)
							
			#####
			# makePlots:
			# Make some plots of the point to point datasets.
			# MakePlot runs a series of analysis, comparing every pair in DataVars and ModelVars
			#	 under a range of different masks. For instance, only data from Antarctic Ocean, or only data from January.
			# The makePlot produces a shelve file in workingDir containing all results of the analysis.
			m = makePlots(	b.MatchedDataFile, 
					b.MatchedModelFile, 
					name, 
					model, 
					region 		= region,
					year 		= years[model], 
					plotallcuts	= plotallcuts, 
					shelveDir 	= folder(postprocFolder+name+region+'/shelves'),
					imageDir	= folder(imageFolder +'P2Pplots/'+years[model]+'/'+name+region),
					compareCoords	= True )

			shelvesAV[model][name.replace(region,'')][region] = m.shelvesAV
			
			#####
			# make CSV from Shelves:
			# Make some CSV files containing the results of the point to point plots.
			# the option "check" allows to test if a line already exists in the output file. 
			# Alternatives are "over write" and "append", neither are recommnended here.
			
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

						
			
			#####
			# Ocean and month targets for this particular dataset.
			MonthShelves = []
			OceanShelves = []
			for newSlice in m.shelvesAV.keys(): 
			   for xkey in m.shelvesAV[newSlice].keys():
				for ykey in m.shelvesAV[newSlice][xkey].keys():        	      
				    shelve = m.shelvesAV[newSlice][xkey][ykey]			  
				    if newSlice in month_name: 	MonthShelves.append(shelve)
				    if newSlice in Ocean_names:	OceanShelves.append(shelve)
			if len(MonthShelves):	    
			  	filename = folder(imageFolder+'/Targets/'+years[model]+'/Months')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_Months.png'
				makeTargets(	MonthShelves, 
						filename,
						legendKeys = ['newSlice',],					
						)
			if len(OceanShelves):	
				filename = folder(imageFolder+'/Targets/'+years[model]+'/Oceans')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_Oceans.png'
				makeTargets(	OceanShelves, 
						filename,
						legendKeys = ['newSlice',],					
						)

		#####				
		# Here are some fields for comparing fields in the same model
		Summary= {}		
		Summary['MaredatAll'] = []
		Summary['MaredatStandard'] = []			
		Summary['WOAAll'] = []
		Summary['WOAStandard'] = []	

		Summary['SurfaceMetricsAll'] = []	
		Summary['SurfaceMetricsStandard'] = []			
		surfacemetrics = ['chl', 'pCO2', 'nitrate',]
							
		for name in shelvesAV[model].keys():
		  for region in shelvesAV[model][name].keys():
		    for newSlice in shelvesAV[model][name][region].keys(): 
		      for xkey in shelvesAV[model][name][region][newSlice].keys():
			for ykey in shelvesAV[model][name][region][newSlice][xkey].keys():        	      
			  	shelve = shelvesAV[model][name][region][newSlice][xkey][ykey]
				if name in MaredatTypes:
	        		  	if newSlice == 'All':		Summary['MaredatAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['MaredatStandard'].append(shelve)
	        		if name in WOATypes:
	        		  	if newSlice == 'All':		Summary['WOAAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['WOAStandard'].append(shelve)	
	        		
	        		if name in surfacemetrics:
	        		  	if newSlice == 'All':		Summary['SurfaceMetricsAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['SurfaceMetricsStandard'].append(shelve)	
	        		  		        				
	        		for woa in ['silicate','nitrate','phosphate','salinity','temperature','iron',]:
	        		   for ns in ['All', 'Standard']:
	        		   	if ns == newSlice and woa == name.lower():
	        		   		try: 	Summary[woa+ns].append(shelve)
	        		   		except:	Summary[woa+ns]= [shelve,]
		for k in Summary.keys():
			filename = folder(imageFolder+'/Targets/'+years[model]+'/Summary')+model+'_'+years[model]+'_'+k+'.png'
			
	  		makeTargets(Summary[k], 
					filename,
					legendKeys = ['name',],
					debug=True)
		
					
	print "Working dir:",postprocFolder
	

			
			


	
	
	
if __name__=="__main__":

	"""
	./maverick/lib/model_run.py documentation
	
	Module for the Run class.  Each run object holds information about a
	model run, primarily for passing to the assessment area code.

	Keys available in run dictionary for use by assessment code:
	    General
	    -------
	    runid            : 5-character UM runid
	    summary_file     : RMS stats from valnote (csv)
	    radiation_table  : rad means from valnote (plain text)
	    ocean_model      : UM or NEMO
	    run_type         : type of run:
		                      PDC = present-day control
		                      PIC = pre-industrial control
		                      AMIP = forced-SST atmos-only

	    Obs and ancillary data
	    ----------------------
	    clim_root        : Observational climatologies (/project/cma/clim)
	    ancil_root       : Ancillary input files (/project/cma/ancil)

	    Data locations (stashsplit directories except supermeans)
	    ---------------------------------------------------------
	    data_root        : top dir for most data
	    supermean_root   : suprmn filenames up to season
	    ss_daily         : daily means
	    ss_monthly       : monthly means
	    ss_seasonal      : seasonal means
	    ss_annual        : annual means
	    ss_spatiot       : diags for spatio-temp variability

	    Dates specifying standard periods for assessment
	    Values are datetime.date objects.
	    ------------------------------------------------
	    from_daily       : period start for daily data
	    to_daily         : period end for daily data
	    from_monthly     : period start for monthly data
	    to_monthly       : period end for monthly data
	    from_seasonal    : period start for seasonal data
	    to_seasonal      : period end for seasonal data
	    from_annual      : period start for annual data
	    to_annual        : period end for annual data

	    The 'Runs' class is a dictionary of 'Run' instances.
	    Both Run and Runs have an id attribute, for file/dir naming purposes.
	    The id attr for Run instances is the runid; for the Runs group, it is
	    built from the id attrs of the indvidual Run instances.

	    Note difference from areas definition.  The 'runs' section in the config
	    is of the form name=runid.  Name is e.g. 'control', 'exper'.
	    The run is selected from the run_defs config using runid, not name.

	"""
	

	run = {}
	# Info about the run.
	# Commented out those parts of the run dictionary that aren't needed.
	
	run['runid']=		'xhonp' 
	#run['run_type']=	'AMIP' 		#AMIP = forced-SST atmos-only
	run['ocean_model']=	'NEMO'  	# Or MEDUSA
	run['grid']=		'ORCA1'		# ORCA1 or ORCA025
	# Various start and end times:
	#run['start']=		'1982.0' 
	#run['start_year']=	'1982.0' 
	run['end_year']=	'1998'
	#run['nyear']=		'25.0' 
	#run['from_annual']=	'1981-12-01'
	#run['from_daily']=	#'1981-12-01'
	#run['from_monthly']=	'1981-12-01'
	#run['from_seasonal']=	#'1981-12-01'
	#run['to_annual']=	'2005-12-01'
	#run['to_daily']=	'1986-11-30'
	#run['to_monthly']=	'2006-11-01'
	#run['to_seasonal']=	'2006-09-01'
	
	
	
	
	# Data locations:
	run['ss_annual']=	'/group_workspaces/jasmin/esmeval/example_data/bgc/NEMO/'  	# Assume that there is one file per year
	#run['ss_daily']=	'' 
	#run['ss_monthly']=	'/group_workspaces/jasmin/esmeval/example_data/bgc/NEMO/'  	# Assume that there is one file per month
	#run['ss_seasonal']=	'/group_workspaces/jasmin/esmeval/example_data/autoassess/model_data/amzgg/amzgg.splitlev.000003' 
	#run['ss_spatiot']=	'' 
	#run['ancil_root']=	'/project/cma/ancil' 
	#run['clim_root']=	'/project/cma/clim' 
	run['data_root']=	'/group_workspaces/jasmin/esmeval/example_data/bgc' 
	
	
	# Specific files needed to output for AutoAssess
	#run['metrics_model_uncertainty']='amzgg_m_stdev' 
	#run['radiation_table']=		 '/home/users/simonread/dev_workspace/CMA_rev1097_sandbox/example_output/valnote_results/antia_v_amzgg/rad.txt' 	
	run['summary_file']=		 'CSV/summary_global.csv' 
	#run['supermean_root']=		 '/group_workspaces/jasmin/esmeval/example_data/autoassess/model_data/amzgg/supermeans/amzgga.ms2006' 
	
	
	# This is a new one that doesn't exist in the current set up:
	run['LOCALDATA'] = '/group_workspaces/jasmin/esmeval/example_data/bgc'




	# Can use command line arguments to choose a model.
	#models 	= []
	#years 	= []
	#jobIDs 	= []
	
	#####
	# Determine command line arguments
	#for a in argv[1:]:	
	#	try:	
	#		y = int(a)
	#		years.append(a)
	#		continue			
	#	except:pass
		
	#	if str(a).upper() in ['MEDUSA','NEMO']:
	#		models.append(str(a).upper())
	#		continue			
			
	#	print "Command line argument not understood:",a
	

	#####
	#Set Defaults:
	#if not len(years): 	years = ['1998',]
	#if not len(models): 	models = ['MEDUSA','NEMO']

	#print "#############################"
	#print "__main__ arguments: "
	#print "models:        ",models
	#print "year:          ",years
	#print "#############################"

	
	
#	for year in years:
	testsuite_AutoAssess(run)
	
	print 'The end.'
	









