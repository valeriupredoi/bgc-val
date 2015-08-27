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
from UKESMpython import folder,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict
from p2p import matchDataAndModel,makePlots,makeTargets, csvFromShelves, makePatternStatsPlots
#from 
from pftnames import MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names,getmt

###	Potential problems?
###		Reliance on ORCA1 grid



def testsuite_p2p(	models=['MEDUSA','ERSEM','NEMO'],
			year=1998,
			jobID='xhonp',
			plotallcuts = False,):

	#####
	# Can use command line arguments to choose a model.
	#if len(argv[1:]): models  = argv[1:]
	#else:	models = ['MEDUSA','ERSEM','NEMO']
    	
    	#####
    	# Which jobs to look at. 
	#ERSEMjobID = 'xhonp'
	jobIDs={}
	jobIDs['ERSEM'] 	= jobID
	jobIDs['NEMO'] 		= jobID
	jobIDs['MEDUSA'] 	= jobID
	
	#####
	# Plot p2p for all regions/oceans, or just everything and "standard" cuts.
	#plotallcuts = True
	
	#####
	# Which Year to investigate for each model.
	# In an ideal world, they would all be the same, except that my current run is stuck in the queue.
	year = str(year)	
	years = {m:year for m in ['MEDUSA','ERSEM','NEMO']}


	
	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
	#####
	# Location of model files.	
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobIDs['MEDUSA']+'-' + years['MEDUSA']+'/'
	ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobIDs['ERSEM']+'/'+years['ERSEM']+'/'+jobIDs['ERSEM']+'_'+years['ERSEM']
	#NEMOFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobIDs['NEMO'] +'/'+years['NEMO'] +'/'+jobIDs['NEMO'] +'_'+years['NEMO']
	NEMOFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobIDs['MEDUSA']+'-' + years['MEDUSA']+'/'	

	regions = ['Surface','100m','500m',]#'200m'
	
	#####
	# Which analysis to run
	doCHL 		= 0#True
	doDMS_clim	= 0#True
	doDMS_pixels	= True
	doDMS_pixels2	= 0#True
	doMAREDAT 	= 0#True
	doN		= 0#True
	doPSF		= 0#True	
	doSalTemp	= 0#True
	doMLD		= 0#True
	doPCO2		= 0#True
	doLight		= 0
	#####
	# getmt
	#mt = getmt()
	#def getVarsFromMT(var, model):
	#	try:	outs = mt[model][var]['vars']
	#	except: outs = mt[model][var]
	#	#print model,var,outs
	#	if len(outs):return outs
	#	return []
		
		
	#####
	# AutoVivification is a form of nested dictionary.
	# we use av here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	if doCHL:
		av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
		#av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+jobIDs['MEDUSA']+'_' + years['MEDUSA']+"_MEDUSA_bio.nc"
		av['chl']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMMisc.nc'	
		av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
		av['chl']['MEDUSA']['Vars'] 		= ['CHD','CHN']	
		av['chl']['ERSEM']['Vars'] 		= ['chl',]
		av['chl']['regions'] 			= ['',]
		av['chl']['MEDUSA']['grid']		= 'ORCA1'		
		av['chl']['ERSEM']['grid']		= 'ORCA1'				

	if doDMS_clim:
		dmsd= {'dms_and':'anderson','dms_ara':'aranamit','dms_hal':'halloran','dms_sim':'simodach'} 
		for dms in ['dms_and','dms_ara','dms_hal','dms_sim',]:
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'		
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['Data']['Vars'] 		= ['lanaetal',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]
			av[dms]['regions'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] = 'Flat1deg'
		plotallcuts =  True # False#

	if doDMS_pixels:
		dmsd= {'dms_p_and':'anderson','dms_p_ara':'aranamit','dms_p_hal':'halloran','dms_p_sim':'simodach'} 
		for dms in ['dms_p_sim','dms_p_and','dms_p_ara','dms_p_hal',]:
			av[dms]['Data']['File'] 		= '/data/euryale7/scratch/ledm/DMS_Lana2011nc/DMSpixels.nc'		
			av[dms]['MEDUSA']['File'] 		= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrum/dms_xkrum.nc'	
			av[dms]['Data']['Vars'] 		= ['DMS',]
			av[dms]['MEDUSA']['Vars'] 		= [dmsd[dms],]			
			av[dms]['regions'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] = 'Flat1deg'
		plotallcuts =   True # False#

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
			av[dms]['regions'] 			= ['Surface',]	
			av[dms]['MEDUSA']['grid'] = 'Flat1deg'
		plotallcuts =  0#True # False#
		
						
	if doMAREDAT:
		av['diatoms']['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
		av['diatoms']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['diatoms']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMphytoBm.nc'				
		av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
		av['diatoms']['MEDUSA']['Vars'] 	= ['PHD',]	
		av['diatoms']['ERSEM']['Vars'] 		= ['P1c',]
		av['diatoms']['regions'] 		= ['',]	
	
		av['bac']['Data']['File'] 		= MAREDATFolder+"MarEDat20120214Bacteria.nc"	
		av['bac']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMbac.nc'			
		av['bac']['Data']['Vars'] 		= ['BIOMASS',]
		av['bac']['ERSEM']['Vars'] 		= ['B1c',]
		av['bac']['regions'] 			= ['',]
	
		av['picophyto']['Data']['File'] 	= MAREDATFolder+"MarEDat20111206Picophytoplankton.nc"	
		av['picophyto']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMphytoBm.nc'			
		av['picophyto']['Data']['Vars'] 	= ['BIOMASS',]
		av['picophyto']['ERSEM']['Vars'] 	= ['P3c',]
		av['picophyto']['regions'] 		= ['',]	
		
		av['microzoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
		av['microzoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['microzoo']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMzoo.nc'			
		av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['microzoo']['MEDUSA']['Vars'] 	= ['ZMI',]	
		av['microzoo']['ERSEM']['Vars'] 	= ['Z5c',]
		av['microzoo']['regions'] 		= ['',]	
	
		av['mesozoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
		av['mesozoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['mesozoo']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMzoo.nc'			
		av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['mesozoo']['MEDUSA']['Vars'] 	= ['ZME',]	
		av['mesozoo']['ERSEM']['Vars'] 		= ['Z4c',]
		av['mesozoo']['regions'] 		= ['',]

	if doN:
		av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
		av['nitrate']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'	
		#av['nitrate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
		av['nitrate']['MEDUSA']['File'] 	= MEDUSAFolder+jobIDs['MEDUSA']+'_'+ years['MEDUSA']+"_MEDUSA_bio.nc"		
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
		av['nitrate']['ERSEM']['Vars'] 		= ['N3n','N4n',]
		av['nitrate']['MEDUSA']['Vars'] 	= ['DIN',]									
		av['nitrate']['regions'] 		= regions
		av['nitrate']['MEDUSA']['grid']		= 'ORCA1'		
		av['nitrate']['ERSEM']['grid']		= 'ORCA1'				
			
		plotallcuts =  True # False#					
	if doPSF:
		av['silicate']['Data']['File'] 		= WOAFolder+'silicate_monthly_1deg.nc'	
		av['silicate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
		av['silicate']['Data']['Vars'] 		= ['i_an',] 		#l+'_mn',
		av['silicate']['ERSEM']['Vars'] 	= ['N5s',]
		av['silicate']['regions'] 		= regions
		av['silicate']['MEDUSA']['Vars'] 	= ['SIL',]									
		av['silicate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
			
		av['phosphate']['Data']['File'] 	= WOAFolder+'phosphate_monthly_1deg.nc'	
		av['phosphate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
		av['phosphate']['Data']['Vars'] 	= ['p_an',] 		#l+'_mn',
		av['phosphate']['ERSEM']['Vars'] 	= ['N1p',]
		av['phosphate']['regions'] 		= regions		
					
		av['iron']['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		av['iron']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['iron']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'			
		av['iron']['Data']['Vars'] 		= ['Fe_D_CONC_BOTTLE',]
		av['iron']['MEDUSA']['Vars'] 		= ['FER',]	
		av['iron']['ERSEM']['Vars'] 		= ['N7f',]
		av['iron']['regions'] 			= ['',]
		
	if doSalTemp:
		av['salinity']['Data']['File'] 		= WOAFolder+'salinity_monthly_1deg.nc'	
		av['salinity']['NEMO']['File'] 		= NEMOFolder+'_NEMO.nc'	
		av['salinity']['Data']['Vars'] 		= ['s_an',]
		av['salinity']['NEMO']['Vars'] 		= ['vosaline',]
		av['salinity']['regions'] 		= regions	 

		av['temperature']['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
		av['temperature']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'	
		av['temperature']['Data']['Vars'] 	= ['t_an',]	
		av['temperature']['NEMO']['Vars'] 	= ['votemper',]
		av['temperature']['regions'] 		= regions	
				   
	if doMLD:	
		av['mld']['Data']['File'] 		= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['mld']['NEMO']['File'] 		= NEMOFolder+'_NEMO.nc'	
		av['mld']['NEMO']['File'] 		= NEMOFolder+jobIDs['MEDUSA']+'_'+ years['MEDUSA']+"_MEDUSA_MLD.nc"			
		av['mld']['Data']['Vars'] 		= ['mld','mask',]
		av['mld']['NEMO']['Vars'] 		= ['somxl010',]	
		av['mld']['regions'] 			= ['',]
		av['mld']['NEMO']['grid'] 		= 'ORCA1'


	if doLight:	
		av['irradiation']['Data']['File'] 		= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		av['irradiation']['NEMO']['File'] 		= NEMOFolder+jobIDs['MEDUSA']+'_'+ years['MEDUSA']+"_MEDUSA_Light.nc"			
		av['irradiation']['Data']['Vars'] 		= []
		av['irradiation']['NEMO']['Vars'] 		= ['MED_QSR',]	
		av['irradiation']['regions'] 			= ['',]
		av['irradiation']['NEMO']['grid'] 		= 'ORCA1'
		
				
		#av['mld_DR003']['Data']['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['mld_DR003']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'			
		#av['mld_DR003']['Data']['Vars'] 	= ['mld','mask',]
		#av['mld_DR003']['NEMO']['Vars'] 	= ['somxl010',]	
		#av['mld_DR003']['regions'] 		= ['',]

		#av['mld_DReqDTm02']['Data']['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		#av['mld_DReqDTm02']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'			
		#av['mld_DReqDTm02']['Data']['Vars'] 	= ['mld','mask',]
		#av['mld_DReqDTm02']['NEMO']['Vars'] 	= ['somxl010',]	
		#av['mld_DReqDTm02']['regions'] 		= ['',]
		
	if doPCO2:
		av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
		av['pCO2']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMMisc.nc'	
		av['pCO2']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
		av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
		av['pCO2']['ERSEM']['Vars'] 	= ['pCO2w',]
		av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]	
		av['pCO2']['regions'] 		= ['',]
	
	#for var in av.keys():
	#	for m in models: 
	#		#keys = getVarsFromMT(var,m)
	#		av[var][m]['Vars'] = getVarsFromMT(var,m)
	#		#print var,m, keys, av[var][m]['Vars']
	#		#if len(keys) and keys != av[var][m]['Vars']:print "Dont' match."
	#assert False	
	
	#AutoVivToYaml(av,folder('yaml')+'P2P_Settings.yaml')	
	#av = 0
	#print av
	#av = YamlToDict(folder('yaml')+'P2P_Settings.yaml',)
	#print av.keys(), av['chl'].keys(),av['chl']['MEDUSA']
	
	
	
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
			imageFolder 	= folder('images/'+model+'-'+jobIDs[model])
			workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/"+model+'-'+jobIDs[model]+'-'+years[model])		

		
			try:
			    if not exists(av[name]['Data']['File']):
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name]['Data']['File']
				continue
			except:
				print "testsuite_p2p.py:\tWARNING:\tDict entry does not exist\tav[",name,"][",model,'][File]'
				continue			    	
			try:
			    if not exists(av[name][model]['File']):
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name][model]+'[File]'
				continue
			except:
				print "testsuite_p2p.py:\tWARNING:\tDict entry exist:\tav[",name,"][",model,'][File]'
				continue			
			print "\n\n\ntestsuite_p2p.py:\tINFO:\tRunning:",name
			
			
			#####
			# matchDataAndModel:
			# Match (real) Data and Model. 
			# Does not produce and plots.
			try:	grid = av[name][model]['grid']
			except: grid = 'ORCA1'
			b = matchDataAndModel(av[name]['Data']['File'], 
								av[name][model]['File'],
								name,
								DataVars  	= av[name]['Data']['Vars'],
								ModelVars 	= av[name][model]['Vars'],
								model 		= model,
								jobID		= jobIDs[model],
								year		= years[model],
								workingDir 	= folder(workingDir+name),
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
					model, 
					region 		= region,
					year 		= years[model], 
					plotallcuts	= plotallcuts, 
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
			#print shelvesAV.keys(),m.shelvesAV.keys()
			#assert False	    				
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
							months,#xkeysLabels=
							filenamebase,	# filename base						
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
							Ocean_names,#xkeysLabels=
							filenamebase,	# filename base						
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
							months,#xkeysLabels=
							filenamebase,	# filename base						
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
							seasons,	#xkeysLabels=
							filenamebase,	# filename base						
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
							months,			#xkeysLabels=
							filenamebase,		# filename base						
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
		Summary['DMSAll'] = []
		Summary['DMSStandard']=[]	
		Summary['DMS_p_All'] = []
		Summary['DMS_p_Standard']=[]

		MonthsPatterns = {'DMS_p':{},'DMS_e':{},'Maredat':{},'WOA':{},}#'DMS_p':{},}
		OceansPatterns = {'DMS_p':{},'DMS_e':{},'Maredat':{},'WOA':{},}#'DMS_p':{},}		
		
		for name in shelvesAV[model].keys():
   		  for p in MonthsPatterns.keys(): MonthsPatterns[p][name] = []
   		  for p in OceansPatterns.keys(): OceansPatterns[p][name] = []   		  
		  for region in shelvesAV[model][name].keys():
		    for newSlice in shelvesAV[model][name][region].keys(): 
		      for xkey in shelvesAV[model][name][region][newSlice].keys():
			for ykey in shelvesAV[model][name][region][newSlice][xkey].keys():        	      
			  	shelve = shelvesAV[model][name][region][newSlice][xkey][ykey]
	       		  	if newSlice == 'All':		Summary['AllAll'].append(shelve)
	       		  	if newSlice == 'Standard':	Summary['AllStandard'].append(shelve)
				if name in MaredatTypes:
	        		  	if newSlice == 'All':		Summary['MaredatAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['MaredatStandard'].append(shelve)
	        		if name in WOATypes:
	        		  	if newSlice == 'All':		Summary['WOAAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['WOAStandard'].append(shelve)	
	        		
	        		if name in dmsmetrics:
	        		  	if newSlice == 'All':		Summary['DMSAll'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary['DMSStandard'].append(shelve)	
					if newSlice in Ocean_names:	OceansPatterns['DMS_e'][name].append(shelve)
					if newSlice in months:		MonthsPatterns['DMS_e'][name].append(shelve)
					
	        		if name in dmspmetrics:
	        		  	if newSlice == 'All':  		Summary[ 'DMS_p_All'].append(shelve)
	        		  	if newSlice == 'Standard':	Summary[ 'DMS_p_Standard'].append(shelve)
					if newSlice in Ocean_names:	OceansPatterns['DMS_p'][name].append(shelve)
					if newSlice in months:		MonthsPatterns['DMS_p'][name].append(shelve)
										
	        		  	#Patterns['DMS_p'][name].append(shelve)
	        		  	#if newSlice not in patternShelves: patternShelves.append(newSlice)
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
					
		for k in OceansPatterns.keys():					
			filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/OceansPatterns/')+k+'_'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
			print k, OceansPatterns[k],filenamebase
			
			makePatternStatsPlots(	OceansPatterns[k], # {legend, shelves}
					k,	#xkeysname
					Ocean_names,			#xkeysLabels=
					filenamebase,		# filename base						
					)
		for k in MonthsPatterns.keys():					
			filenamebase = folder(imageFolder+'/Patterns/'+years[model]+'/MonthsPatterns/')+k+'_'+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region
			print k, MonthsPatterns[k],filenamebase
			makePatternStatsPlots(	
					MonthsPatterns[k], # {legend, shelves}
					k,	#xkeysname
					months,			#xkeysLabels=
					filenamebase,		# filename base						
					)	
																						
		#if model=='ERSEM':
		AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+jobIDs[model]+'.yaml')
		#else:
		#	AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+'.yaml')				

			

					
	print "Working dir:",workingDir
	

			
			


	
	
	
if __name__=="__main__":

	# Can use command line arguments to choose a model.
	models 		= []
	years 		= []
	jobIDs 		= []
	
	#####
	# Determine command line arguments
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
			if 'ERSEM' not in models:models.append('ERSEM')
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
	



















	
