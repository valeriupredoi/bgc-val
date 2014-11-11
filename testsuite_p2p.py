#!/usr/bin/ipython 
#Standard Python modules:
from sys import argv
from os.path import exists
from calendar import month_name

#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict
from p2p import matchDataAndModel,makePlots,makeTargets

from pftnames import MaredatTypes,WOATypes,Ocean_names

###	Potential problems?
###		Reliance on ORCA1 grid



def testsuite_p2p(	models=['MEDUSA','ERSEM','NEMO'],
			year=1998,
			ERSEMjobID='xhonp',
			plotallcuts = False,):

	#####
	# Can use command line arguments to choose a model.
	#if len(argv[1:]): models  = argv[1:]
	#else:	models = ['MEDUSA','ERSEM','NEMO']
    	
    	#####
    	# Which jobs to look at. 
	#ERSEMjobID = 'xhonp'
	jobIDs={}
	jobIDs['ERSEM'] 	= ERSEMjobID
	jobIDs['NEMO'] 		= ERSEMjobID
	jobIDs['MEDUSA'] 	= 'iMarNet'
	
	#####
	# Plot p2p for all regions/oceans, or just everything and "standard" cuts.
	#plotallcuts = True
	
	#####
	# Which Year to investigate for each model.
	# In an ideal world, they would all be the same, except that my current run is stuck in the queue.
	year = str(year)	
	years = {m:year for m in ['MEDUSA','ERSEM','NEMO']}
	#for 
	#years['NEMO']	= year
	#years['ERSEM']	= year	
	#years['MEDUSA']	= year
#	years['NEMO']	= '1893'
#	years['ERSEM']	= '1893'	
#	years['MEDUSA']	= '2007'
		


	
	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	
	#####
	# Location of model files.	
	MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobIDs['ERSEM']+'/'+years['ERSEM']+'/'+jobIDs['ERSEM']+'_'+years['ERSEM']
	NEMOFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobIDs['NEMO'] +'/'+years['NEMO'] +'/'+jobIDs['NEMO'] +'_'+years['NEMO']
	

	
	#####
	# Which analysis to run
	doCHL 		= True
	doMAREDAT 	= True
	doNPSF		= True
	doSalTemp	= True
	doMLD		= True

	
	#####
	# AutoVivification is a form of nested dictionary.
	# we use av here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	if doCHL:
		av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
		av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['chl']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMMisc.nc'			
		av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
		av['chl']['MEDUSA']['Vars'] 		= ['CHL',]	
		av['chl']['ERSEM']['Vars'] 		= ['chl',]
		av['chl']['regions'] 			= ['',]
		
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

	if doNPSF:
		av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
		av['nitrate']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'	
		av['nitrate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
		av['nitrate']['ERSEM']['Vars'] 		= ['N3n','N4n',]
		av['nitrate']['MEDUSA']['Vars'] 	= ['DIN',]									
		av['nitrate']['regions'] 		= ['Surface','100m','200m','500m',]

		av['silicate']['Data']['File'] 		= WOAFolder+'silicate_monthly_1deg.nc'	
		av['silicate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
		av['silicate']['Data']['Vars'] 		= ['i_an',] 		#l+'_mn',
		av['silicate']['ERSEM']['Vars'] 	= ['N5s',]
		av['silicate']['regions'] 		= ['Surface','100m','200m','500m',]
		av['silicate']['MEDUSA']['Vars'] 	= ['SIL',]									
		av['silicate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
			
		av['phosphate']['Data']['File'] 	= WOAFolder+'phosphate_monthly_1deg.nc'	
		av['phosphate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
		av['phosphate']['Data']['Vars'] 	= ['p_an',] 		#l+'_mn',
		av['phosphate']['ERSEM']['Vars'] 	= ['N1p',]
		av['phosphate']['regions'] 		= ['Surface','100m','200m','500m',]		
					
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
		av['salinity']['regions'] 		= ['Surface','500m','100m','200m','1000m',]	 

		av['temperature']['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
		av['temperature']['NEMO']['File'] 	= NEMOFolder+'_NEMO.nc'	
		av['temperature']['Data']['Vars'] 	= ['t_an',]	
		av['temperature']['NEMO']['Vars'] 	= ['votemper',]
		av['temperature']['regions'] 		= ['Surface','500m','100m','200m','1000m',]	
						
				   
	if doMLD:	
		av['mld']['Data']['File'] 		= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		av['mld']['NEMO']['File'] 		= NEMOFolder+'_NEMO.nc'			
		av['mld']['Data']['Vars'] 		= ['mld','mask',]
		av['mld']['NEMO']['Vars'] 		= ['somxl010',]	
		av['mld']['regions'] 			= ['',]

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
	
	
	AutoVivToYaml(av,folder('yaml')+'P2P_Settings.yaml')	
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
				print "testsuite_p2p.py:\tWARNING:\tFile does not exist\tav[",name,"][",model,'][File]'
				continue			    	
			try:
			    if not exists(av[name][model]['File']):
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
								workingDir 	= folder(workingDir+name),
								region 		= region)
							
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
					workingDir 	= folder(workingDir+name+region),
					imageDir	= imageDir,
					compareCoords	=True)
			#Get an autoviv of the shelves.
			shelvesAV[model][name.replace(region,'')][region] = m.shelvesAV
								
										
			#####
			# makeTargets:
			# Make a target diagram of all matches for this particular dataset. # not a great idea if plotAllcuts == True
			filename = folder(imageFolder+'/Targets/'+years[model]+'/AllSlices')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'.png'
			#if model=='ERSEM':filename = filename.replace('ERSEM','ERSEM-'+ERSEMjobID)			
			t = makeTargets(	m.shelves, 
							filename,
							#name=name,
							legendKeys = ['newSlice','ykey',],
							debug=True)
							#imageDir='', 
						
			
			#####
			# Ocean and month targets for this particular dataset.
			Months = []
			Oceans = []
			for newSlice in m.shelvesAV.keys(): 
			   for xkey in m.shelvesAV[newSlice].keys():
				for ykey in m.shelvesAV[newSlice][xkey].keys():        	      
				    shelve = m.shelvesAV[newSlice][xkey][ykey]			  
				    if newSlice in month_name: 	Months.append(shelve)
				    if newSlice in Ocean_names:	Oceans.append(shelve)
				    
	          	filename = folder(imageFolder+'/Targets/'+years[model]+'/Months')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_Months.png'
			#if model=='ERSEM':filename = filename.replace('ERSEM','ERSEM-'+ERSEMjobID)				          	
			makeTargets(	Months, 
					filename,
					legendKeys = ['newSlice',],					
					)
					
			filename = folder(imageFolder+'/Targets/'+years[model]+'/Oceans')+model+'-'+jobIDs[model]+'_'+years[model]+'_'+name+region+'_Oceans.png'
			#if model=='ERSEM':filename = filename.replace('ERSEM','ERSEM-'+ERSEMjobID)						
			makeTargets(	Oceans, 
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
	        		  	
	        		for woa in ['silicate','nitrate','phosphate','salinity','temperature','iron',]:
	        		   for ns in ['All', 'Standard']:
	        		   	if ns == newSlice and woa == name.lower():
	        		   		try: 	Summary[woa+ns].append(shelve)
	        		   		except:	Summary[woa+ns]= [shelve,]
		for k in Summary.keys():
			filename = folder(imageFolder+'/Targets/'+years[model]+'/Summary')+model+'_'+years[model]+'_'+k+'.png'
			#if model=='ERSEM':filename = filename.replace('ERSEM','ERSEM-'+ERSEMjobID)
	  		makeTargets(Summary[k], 
						filename,
						legendKeys = ['name',],#'newSlice',
						debug=True)#imageDir='', diagramTypes=['Target',]
		if model=='ERSEM':
			AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+ERSEMjobID+'.yaml')
		else:
			AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+'.yaml')				
	#####
	# Here are some fields for comparing fields between models.
	# Only works if the same "name" is matched with multiple models.
	
	modelIC = {} 	#model intercomparison shelve files dictionary
		
	for model in shelvesAV.keys():	
	  for name in shelvesAV[model].keys():
	    for region in shelvesAV[model][name].keys():
	      for newSlice in shelvesAV[model][name][region].keys(): 
	        for xkey in shelvesAV[model][name][region][newSlice].keys():
		  for ykey in shelvesAV[model][name][region][newSlice][xkey].keys():
	 	    
			shelve = shelvesAV[model][name][region][newSlice][xkey][ykey]			
	    	   	for ns in ['All', 'Standard']:
	    	   		if ns != newSlice: continue
	    	   		if name in MaredatTypes:
	    	   			try:	modelIC['Maredat'+ns].append(shelve)
	    	   			except: modelIC['Maredat'+ns] = [shelve,]
	        	   		try: 	modelIC[name+ns].append(shelve)
	        	   		except:	modelIC[name+ns]= [shelve,]	    	   			
    	   			if name in WOATypes:
	    	   			try:	modelIC['WOA'+ns].append(shelve)
	    	   			except: modelIC['WOA'+ns] = [shelve,]
	        	   		try: 	modelIC[name+ns].append(shelve)
	        	   		except:	modelIC[name+ns]= [shelve,]
				if name in ['iron',]:
	    	   			try:	modelIC['Misc'+ns].append(shelve)
	    	   			except: modelIC['Misc'+ns] = [shelve,]
	        	   		try: 	modelIC[name+ns].append(shelve)
	        	   		except:	modelIC[name+ns]= [shelve,]					
	

				        	   		
	for k in modelIC.keys():
		for i in modelIC[k]:print k, i
		if len(modelIC[k]) <2: continue
		
		makeTargets(modelIC[k], 
					folder(imageFolder+'/ModelIntercomparison/Targets/')+'intercomparison_'+k+'.png',
					legendKeys = ['xtype','name',],					
					)		
			
			
			

					
	print "Working dir:",workingDir
	
def multiERSEMtargets(shelvesAVs,ERSEMjobIDs):
	#####
	# Here are some fields for comparing fields ERSEM runs.
	print "multiERSEMtargets",ERSEMjobIDs

	shAVs = AutoVivification()	
	modelIC = {} 	#model intercomparison shelve files dictionary
	
	for e,avfn in zip(ERSEMjobIDs,shelvesAVs):
		shAVs[e] = YamlToDict(avfn)
		print "loaded:",e,av
		
	for jobID in shAVs.keys():
	  for model in shAVs[jobID].keys():	
	    for name in shAVs[jobID][model].keys():
	      for region in shAVs[jobID][model][name].keys():
	        for newSlice in shAVs[jobID][model][name][region].keys(): 
	          for xkey in shAVs[jobID][model][name][region][newSlice].keys():
		    for ykey in shAVs[jobID][model][name][region][newSlice][xkey].keys():
	 	    
			shelve = shAVs[jobID][model][name][region][newSlice][xkey][ykey]			
	    	   	for ns in ['All', 'Standard']:
	    	   		if ns != newSlice: continue
	    	   		if name in MaredatTypes:
	    	   			try:	modelIC['Maredat'+ns].append(shelve)
	    	   			except: modelIC['Maredat'+ns] = [shelve,]
	        	   		try: 	modelIC[name+ns].append(shelve)
	        	   		except:	modelIC[name+ns]= [shelve,]	    	   			
    	   			if name in WOATypes:
	    	   			try:	modelIC['WOA'+ns].append(shelve)
	    	   			except: modelIC['WOA'+ns] = [shelve,]
	        	   		try: 	modelIC[name+ns].append(shelve)
	        	   		except:	modelIC[name+ns]= [shelve,]
				if name in ['iron',]:
	    	   			try:	modelIC['Misc'+ns].append(shelve)
	    	   			except: modelIC['Misc'+ns] = [shelve,]
	        	   		try: 	modelIC[name+ns].append(shelve)
	        	   		except:	modelIC[name+ns]= [shelve,]					
	

				        	   		
	for k in modelIC.keys():
		for i in modelIC[k]:print k, i
		if len(modelIC[k]) <2: continue
		
		makeTargets(modelIC[k], 
					folder(imageFolder+'/ModelIntercomparison/Targets/')+'intercomparison_'+k+'.png',
					legendKeys = ['xtype','name',],					
					)		
			
			
			


	
	
	
if __name__=="__main__":
	# Can use command line arguments to choose a model.
	models 		= []
	years 		= []
	ERSEMjobIDs 	= []
	
	#####
	# Determine command line arguments
	for a in argv[1:]:	
		try:	
			y = int(a)
			years.append(a)
			continue			
		except:pass
		
		if str(a).upper() in ['MEDUSA','ERSEM','NEMO']:
			models.append(str(a).upper())
			continue			
		if a[:4] in ['xhon','xjez']:
			ERSEMjobIDs.append(a)
			if 'ERSEM' not in models:models.append('ERSEM')
			continue
			
		print "Command line argument not understood:",a
	

	#####
	#Set Defaults:
	if not len(years): 	years = ['1998',]
	if not len(models): 	models = ['MEDUSA','ERSEM','NEMO']
	if not len(ERSEMjobIDs):ERSEMjobIDs = ['xhonp',]	

	print "#############################"
	print "__main__ arguments: "
	print "models:        ",models
	print "year:          ",years
	print "ERSEM jobID:   ",ERSEMjobIDs
	print "#############################"
	#sleep(20)


			
	for year in years:
		#multiERSEMtargets(['yaml/shelvesAVERSEM'+year+e+'.yaml' for e in ERSEMjobIDs],ERSEMjobIDs)	
		#continue
		testsuite_p2p(models = models,	year=year,ERSEMjobID=ERSEMjobIDs[0] ) 
		if len(ERSEMjobIDs)==1:continue
		for e in ERSEMjobIDs[1:]:
			testsuite_p2p(models = ['ERSEM',],year=year,ERSEMjobID=e ) 
		multiERSEMtargets(['shelvesAVERSEM'+year+e+'.yaml' for e in ERSEMjobIDs],ERSEMjobIDs)
	
	print 'The end.'
	



















	
