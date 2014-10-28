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



def testsuite_p2p():

	#####
	# Can use command line arguments to choose a model.
	if len(argv[1:]): models  = argv[1:]
	else:	models = ['MEDUSA','ERSEM','NEMO']
    	
    	#####
    	# Which ERSEM job to look at. 
	ERSEMjobID = 'xhonc'

	#####
	# Which Year to investigate for each model.
	# In an ideal world, they would all be the same, except that my current run is stuck in the queue.
	years = {}
	years['NEMO']	= '1893'
	years['ERSEM']	= '1893'	
	years['MEDUSA']	= '2007'
	


	
	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	
	#####
	# Location of model files.	
	MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ ERSEMjobID+'/'+years['ERSEM']+'/'+ERSEMjobID+'_'+years['ERSEM']
	NEMOFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ ERSEMjobID+'/'+years['NEMO'] +'/'+ERSEMjobID+'_'+years['NEMO']
	
	#####
	# Location of image Output files
	imageFolder 	= folder('images')
	
	#####
	# Which analysis to run
	doCHL 		= True
	doMAREDAT 	= True
	doSalTemp	= True
	doMLD		= True
	doNPSF		= True
	
	#####
	# AutoVivification is a form of nested dictionary.
	# we use av here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	
	av = AutoVivification()
	if doCHL:
		av['chl']['Data'  ]['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
		av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['chl']['ERSEM' ]['File'] 		= ERSEMFolder+'_ERSEMMisc.nc'			
		av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
		av['chl']['MEDUSA']['Vars'] 		= ['CHL',]	
		av['chl']['ERSEM']['Vars'] 		= ['chl',]
		av['chl']['region'] 			= ''
		
	if doMAREDAT:
		av['diatoms']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
		av['diatoms']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['diatoms']['ERSEM' ]['File'] 	= ERSEMFolder+'_ERSEMphytoBm.nc'				
		av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
		av['diatoms']['MEDUSA']['Vars'] 	= ['PHD',]	
		av['diatoms']['ERSEM']['Vars'] 		= ['P1c',]
		av['diatoms']['region'] 		= ''	
	
		av['bac']['Data'  ]['File'] 		= MAREDATFolder+"MarEDat20120214Bacteria.nc"	
		av['bac']['ERSEM' ]['File'] 		= ERSEMFolder+'_ERSEMbac.nc'			
		av['bac']['Data']['Vars'] 		= ['BIOMASS',]
		av['bac']['ERSEM']['Vars'] 		= ['B1c',]
		av['bac']['region'] 			= ''	
	
		av['picophyto']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20111206Picophytoplankton.nc"	
		av['picophyto']['ERSEM' ]['File'] 	= ERSEMFolder+'_ERSEMphytoBm.nc'			
		av['picophyto']['Data']['Vars'] 	= ['BIOMASS',]
		av['picophyto']['ERSEM']['Vars'] 	= ['P3c',]
		av['picophyto']['region'] 		= ''	
		
		av['microzoo']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
		av['microzoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['microzoo']['ERSEM' ]['File'] 	= ERSEMFolder+'_ERSEMzoo.nc'			
		av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['microzoo']['MEDUSA']['Vars'] 	= ['ZMI',]	
		av['microzoo']['ERSEM']['Vars'] 	= ['Z5c',]
		av['microzoo']['region'] 		= ''	
	
		av['mesozoo']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
		av['mesozoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['mesozoo']['ERSEM' ]['File'] 	= ERSEMFolder+'_ERSEMzoo.nc'			
		av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['mesozoo']['MEDUSA']['Vars'] 	= ['ZME',]	
		av['mesozoo']['ERSEM']['Vars'] 		= ['Z4c',]
		av['mesozoo']['region'] 		= ''

	if doNPSF:
		for woa in ['nitrate','phosphate','silicate',]:
			if woa == 'silicate':	
				l='i' 
				ERSEMVars  	= ['N5s',]
				MEDVars		= ['SIL',]
			elif woa == 'nitrate':	
				l='n' 
				ERSEMVars  	= ['N3n','N4n',]
				MEDVars		= ['DIN',]
			elif woa == 'phosphate':	
				l='p' 
				ERSEMVars  	= ['N1p',]
		
			for s in ['Surface','100m','200m','500m',]:#'Transect',]:#'All',
				av[woa+s]['Data'  ]['File'] 	= WOAFolder+woa+'_monthly_1deg.nc'	
				av[woa+s]['ERSEM' ]['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
				av[woa+s]['Data']['Vars'] 	= [l+'_an',] 		#l+'_mn',
				av[woa+s]['ERSEM']['Vars'] 	= ERSEMVars
				av[woa+s]['region'] 		= s	    
				if woa != 'phosphate':
					av[woa+s]['MEDUSA']['Vars'] 	= MEDVars									
					av[woa+s]['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"
					
		av['iron']['Data'  ]['File'] 	= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		av['iron']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+years['MEDUSA']+".nc"	
		av['iron']['ERSEM' ]['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'			
		av['iron']['Data']['Vars'] 	= ['Fe_D_CONC_BOTTLE',]
		av['iron']['MEDUSA']['Vars'] 	= ['FER',]	
		av['iron']['ERSEM']['Vars'] 	= ['N7f',]
		av['iron']['region'] 		= ''
		
	if doSalTemp:
		for woa in ['salinity','temperature',]:
			if woa == 'salinity':		NEMOVars  	= ['vosaline',]
			if woa == 'temperature':	NEMOVars  	= ['votemper',]
			for s in ['Surface','500m','100m','200m','1000m',]:
				av[woa+s]['Data'  ]['File'] 	= WOAFolder+woa+'_monthly_1deg.nc'	
				av[woa+s]['NEMO' ]['File'] 	= NEMOFolder+'_NEMO.nc'	
				av[woa+s]['Data']['Vars'] 	= [woa[0]+'_an',]	#woa[0]+'_mn',
				av[woa+s]['NEMO']['Vars'] 	= NEMOVars
				av[woa+s]['region'] 		= s	 
				   
	if doMLD:	
		av['mld']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		av['mld']['NEMO' ]['File'] 	= NEMOFolder+'_NEMO.nc'			
		av['mld']['Data']['Vars'] 	= ['mld','mask',]
		av['mld']['NEMO']['Vars'] 	= ['somxl010',]	
		av['mld']['region'] 		= ''

		av['mld_DR003']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		av['mld_DR003']['NEMO' ]['File'] 	= NEMOFolder+'_NEMO.nc'			
		av['mld_DR003']['Data']['Vars'] 	= ['mld','mask',]
		av['mld_DR003']['NEMO']['Vars'] 	= ['somxl010',]	
		av['mld_DR003']['region'] 		= ''

		av['mld_DReqDTm02']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
		av['mld_DReqDTm02']['NEMO' ]['File'] 	= NEMOFolder+'_NEMO.nc'			
		av['mld_DReqDTm02']['Data']['Vars'] 	= ['mld','mask',]
		av['mld_DReqDTm02']['NEMO']['Vars'] 	= ['somxl010',]	
		av['mld_DReqDTm02']['region'] 		= ''
	
	
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
		#for name in ['chl',]:	
			#####
			# Do some checks to make sure that the files all exist:
			print model,name
			try: 
				if not isinstance(av[name][model],dict): continue
			except KeyError:
				print "No ",name, 'in ',model
				continue
				
			region = str(av[name]['region'])
			
			workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/"+model+'-'+years[model])		
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
								DataVars  = av[name]['Data']['Vars'],
								ModelVars = av[name][model]['Vars'],
								jobID=model,
								year=years[model],
								workingDir = folder(workingDir+name),
								region = region)
								
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
							year = years[model], 
							plotallcuts=False, 
							workingDir = folder(workingDir+name),
							compareCoords=True)
			#Get an autoviv of the shelves.
			shelvesAV[model][name.replace(region,'')][region] = m.shelvesAV
										
			#####
			# makeTargets:
			# Make a target diagram of all matches for this particular dataset. # not a great idea if plotAllcuts == True
			filename = folder(imageFolder+model+'/'+years[model]+'/Targets/')+model+'_'+years[model]+'_'+name+'.png'
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
	          
			makeTargets(Months, 
						folder(imageFolder+model+'/Targets/Months')+model+'_'+years[model]+'_'+name+'_Months.png',
						legendKeys = ['newSlice',],					
						)									
			makeTargets(Oceans, 
						folder(imageFolder+model+'/Targets/'+years[model]+'/Oceans')+model+'_'+years[model]+'_'+name+'_Oceans.png',
						legendKeys = ['newSlice',],					
						)
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
			filename = folder(imageFolder+model+'/Targets/'+years[model]+'/Summary')+model+'_'+years[model]+'_'+k+'.png'
	  		makeTargets(Summary[k], 
						filename,
						legendKeys = ['name', ],#'newSlice',
						debug=True)#imageDir='', diagramTypes=['Target',]
		AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV.yaml')
				
	#####
	# Here are some fields for comparing fields between models.

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
	

	
if __name__=="__main__":
	testsuite_p2p() 
	print 'The end.'
	
	
