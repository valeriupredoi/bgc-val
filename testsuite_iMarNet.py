#!/usr/bin/ipython 
#
# Copyright 2014 Plymouth Marine Laboratory
#
# This file is part of the ukesm-validation library.
#
# ukesm-validation is free software: you can redistribute it and/or modify it
# under the terms of the Lesser GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version. 

# ukesm-validation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU
# General Public License for more details.
# You should have received a copy of the Lesser GNU General
# Public License along with ukesm-validation. If not, see <http://www.gnu.org/licenses/>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.ukesm
#

#Standard Python modules:
from sys import argv
from os.path import exists
from calendar import month_name

#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict
from p2p import matchDataAndModel,makePlots,makeTargets, csvFromShelves
#from 
from pftnames import MaredatTypes,WOATypes,Ocean_names,getmt

###	Potential problems?
###		Reliance on ORCA1 grid
from shelve import open as shOpen

# Will not run this code in MO.
from ncdfView import ncdfView

def testsuite_iMarNet(	models=['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10'],
			year=1998,
			#ERSEMjobID='xhonc',
			plotallcuts = False,):

	
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
	doCHL 		= 0#True
	doN		= 0#True
	doP		= 0#False
	doSi		= 0#True
	doFe		= 0#True		
	doPCO2		= 0#True
	doIntPP		= True

	#doPSF		= 0#True	
	#doSalTemp	= 0#True
	#doMLD		= 0#True
			
	#####
	# AutoVivification is a form of nested dictionary.
	# we use av here to determine which files to analyse and which fields in those files.
	# Region is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	modelSkips = {	'Diat-HadOCC':	[],
			'ERSEM':	[],
			'HadOCC':	[],
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
		av['chl']['regions'] 			= regions
		
	if doN:
		av['nitrate']['Data']['File'] 		=  WOAFolder+'nitrate_monthly_1deg.nc'		
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 
		for m in models:
			if 'chl' in modelSkips[m]:continue		
			av['nitrate'][m]['Vars'] 	= ['no3',]						
			av['nitrate'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
		av['nitrate']['regions'] 		= regions

	if doP:
		av['phosphate']['Data']['File'] 	=  WOAFolder+'phosphate_monthly_1deg.nc'		
		av['phosphate']['Data']['Vars'] 	= ['p_an',] 
		for m in models:
			if 'chl' in modelSkips[m]:continue		
			av['phosphate'][m]['Vars'] 	= ['po4',]						
			av['phosphate'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
		av['nitrate']['regions'] 		= regions
				
	if doSi:
		av['silicate']['Data']['File'] 		=  WOAFolder+'silicate_monthly_1deg.nc'		
		av['silicate']['Data']['Vars'] 		= ['i_an',] 
		for m in models:
			if 'chl' in modelSkips[m]:continue
			av['silicate'][m]['Vars'] 	= ['si',]
			av['silicate'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
		av['silicate']['regions'] 		= regions
		
	if doIntPP:
		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'		
		av['intpp']['Data']['Vars'] 	= ['PPint',]
		#av['pp']['Data']['File'] 	=  MAREDATFolder+'PP100108.nc'		
		#av['pp']['Data']['Vars'] 	= ['PP',] 
			 
		for m in models:
			av['intpp'][m]['Vars'] 	= ['intpp',]						
			av['intpp'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
		av['intpp']['regions'] 		= regions
		
	if doFe:
					
		av['iron']['Data']['File'] 	=  GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		av['iron']['Data']['Vars'] 	= ['Fe_D_CONC_BOTTLE',] 
		for m in models:
			av['iron'][m]['Vars'] 	= ['dfe',]						
			av['iron'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
		av['iron']['regions'] 		= ['',]
	
	if doPCO2:
		av['pCO2']['Data']['File'] 	=  TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'		
		av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 
		for m in models:
			av['pCO2'][m]['Vars'] 	= ['spco2',]						
			av['pCO2'][m]['File']	= iMarNetFolder+iMarNetFiles[m]
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
			imageFolder 	= folder('images/iMarNet/'+model+'-'+jobIDs[model])
			workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/iMarNet/"+model+'-'+jobIDs[model]+'-'+years[model])		

		
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
			b = matchDataAndModel(av[name]['Data']['File'], 
								av[name][model]['File'],
								name,
								DataVars  	= av[name]['Data']['Vars'],
								ModelVars 	= av[name][model]['Vars'],
								model 		= model,
								jobID		= jobIDs[model],
								year		= years[model],
								workingDir 	= folder(workingDir+name+region),
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
					#model, 
					jobID 		= 'IMARNET_'+model,
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

		Summary['AllAll'] = []	
		Summary['AllStandard'] = []			
		Summary['SurfaceMetricsAll'] = []	
		Summary['SurfaceMetricsStandard'] = []			
		surfacemetrics = ['chl', 'pCO2', 'nitrate',]
							
		for name in shelvesAV[model].keys():
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
		#if model=='ERSEM':
		AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+jobIDs[model]+'.yaml')
		#else:
		#	AutoVivToYaml(shelvesAV, folder('yaml')+'shelvesAV'+model+years[model]+'.yaml')				

			

					
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
	



















	
