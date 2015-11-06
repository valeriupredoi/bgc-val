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
from UKESMpython import populateSlicesList, AutoVivification, folder, reducesShelves,listShelvesContents,mnStr
from testsuite_p2p import testsuite_p2p
from p2p import makePatternStatsPlots
from pftnames import months

#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	

# 	from 
def analysis():
	# DMS model:
	#model= 'MEDUSA'
	#jobID = 'xkrum'
	#year = 'clim'		
	#Grid = 'ORCA1'
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobID+'-' + year+'/'

	# DMS model:
	models= ['MEDUSA','NEMO']
	jobID = 'xkrus'
	years = ['2075','2076','2077']
	Grid = 'ORCA1'
	MEDUSAFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	NEMOFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	
	# ORCA025:
	#model= 'MEDUSA'
	#jobID = 'xjwki'
	#year = '1979'		
	#Grid = 'ORCA025'	
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/"+jobID+"_postProc/"+year+"/"
	
		
	#####
	# Which analysis to run
	doCHL 		= 0#True
	doMAREDAT 	= 0#True
	doN		= 0#True
	doSi		= 0#True	
	doFe		= 0#True		
	doPCO2		= 0#True
	doIntPP		= 0#True
	doO2		= True	
	
	doSal		= True
	doTemp		= True
	doMLD		= True
	
	
	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)	
	depthLevels 	= ['Transect','Surface',]#'100m','200m','500m',]#'transect',]#'		

	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
	LesterFolder 	= "/data/euryale7/scratch/ledm/LestersReportData/"			
	
	shelvesAV = []
	
	for year in years:		
		#####
		# Location of model files.
		MEDUSAFolder = MEDUSAFolder_pref+year+'/'
		NEMOFolder  	= NEMOFolder_pref+year+'/'
		#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
		#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobID+'-' + year+'/'
	

		#####
		# AutoVivification is a form of nested dictionary.
		# We use AutoVivification here to determine which files to analyse and which fields in those files.
		# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
		av = AutoVivification()
		if doCHL:
			av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
			if Grid == 'ORCA1':	
				av['chl']['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_MEDUSA_bio.nc"
				av['chl']['MEDUSA']['Vars'] 	= ['CHD','CHN']				
			if Grid == 'ORCA025':
				av['chl']['MEDUSA']['File']	= MEDUSAFolder+"xjwki_1979_CH.nc"
				av['chl']['MEDUSA']['Vars'] 	= ['CHL',]			
			av['chl']['Data']['Vars'] 		= ['Chlorophylla',]

			av['chl']['depthLevels'] 		= ['',]
			av['chl']['MEDUSA']['grid']		= Grid		
						
		if doMAREDAT:
			av['diatoms']['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
			av['diatoms']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"	
			av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
			av['diatoms']['MEDUSA']['Vars'] 	= ['PHD',]	
			av['diatoms']['depthLevels'] 		= ['',]	
			av['diatoms']['MEDUSA']['grid']		= Grid						
		
			av['microzoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
			av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
			av['microzoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"	
			av['microzoo']['MEDUSA']['Vars'] 	= ['ZMI',]	
			av['microzoo']['MEDUSA']['grid']	= Grid		
			av['microzoo']['depthLevels'] 		= ['',]	
	
			av['mesozoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
			av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
			av['mesozoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"	
			av['mesozoo']['MEDUSA']['Vars'] 	= ['ZME',]	
			av['mesozoo']['MEDUSA']['grid']		= Grid		
			av['mesozoo']['depthLevels'] 		= ['',]

		if doN:
			av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
			av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
			if Grid == 'ORCA1':			
				av['nitrate']['MEDUSA']['File'] = MEDUSAFolder+jobID+'_'+ year+"_MEDUSA_bio.nc"		
			if Grid == 'ORCA025':
				av['nitrate']['MEDUSA']['File'] = MEDUSAFolder+jobID+'_'+ year+"_DIN.nc"							
			av['nitrate']['MEDUSA']['Vars'] 	= ['DIN',]									
			av['nitrate']['MEDUSA']['grid']		= Grid		
			av['nitrate']['depthLevels'] 		= depthLevels
			
		if doSi:
			av['silicate']['Data']['File'] 		= WOAFolder+'silicate_monthly_1deg.nc'	
			av['silicate']['Data']['Vars'] 		= ['i_an',] 		#l+'_mn',
			av['silicate']['MEDUSA']['Vars'] 	= ['SIL',]									
			av['silicate']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"
			av['silicate']['MEDUSA']['grid']	= Grid		
			av['silicate']['depthLevels'] 		= depthLevels
						
		if doFe:							
			av['iron']['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
			av['iron']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+year+".nc"	
			av['iron']['Data']['Vars'] 		= ['Fe_D_CONC_BOTTLE',]
			av['iron']['MEDUSA']['Vars'] 		= ['FER',]	
			av['iron']['depthLevels'] 		= ['',]
			av['iron']['MEDUSA']['grid']		= Grid		

		if doO2:
			av['oxygen']['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
			av['oxygen']['Data']['Vars'] 	= ['o_an',] 
			av['oxygen']['MEDUSA']['Vars'] 	= ['OXY',]
			av['oxygen']['MEDUSA']['File']	= MEDUSAFolder+jobID+"_"+year+"_OXY.nc"
			av['oxygen']['MEDUSA']['grid']	= Grid
			av['oxygen']['depthLevels'] 	= depthLevels
		
		if doPCO2:
			av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
			av['pCO2']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"
			av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
			av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]	
			av['pCO2']['depthLevels'] 	= ['',]
			av['pCO2']['MEDUSA']['grid']	= Grid				


	#	if doIntPP:
	#		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
	#		av['intpp']['Data']['Vars'] 	= ['PPint',]
	#		#av['intpp']['ERSEM']['Vars'] 	= ['netPP',] This field is net, not integrated.
	#		#av['intpp']['ERSEM']['File']	= ERSEMFolder+'_ERSEMMisc.nc' This file will need to be remade.
	#		av['intpp']['ERSEM']['grid']	= Grid
	#		av['intpp']['depthLevels'] 		= ['',]
		
		
		if doSal:
			av['salinity']['Data']['File'] 		= WOAFolder+'salinity_monthly_1deg.nc'	
			av['salinity']['Data']['Vars'] 		= ['s_an',]
			av['salinity']['NEMO']['File'] 		= NEMOFolder+jobID+"_"+year+'_SAL.nc'	
			av['salinity']['NEMO']['Vars'] 		= ['vosaline',]
			av['salinity']['NEMO']['grid'] 		= Grid
			av['salinity']['depthLevels'] 		= depthLevels	 

		
		if doTemp:
			av['temperature']['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
			av['temperature']['Data']['Vars'] 	= ['t_an',]	
			av['temperature']['NEMO']['File'] 	= NEMOFolder+jobID+"_"+year+'_TEMP.nc'	
			av['temperature']['NEMO']['Vars'] 	= ['votemper',]
			av['temperature']['NEMO']['grid'] 	= Grid	
			av['temperature']['depthLevels'] 	= depthLevels	

						   
		if doMLD:	
			av['mld']['Data']['File'] 		= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
			av['mld']['Data']['Vars'] 		= ['mld','mask',]
			av['mld']['NEMO']['File'] 		= NEMOFolder+jobID+"_"+year+'_MLD.nc'			
			av['mld']['NEMO']['Vars'] 		= ['somxl010',]	
			av['mld']['NEMO']['grid'] 		= Grid
			av['mld']['depthLevels'] 		= ['',]

		
				
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
			 plotOceans		=0,#True,
			 plotHemispheres	=0,# True
			 plotSeasons		=0,# True
			 plotOceanSeasons	=0,# True		 		 
			 plotOceanMonths   	=0,#True	
			 plotHemispheresMonths  =0,#True,
		)

		for model in models:
			workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/"+model+'-'+jobID+'-'+year)
			imageFolder 	= folder('images/'+model+'-'+jobID)

	
			shelvesAV.extend(
		    		testsuite_p2p(
					model = model,
					jobID = jobID,
					year  = year,
					av = av,
					plottingSlices= plottingSlices,
					workingDir = workingDir,
					imageFolder= imageFolder
			 	)
			)
		
	
	#####
	# InterAnnual analysis
	allshelves = listShelvesContents(shelvesAV)
	
	for model in allshelves.models:
		for name in allshelves.names:
		    for depthLevel in allshelves.depthLevels:
		    	print 'InterAnnual analysis:',model,name,depthLevel
			#####
			# Each year is a different colour.
			yearslyShelves = {}
			for y in years:
				yearslyShelves[y] = reducesShelves(shelvesAV,models=[model],names=[name,],years=[y,],depthLevels=[depthLevel,],sliceslist=months,)
	
			filenamebase = folder('images/'+model+'-'+jobID+'/Patterns/Interannual/'+name+depthLevel)+model+jobID+name+'_'+years[0]+'-'+years[-1]+depthLevel
			makePatternStatsPlots(	
						yearslyShelves, # {legend, shelves}
						model+' '+name,				# title
						months,					# xkeysLabels
						filenamebase,				# filename base
						grid	= Grid,
						)	
					
			#####			
			# One long line showing multiple years - monthly. (looks best with up to 4 years.)
			longMonths, longShelves = [],[]
			for y in years:
				shelves = yearslyShelves[y]
				longMonths.extend([(mn,y) for mn in months])
				if len(shelves)==12:	longShelves.extend(shelves)
				else: 
					print "Shelves are not 12 long:",len(shelves)
					continue
					assert False
			if len(longMonths) != len(longShelves):
				print "len(longMonths) != len(longShelves):",len(longMonths),' != ',len(longShelves)
				continue
				assert False
		
			filenamebase = folder('images/'+model+'-'+jobID+'/Patterns/Interannual/'+name)+model+jobID+name+'_'+years[0]+'-'+years[-1]+depthLevel+'_longtimeseries'
			print "makePatternStatsPlots:",{name:longShelves}, model+' '+name,	longMonths,filenamebase,Grid
			makePatternStatsPlots(	
						{name:longShelves}, 		# {legend, shelves}
						model+' '+name+' ' +depthLevel,			# title
						longMonths,			# xkeysLabels
						filenamebase,			# filename base
						grid	= Grid,
						)		

			#####			
			# One long line showing multiple years - annual.
			for sl in ['All','Standard']:
				yearslyShelves = reducesShelves(shelvesAV,models=[model],names=[name,],years=years,depthLevels=[depthLevel,],sliceslist=[sl,],)	

				if len(yearslyShelves) != len(years):
					print "len(yearslyShelves) != len(years):",len(yearslyShelves),' != ',len(years)
					assert False
		
				filenamebase = folder('images/'+model+'-'+jobID+'/Patterns/Interannual/'+name)+model+jobID+name+'_'+years[0]+'-'+years[-1]+depthLevel+'_annual'+sl
				print "makePatternStatsPlots - Annual :",sl,{name:sorted(yearslyShelves)},years, model+' '+name,filenamebase,Grid
				makePatternStatsPlots(	
							{name:sorted(yearslyShelves)}, 		# {legend, shelves}
							model+' '+name+' ' +depthLevel,		# title
							sorted(years),				# xkeysLabels
							filenamebase,				# filename base
							grid	= Grid,
							)	
			
if __name__=="__main__":
	analysis()	
	
	
			

