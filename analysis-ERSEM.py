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
#from os.path import exists
#from calendar import month_name

#Specific local code:
from UKESMpython import populateSlicesList, AutoVivification, folder, reducesShelves,mnStr
	#,getFileList, , AutoVivToYaml,YamlToDict, slicesDict,reducesShelves
from p2p import makePatternStatsPlots
from pftnames import months#MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names,getmt

from testsuite_p2p import testsuite_p2p
#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	
#	This code is the 

# 	from 
def analysis():
	model= 'ERSEM'
	jobID = 'xhonc'
	years = ['1996','1997','1998','1999','2000']		# arbitrary year chosen for 
	grid = 'ORCA1'
	#####
	# Which analysis to run
	doCHL 		= True
	doMAREDAT 	= True
	doN		= 0#True
	doP		= 0#True
	doSi		= 0#True	
	doO2		= 0#True	
	doFe		= True
	doPCO2		= 0#True
	doIntPP		= 0#True


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

	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)
	depthLevels 	= ['Surface','100m','200m','500m',]
		
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
		ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobID+'/'+year+'/'+jobID+'_'+year	
	
	
		#####
		# AutoVivification is a form of nested dictionary.
		# We use AutoVivification here to determine which files to analyse and which fields in those files.
		# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
		av = AutoVivification()
		if doCHL:
			av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
			av['chl']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMMisc.nc'	
			av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
			av['chl']['ERSEM']['Vars'] 		= ['chl',]
			av['chl']['depthLevels'] 		= ['',]
			av['chl']['ERSEM']['grid']		= grid		
						
		if doMAREDAT:
			av['diatoms']['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
			av['diatoms']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMphytoBm.nc'				
			av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
			av['diatoms']['ERSEM']['Vars'] 		= ['P1c',]
			av['diatoms']['depthLevels'] 		= ['',]	
			av['diatoms']['ERSEM']['grid']		= grid						
	
			av['bac']['Data']['File'] 		= MAREDATFolder+"MarEDat20120214Bacteria.nc"	
			av['bac']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMbac.nc'			
			av['bac']['Data']['Vars'] 		= ['BIOMASS',]
			av['bac']['ERSEM']['Vars'] 		= ['B1c',]
			av['bac']['depthLevels'] 		= ['',]
			av['bac']['ERSEM']['grid']		= grid						
	
			av['picophyto']['Data']['File'] 	= MAREDATFolder+"MarEDat20111206Picophytoplankton.nc"	
			av['picophyto']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMphytoBm.nc'			
			av['picophyto']['Data']['Vars'] 	= ['BIOMASS',]
			av['picophyto']['ERSEM']['Vars'] 	= ['P3c',]
			av['picophyto']['depthLevels'] 		= ['',]	
			av['picophyto']['ERSEM']['grid']		= grid						
		
			av['microzoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
			av['microzoo']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMzoo.nc'			
			av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
			av['microzoo']['ERSEM']['Vars'] 	= ['Z5c',]
			av['microzoo']['depthLevels'] 		= ['',]	
			av['microzoo']['ERSEM']['grid']		= grid						
	
			av['mesozoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
			av['mesozoo']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMzoo.nc'			
			av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
			av['mesozoo']['ERSEM']['Vars'] 		= ['Z4c',]
			av['mesozoo']['depthLevels'] 		= ['',]
			av['mesozoo']['ERSEM']['grid']		= grid						

		if doN:
			av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
			av['nitrate']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'	
			av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
			av['nitrate']['ERSEM']['Vars'] 		= ['N3n','N4n',]
			av['nitrate']['depthLevels'] 		= depthLevels
			av['nitrate']['ERSEM']['grid']		= grid				
			
		if doSi:
			av['silicate']['Data']['File'] 		= WOAFolder+'silicate_monthly_1deg.nc'	
			av['silicate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
			av['silicate']['Data']['Vars'] 		= ['i_an',] 		#l+'_mn',
			av['silicate']['ERSEM']['Vars'] 	= ['N5s',]
			av['silicate']['depthLevels'] 		= depthLevels
			av['silicate']['ERSEM']['grid']		= grid	
			
		if doP:						
			av['phosphate']['Data']['File'] 	= WOAFolder+'phosphate_monthly_1deg.nc'	
			av['phosphate']['Data']['Vars'] 	= ['p_an',] 		#l+'_mn',
			av['phosphate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
			av['phosphate']['ERSEM']['Vars'] 	= ['N1p',]
			av['phosphate']['ERSEM']['grid']	= grid	
			av['phosphate']['depthLevels'] 		= depthLevels
			
		if doFe:								
			av['iron']['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
			av['iron']['Data']['Vars'] 		= ['Fe_D_CONC_BOTTLE',]
			av['iron']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'			
			av['iron']['ERSEM']['Vars'] 		= ['N7f',]
			av['iron']['ERSEM']['grid']		= grid				
			av['iron']['depthLevels'] 		= ['',]
			
		if doO2:
			av['oxygen']['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
			av['oxygen']['Data']['Vars'] 	= ['o_an',] 
			av['oxygen']['ERSEM']['Vars'] 	= ['O2o',]
			av['oxygen']['ERSEM']['File']	= ERSEMFolder+'_ERSEMO2.nc'
			av['oxygen']['ERSEM']['grid']	= grid
			av['oxygen']['depthLevels'] 	= depthLevels
		
		if doPCO2:
			av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
			av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
			av['pCO2']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMMisc.nc'	
			av['pCO2']['ERSEM']['Vars'] 	= ['pCO2w',]
			av['pCO2']['ERSEM']['grid'] 	= grid		
			av['pCO2']['depthLevels'] 	= ['',]

		if doIntPP:
			av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
			av['intpp']['Data']['Vars'] 	= ['PPint',]
			#av['intpp']['ERSEM']['Vars'] 	= ['netPP',] This field is net, not integrated.
			#av['intpp']['ERSEM']['File']	= ERSEMFolder+'_ERSEMMisc.nc' This file will need to be remade.
			av['intpp']['ERSEM']['grid']	= grid
			av['intpp']['depthLevels'] 	= ['',]
		


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
			imageFolder= imageFolder)
		    )
	
	#####
	# InterAnnual analysis
	for name in av.keys():
	   for depthLevel in av[name]['depthLevels']:
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
					grid	= grid,
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
				assert False
		if len(longMonths) != len(longShelves):
			print "len(longMonths) != len(longShelves):",len(longMonths),' != ',len(longShelves)
			assert False
		
		filenamebase = folder('images/'+model+'-'+jobID+'/Patterns/Interannual/'+name)+model+jobID+name+'_'+years[0]+'-'+years[-1]+depthLevel+'_longtimeseries'
		print "makePatternStatsPlots:",{name:longShelves}, model+' '+name,	longMonths,filenamebase,grid
		makePatternStatsPlots(	
					{name:longShelves}, 		# {legend, shelves}
					model+' '+name+' ' +depthLevel,			# title
					longMonths,			# xkeysLabels
					filenamebase,			# filename base
					grid	= grid,
					)		

		#####			
		# One long line showing multiple years - annual.
		for sl in ['All','Standard']:
			yearslyShelves = reducesShelves(shelvesAV,models=[model],names=[name,],years=years,depthLevels=[depthLevel,],sliceslist=[sl,],)	

			if len(yearslyShelves) != len(years):
				print "len(yearslyShelves) != len(years):",len(yearslyShelves),' != ',len(years)
				assert False
		
			filenamebase = folder('images/'+model+'-'+jobID+'/Patterns/Interannual/'+name)+model+jobID+name+'_'+years[0]+'-'+years[-1]+depthLevel+'_annual'+sl
			print "makePatternStatsPlots - Annual :",sl,{name:sorted(yearslyShelves)},years, model+' '+name,filenamebase,grid
			makePatternStatsPlots(	
						{name:sorted(yearslyShelves)}, 		# {legend, shelves}
						model+' '+name+' ' +depthLevel,		# title
						sorted(years),				# xkeysLabels
						filenamebase,				# filename base
						grid	= grid,
						)	
								
if __name__=="__main__":
	analysis()	
	print "The end."
	
	
	
	
	
