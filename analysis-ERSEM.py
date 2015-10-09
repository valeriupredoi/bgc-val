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
from UKESMpython import populateSlicesList, AutoVivification, folder#,getFileList, , AutoVivToYaml,YamlToDict, slicesDict,reducesShelves
#from p2p import matchDataAndModel,makePlots,makeTargets, csvFromShelves, makePatternStatsPlots
#from pftnames import MaredatTypes,WOATypes,Ocean_names,OceanMonth_names,months, Seasons,Hemispheres,HemispheresMonths, OceanSeason_names,getmt

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
	year = '1997'		# arbitrary year chosen for 
	
	#####
	# Which analysis to run
	doCHL 		= 0#True
	doMAREDAT 	= 0#True
	doN		= 0#True
	doP		= 0#True
	doSi		= 0#True	
	doFe		= 0#True
	doPCO2		= 0#True
	doIntPP		= 0#True
	doO2		= True	


	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
	LesterFolder 	= "/data/euryale7/scratch/ledm/LestersReportData/"			
	#####
	# Location of model files.	
	ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"+ jobID+'/'+year+'/'+jobID+'_'+year

	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)
	depthLevels 	= ['Surface','100m','200m','500m',]
	
	
	
	
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
		av['chl']['depthLevels'] 			= ['',]
		av['chl']['ERSEM']['grid']		= 'ORCA1'				


		
						
	if doMAREDAT:
		av['diatoms']['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
		av['diatoms']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMphytoBm.nc'				
		av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
		av['diatoms']['ERSEM']['Vars'] 		= ['P1c',]
		av['diatoms']['depthLevels'] 		= ['',]	
		av['diatoms']['ERSEM']['grid']		= 'ORCA1'						
	
		av['bac']['Data']['File'] 		= MAREDATFolder+"MarEDat20120214Bacteria.nc"	
		av['bac']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMbac.nc'			
		av['bac']['Data']['Vars'] 		= ['BIOMASS',]
		av['bac']['ERSEM']['Vars'] 		= ['B1c',]
		av['bac']['depthLevels'] 		= ['',]
		av['bac']['ERSEM']['grid']		= 'ORCA1'						
	
		av['picophyto']['Data']['File'] 	= MAREDATFolder+"MarEDat20111206Picophytoplankton.nc"	
		av['picophyto']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMphytoBm.nc'			
		av['picophyto']['Data']['Vars'] 	= ['BIOMASS',]
		av['picophyto']['ERSEM']['Vars'] 	= ['P3c',]
		av['picophyto']['depthLevels'] 		= ['',]	
		av['picophyto']['ERSEM']['grid']		= 'ORCA1'						
		
		av['microzoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
		av['microzoo']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMzoo.nc'			
		av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['microzoo']['ERSEM']['Vars'] 	= ['Z5c',]
		av['microzoo']['depthLevels'] 		= ['',]	
		av['microzoo']['ERSEM']['grid']		= 'ORCA1'						
	
		av['mesozoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
		av['mesozoo']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMzoo.nc'			
		av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
		av['mesozoo']['ERSEM']['Vars'] 		= ['Z4c',]
		av['mesozoo']['depthLevels'] 		= ['',]
		av['mesozoo']['ERSEM']['grid']		= 'ORCA1'						

	if doN:
		av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
		av['nitrate']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'	
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
		av['nitrate']['ERSEM']['Vars'] 		= ['N3n','N4n',]
		av['nitrate']['depthLevels'] 		= depthLevels
		av['nitrate']['ERSEM']['grid']		= 'ORCA1'				
			
	if doSi:
		av['silicate']['Data']['File'] 		= WOAFolder+'silicate_monthly_1deg.nc'	
		av['silicate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
		av['silicate']['Data']['Vars'] 		= ['i_an',] 		#l+'_mn',
		av['silicate']['ERSEM']['Vars'] 	= ['N5s',]
		av['silicate']['depthLevels'] 		= depthLevels
		av['silicate']['ERSEM']['grid']		= 'ORCA1'				
	if doP:						
		av['phosphate']['Data']['File'] 	= WOAFolder+'phosphate_monthly_1deg.nc'	
		av['phosphate']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMNuts.nc'	
		av['phosphate']['Data']['Vars'] 	= ['p_an',] 		#l+'_mn',
		av['phosphate']['ERSEM']['Vars'] 	= ['N1p',]
		av['phosphate']['depthLevels'] 		= depthLevels		
		av['phosphate']['ERSEM']['grid']	= 'ORCA1'				
	if doFe:								
		av['iron']['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		av['iron']['ERSEM']['File'] 		= ERSEMFolder+'_ERSEMNuts.nc'			
		av['iron']['Data']['Vars'] 		= ['Fe_D_CONC_BOTTLE',]
		av['iron']['ERSEM']['Vars'] 		= ['N7f',]
		av['iron']['depthLevels'] 		= ['',]
		av['iron']['ERSEM']['grid']		= 'ORCA1'				

	if doO2:
		av['oxygen']['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
		av['oxygen']['Data']['Vars'] 	= ['o_an',] 
		av['oxygen']['ERSEM']['Vars'] 	= ['O2o',]
		av['oxygen']['ERSEM']['File']	= ERSEMFolder+'_ERSEMO2.nc'
		av['oxygen']['ERSEM']['grid']	= 'ORCA1'
		av['oxygen']['depthLevels'] 	= depthLevels
		
	if doPCO2:
		av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
		av['pCO2']['ERSEM']['File'] 	= ERSEMFolder+'_ERSEMMisc.nc'	
		av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
		av['pCO2']['ERSEM']['Vars'] 	= ['pCO2w',]
		av['pCO2']['depthLevels'] 		= ['',]
		av['pCO2']['ERSEM']['grid'] 	= 'ORCA1'		


	if doIntPP:
		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
		av['intpp']['Data']['Vars'] 	= ['PPint',]
		#av['intpp']['ERSEM']['Vars'] 	= ['netPP',] This field is net, not integrated.
		#av['intpp']['ERSEM']['File']	= ERSEMFolder+'_ERSEMMisc.nc' This file will need to be remade.
		av['intpp']['ERSEM']['grid']	= 'ORCA1'
		av['intpp']['depthLevels'] 	= ['',]
		
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

	
	shelvesAV = testsuite_p2p(
		model = model,
		jobID = jobID,
		year  = year,
		av = av,
		plottingSlices= plottingSlices,
		workingDir = workingDir,
		imageFolder= imageFolder)
	
	
			
if __name__=="__main__":
	analysis()	
	
	
