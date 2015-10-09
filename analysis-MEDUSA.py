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
from UKESMpython import populateSlicesList, AutoVivification, folder
from testsuite_p2p import testsuite_p2p


#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	

# 	from 
def analysis():
	model= 'MEDUSA'
	jobID = 'xkrum'
	year = 'clim'		
	
	Grid = 'ORCA1'
	
	#####
	# Which analysis to run
	doCHL 		= True
	doMAREDAT 	= True
	doN		= True
	doSi		= True	
	doFe		= True		
	doPCO2		= True
	doIntPP		= 0#True
	doO2		= True	

	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)	
	depthLevels 	= ['Surface','200m','500m',]#'100m',		

	#####
	# Location of data files.
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
	TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
	LesterFolder 	= "/data/euryale7/scratch/ledm/LestersReportData/"			
	#####
	# Location of model files.	
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobID+'-' + year+'/'
	

	#####
	# AutoVivification is a form of nested dictionary.
	# We use AutoVivification here to determine which files to analyse and which fields in those files.
	# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
	av = AutoVivification()
	if doCHL:
		av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
		av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+jobID+'_' + year+"_MEDUSA_bio.nc"
		av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
		av['chl']['MEDUSA']['Vars'] 		= ['CHD','CHN']	
		av['chl']['depthLevels'] 			= ['',]
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
		av['mesozoo']['MEDUSA']['grid']	= Grid		
		av['mesozoo']['depthLevels'] 		= ['',]

	if doN:
		av['nitrate']['Data']['File'] 		= WOAFolder+'nitrate_monthly_1deg.nc'	
		av['nitrate']['Data']['Vars'] 		= ['n_an',] 		#l+'_mn',
		av['nitrate']['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_'+ year+"_MEDUSA_bio.nc"		
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
		av['iron']['depthLevels'] 			= ['',]
		av['iron']['MEDUSA']['grid']		= Grid		

#	if doO2:
#		av['oxygen']['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
#		av['oxygen']['Data']['Vars'] 	= ['o_an',] 
#		av['oxygen']['ERSEM']['Vars'] 	= ['o2c',]						
#		av['oxygen']['ERSEM']['File']	= ERSEMFolder+'_ERSEMO2.nc'
#		av['oxygen']['ERSEM']['grid']	= Grid
#		av['oxygen']['depthLevels'] 	= depthLevels
		
	if doPCO2:
		av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
		av['pCO2']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"
		av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
		av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]	
		av['pCO2']['depthLevels'] 		= ['',]
		av['pCO2']['MEDUSA']['grid']	= Grid				


#	if doIntPP:
#		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
#		av['intpp']['Data']['Vars'] 	= ['PPint',]
#		#av['intpp']['ERSEM']['Vars'] 	= ['netPP',] This field is net, not integrated.
#		#av['intpp']['ERSEM']['File']	= ERSEMFolder+'_ERSEMMisc.nc' This file will need to be remade.
#		av['intpp']['ERSEM']['grid']	= Grid
#		av['intpp']['depthLevels'] 		= ['',]
		
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
	
	
			

