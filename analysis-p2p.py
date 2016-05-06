#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
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
from sys import argv,exit
from os.path import exists
from calendar import month_name
from socket import gethostname
from glob import glob

#Specific local code:
import UKESMpython as ukp
from p2p import makePatternStatsPlots, testsuite_p2p
from p2p.patternAnalyses import InterAnnualPatterns,BGCvsPhysics
from pftnames import months
from p2p.shelveToDictionary import shelveToDictionary
#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#	

# 	from 




def analysis_jasmin(
		models	= ['NEMO','MEDUSA',],
		jobID 	= 'xkrus',
		years 	= ['2077'], #'2075','2076',
		modelGrid = 'ORCA1',
		annual 	= False,
		noPlots = False,
		analysisSuite='bio',
		):
	
	# DMS model:
	#model= 'MEDUSA'
	#jobID = 'xkrum'
	#year = 'clim'		
	#modelGrid = 'ORCA1'
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM_postProcessed/MEDUSA/outNetCDF/"+jobID+'-' + year+'/'

	# ORCA1:
	#models= ['MEDUSA','NEMO']
	#jobID = 'xkrus'
	#years = ['2077-annual'] #'2075','2076',
	#modelGrid = 'ORCA1'
	#MEDUSAFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	#NEMOFolder_pref= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/xkrus_postProc/"
	#annual = True


	
	# ORCA025:
	#model= 'MEDUSA'
	#jobID = 'xjwki'
	#year = '1979'		
	#modelGrid = 'ORCA025'	
	#MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/"+jobID+"_postProc/"+year+"/"
	
	
	
	#####
	# Which analysis to run
	
	if analysisSuite.lower() in ['default',]:
		doCHL 		= True
		doCHL_CCI 	= 0				
		doDiatoms	= 0#True
		doMicrozoo	= 0#True
		doMesozoo	= 0#True
		doN		= True
		doSi		= True	
		doFe		= True		
		doPCO2		= 0#True
		doIntPP		= 0#True
		
		doO2		= True	
		doSal		= True
		doTemp		= True
		doMLD		= True
	elif analysisSuite.lower() in ['debug',]:
		doCHL 		= True
		doCHL_CCI 	= 0				
		doDiatoms	= 0#True
		doMicrozoo	= 0#True
		doMesozoo	= 0#True		
		doN		= 0#True
		doSi		= 0#True	
		doFe		= 0#True		
		doPCO2		= 0#True
		doIntPP		= 0#True
		doO2		= 0#True	
		doSal		= 0#True
		doTemp		= 0#True
		doMLD		= 0#True
	elif analysisSuite.lower() in ['bio',]:
		doCHL 		= True
		doCHL_CCI 	= 0				
		doDiatoms	= True
		doMicrozoo	= True
		doMesozoo	= True		
		doN		= 0#True
		doSi		= 0#True	
		doFe		= 0#True		
		doPCO2		= 0#True
		doIntPP		= 0#True
		doO2		= 0#True	
		doSal		= 0#True
		doTemp		= 0#True
		doMLD		= 0#True			
	elif analysisSuite.lower() in ['annual',]:
		#####
		# Data sets with annual coverage available (instead o monthly)
		doCHL 		= 0
		doCHL_CCI 	= True		
		doDiatoms	= 0
		doMicrozoo	= 0
		doMesozoo	= 0		
		doN		= True
		doSi		= True	
		doFe		= 0#True		
		doPCO2		= 0#True
		doIntPP		= 0#True
		doO2		= True	
		doSal		= True
		doTemp		= True
		doMLD		= True	
	else:		
		doCHL 		= 0#True
		doCHL_CCI 	= 0				
		doDiatoms	= 0#True
		doMicrozoo	= 0#True
		doMesozoo	= 0#True
		doN		= True
		doSi		= 0#True	
		doFe		= 0#True
		doPCO2		= 0#True
		doIntPP		= 0#True
		doO2		= 0#True	
		doSal		= 0#True
		doTemp		= 0#True
		doMLD		= 0#True

	
	
	#####
	# What depth level to investigate, in the case of big 3D files (T,Sal, N,P,Si, etc)	
	#if annual:	
		print "WARNING: annual data not yet tested."
	#	assert False
	#if annual:depthLevels 	= ['',]
	#else:		depthLevels 	= ['Transect','Surface','100m','200m','500m',]


	#####
	# Location of data files.
	if gethostname().find('pmpc')>-1:	
		print "analysis-JASMIN.py:\tBeing run at PML on ",gethostname()
		
		if annual:	
			#####
			# No need to stitch together multiple months into one file:
			WOAFolder 		= "/data/euryale7/scratch/ledm/WOA/annual/"
			MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
			NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"			
		else:		
			WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
			MEDUSAFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"_postProc/"
			NEMOFolder_pref		= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"_postProc/"			
		#MAREDATFolder 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"
		#GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
		#TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
		#MLDFolder	= "/data/euryale7/scratch/ledm/IFREMER-MLD/"
		
		ObsFolder = "/data/euryale7/backup/ledm/Observations/"
		MAREDATFolder 	= ObsFolder+"/MAREDAT/MAREDAT/"
		GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
		TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
		MLDFolder	= ObsFolder+"/IFREMER-MLD/"
		iMarNetFolder	= ObsFolder+"/LestersReportData/"
		GlodapDir	= ObsFolder+"/GLODAP/"
		GLODAPv2Dir	= ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
		OSUDir		= ObsFolder+"OSU/"
		CCIDir		= ObsFolder+"CCI/"
				
		if jobID in ["xkrus",]:
			# Old school ORCA1 grid
			orcaGridfn 	='data/mesh_mask_ORCA1_75.nc'
		else:
			# New eORCA1 grid		
			orcaGridfn 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'
				
		workDir		= "/data/euryale7/scratch/ledm/ukesm_postProcessed/"
		imgDir		= ukp.folder('images')
		
	if gethostname().find('ceda.ac.uk')>-1:
		print "analysis-JASMIN.py:\tBeing run at CEDA on ",gethostname()
			
		esmvalFolder = "/group_workspaces/jasmin/esmeval/example_data/bgc/"
		
		#####
		# Location of model files.	
		MEDUSAFolder_pref	= ukp.folder(esmvalFolder+"MEDUSA/")
		NEMOFolder_pref		= ukp.folder(esmvalFolder+"MEDUSA/")
		
		#####
		# Location of data files.
		if annual:	WOAFolder 	= ukp.folder(esmvalFolder+"WOA/annual")
		else:		WOAFolder 	= ukp.folder(esmvalFolder+"WOA/")
		MAREDATFolder 	= ukp.folder(esmvalFolder+"/MAREDAT/")
		WOAFolder 	= ukp.folder(esmvalFolder+"WOA/")
		GEOTRACESFolder = ukp.folder(esmvalFolder+"GEOTRACES/GEOTRACES_PostProccessed/")
		TakahashiFolder = ukp.folder(esmvalFolder+"Takahashi2009_pCO2/")
		MLDFolder  	= ukp.folder(esmvalFolder+"IFREMER-MLD/")
	
		# Directory for output files:
		workDir 	= ukp.folder(esmvalFolder+"ukesm_postProcessed/")
		imgDir		= ukp.folder('images')		

		if jobID in ["xkrus",]:
			# Old school ORCA1 grid
			orcaGridfn 	='/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_ORCA1_75.nc'
		else:
			# New eORCA1 grid		
			orcaGridfn 	= '/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'
			
						
	#####
	# Set which spatial and temporal limitations to plot.
	transects 	= ['AtlanticTransect', 'PacificTransect',]
	justAll		= ['All',]				# All is not a slice, it has no cut on location, time, or depth.
	AllStandard	= ['All','Standard','ignoreInlandSeas']	
	HighLatWinter	= ['All','HighLatWinter',]
						
	medusaCoords 	= {'t':'index_t', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	maredatCoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	woaCoords 	= {'t':'index_t', 'z':'depth',  'lat': 'lat', 	   'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}		
	cciCoords	= {'t':'index_t', 'z':'index_z','lat': 'lat',      'lon': 'lon',       'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	shelvesAV = []	
	
	for year in years:		
		#####
		# Location of model files.
		if annual:
			MEDUSAFolder 	= MEDUSAFolder_pref+jobID+"/"
			NEMOFolder 	= NEMOFolder_pref+jobID+"/"			
		else:			
			MEDUSAFolder 	= MEDUSAFolder_pref+year+'/'
			NEMOFolder  	= NEMOFolder_pref+year+'/'		
		

		#####
		# AutoVivification is a form of nested dictionary.
		# We use AutoVivification here to determine which files to analyse and which fields in those files.
		# depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
		av = ukp.AutoVivification()
		if doCHL:
			av['chl']['Data']['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
			if modelGrid == 'ORCA1':	av['chl']['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_CHL.nc"
			if modelGrid == 'ORCA025':	av['chl']['MEDUSA']['File']	= MEDUSAFolder+"xjwki_1979_CH.nc"
			
			av['chl']['Data']['coords'] 		= maredatCoords
			av['chl']['MEDUSA']['coords']		= medusaCoords
			
			av['chl']['MEDUSA']['details']		= {'name': 'CHL', 'vars':['CHL',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			
			av['chl']['Data']['details']		= {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}			

			av['chl']['Data']['source'] 		= 'MAREDAT'
			av['chl']['MEDUSA']['source']		= 'MEDUSA'

			av['chl']['depthLevels'] 		= ['',]
			av['chl']['MEDUSA']['grid']		= modelGrid		
			av['chl']['plottingSlices'] 		= AllStandard

		if doCHL_CCI:						
			name = 'Chlorophyll_cci'
			if annual:
				av[name]['Data']['File'] 	= CCIDir+"ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-annual-fv2.0.nc"	
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]					
			else:
				av[name]['Data']['File'] 	= CCIDir+'ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-all-fv2.0.nc'
				av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_CHL.nc"
			
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['',]
			if annual:	av[name]['plottingSlices'] 	= AllStandard
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= cciCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'CCI'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']		= {'name': name, 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}		
			av[name]['Data']['details']		= {'name': name, 'vars':['chlor_a',], 'convert':  ukp.NoChange,'units':'mg C/m^3'}			
		
		if doDiatoms:
			if annual: 
				print "No diatoms iron file",
				assert 0		
			av['diatoms']['Data']['File'] 		= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
			av['diatoms']['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_PHD.nc"
			
			av['diatoms']['depthLevels'] 		= ['',]	
			av['diatoms']['MEDUSA']['grid']		= modelGrid						
			av['diatoms']['plottingSlices']		= AllStandard

			av['diatoms']['Data']['coords'] 	= maredatCoords
			av['diatoms']['MEDUSA']['coords']	= medusaCoords

			av['diatoms']['MEDUSA']['details']	= {'name': 'diatoms', 'vars':['PHD',],     'convert': ukp.N2Biomass,'units': 'mg C/m^3'}						
			av['diatoms']['Data']['details']	= {'name': 'diatoms', 'vars':['BIOMASS',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			

			av['diatoms']['Data']['source'] 	= 'MAREDAT'
			av['diatoms']['MEDUSA']['source']	= 'MEDUSA'

			
		if doMicrozoo:
			if annual: 
				print "No microzoo iron file",
				assert 0		
			av['microzoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
			av['microzoo']['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_ZMI.nc"	
			
			av['microzoo']['MEDUSA']['grid']	= modelGrid		
			av['microzoo']['depthLevels'] 		= ['',]	
			av['microzoo']['plottingSlices'] 	= AllStandard

			av['microzoo']['Data']['coords'] 	= maredatCoords
			av['microzoo']['MEDUSA']['coords']	= medusaCoords
			
			av['microzoo']['MEDUSA']['details']	= {'name': 'microzoo', 'vars':['ZMI',],     'convert': ukp.N2Biomass,'units': 'mg C/m^3'}			
			av['microzoo']['Data']['details']	= {'name': 'microzoo', 'vars':['BIOMASS',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			

			av['microzoo']['Data']['source'] 	= 'MAREDAT'
			av['microzoo']['MEDUSA']['source']	= 'MEDUSA'
			
			
		if doMesozoo:
			if annual: 
				print "No mesozoo iron file",
				assert 0
						
			av['mesozoo']['Data']['File'] 		= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
			av['mesozoo']['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_ZME.nc"	

			av['mesozoo']['MEDUSA']['grid']		= modelGrid		
			av['mesozoo']['depthLevels'] 		= ['',]
			av['mesozoo']['plottingSlices'] 	= AllStandard

			av['mesozoo']['Data']['coords'] 	= maredatCoords
			av['mesozoo']['MEDUSA']['coords']	= medusaCoords
			
			av['mesozoo']['MEDUSA']['details']	= {'name': 'mesozoo', 'vars':['ZME',],     'convert': ukp.N2Biomass,'units': 'mg C/m^3'}			
			av['mesozoo']['Data']['details']	= {'name': 'mesozoo', 'vars':['BIOMASS',], 'convert': ukp.NoChange,'units':'mg C/m^3'}			

			av['mesozoo']['Data']['source'] 	= 'MAREDAT'
			av['mesozoo']['MEDUSA']['source']	= 'MEDUSA'
			
			
		if doN:
			name = 'nitrate'		
			if annual:	
				av[name]['Data']['File'] 	= WOAFolder+'woa13_all_n00_01.nc'
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]	
			else:		
				av[name]['Data']['File'] 	= WOAFolder+'nitrate_monthly_1deg.nc'	
				if modelGrid == 'ORCA1':	av[name]['MEDUSA']['File'] = MEDUSAFolder+jobID+'_' + year+"_DIN.nc"	
				if modelGrid == 'ORCA025':	av[name]['MEDUSA']['File'] = MEDUSAFolder+jobID+'_'+ year+"_DIN.nc"							
			
			av[name]['MEDUSA']['grid']	= modelGrid		
			av[name]['depthLevels'] 	= ['Surface','Transect','PTransect','SOTransect',]
			if annual:	av[name]['plottingSlices'] 	= AllStandard
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['DIN',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['n_an',], 'convert': ukp.NoChange,}	# no units?

						
		if doSi:
			name = 'silicate'
			if annual:
				av[name]['Data']['File'] 	= WOAFolder+'woa13_all_i00_01.nc'
				av[name]['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]					
			else:
				av[name]['Data']['File'] 	= WOAFolder+'silicate_monthly_1deg.nc'	
				av[name]['MEDUSA']['File'] 	= MEDUSAFolder+jobID+'_' + year+"_SIL.nc"
			
			av[name]['MEDUSA']['grid']		= modelGrid		
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']
			if annual:	av[name]['plottingSlices'] 	= AllStandard
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['MEDUSA']['coords']	= medusaCoords
			
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['MEDUSA']['source']	= 'MEDUSA'			
	
			av[name]['MEDUSA']['details']	= {'name': name, 'vars':['SIL',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['i_an',], 'convert': ukp.NoChange,}	# no units?
			
						
		if doFe:
			if annual: 
				print "No annual iron file",
				assert 0
				
			av['iron']['Data']['File'] 		= GEOTRACESFolder+"Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
			av['iron']['MEDUSA']['File'] 		= MEDUSAFolder+jobID+'_' + year+"_FER.nc"	
			
			av['iron']['depthLevels'] 		= ['',]
			av['iron']['MEDUSA']['grid']		= modelGrid		
			av['iron']['plottingSlices']		= justAll
			
			av['iron']['Data']['coords'] 	= {'t': 'MONTH','z':'DEPTH','lat':'Latitude','lon':'Longitude','cal':'standard','tdict': ukp.tdicts['OneToZero']}
			av['iron']['MEDUSA']['coords']	= medusaCoords
	
			
			av['iron']['Data']['source'] 	= 'GEOTRACES'
			av['iron']['MEDUSA']['source']	= 'MEDUSA'			
	
			av['iron']['MEDUSA']['details']	= {'name': 'iron', 'vars':['FER',], 'convert': ukp.mul1000,'units':'umol F/m^3'}			
			av['iron']['Data']['details']	= {'name': 'iron', 'vars':['Fe_D_CONC_BOTTLE',], 'convert': ukp.NoChange,}	# no units?

			
		if doO2:
			if annual:
				av['oxygen']['MEDUSA']['File'] 	= sorted(glob(MEDUSAFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_ptrc_T.nc"))[0]					
				av['oxygen']['Data']['File'] 	=  WOAFolder+'woa13_all_o00_01.nc'
			else:	
				av['oxygen']['Data']['File'] 	=  WOAFolder+'oxygen-woa13.nc'
				av['oxygen']['MEDUSA']['File']	= MEDUSAFolder+jobID+"_"+year+"_OXY.nc"
			
			av['oxygen']['MEDUSA']['grid']		= modelGrid		
			av['oxygen']['depthLevels'] 		= ['Surface','Transect','PTransect']
			if annual:	av[name]['plottingSlices'] 	= AllStandard
			else:		av[name]['plottingSlices'] 	= HighLatWinter
			
			av['oxygen']['Data']['coords'] 		= woaCoords
			av['oxygen']['MEDUSA']['coords']	= medusaCoords
			
			av['oxygen']['Data']['source'] 		= 'WOA'
			av['oxygen']['MEDUSA']['source']	= 'MEDUSA'			
	
			av['oxygen']['MEDUSA']['details']	= {'name': 'oxygen', 'vars':['OXY',], 'convert': ukp.NoChange,}			
			av['oxygen']['Data']['details']		= {'name': 'oxygen', 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol/m^3'}
			
					
#		#if doPCO2:
#			av['pCO2']['Data']['File'] 	= TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"	
#			av['pCO2']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"
#			av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
#			av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]	
#			av['pCO2']['depthLevels'] 	= ['',]
#			av['pCO2']['MEDUSA']['grid']	= modelGrid				
#			#av['pCO2']['plottingSlices'] 	= []
#
#	#	if doIntPP:
#	#		av['intpp']['Data']['File'] 	=  LesterFolder+'PPint_1deg.nc'
#	#		av['intpp']['Data']['Vars'] 	= ['PPint',]
#	#		#av['intpp']['ERSEM']['Vars'] 	= ['netPP',] This field is net, not integrated.
#	#		#av['intpp']['ERSEM']['File']	= ERSEMFolder+'_ERSEMMisc.nc' This file will need to be remade.
#	#		av['intpp']['ERSEM']['grid']	= modelGrid
#	#		av['intpp']['depthLevels'] 		= ['',]
		
		
		if doSal:
			name = 'salinity'
			if annual:
				av[name]['NEMO']['File'] 	= sorted(glob(NEMOFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_grid_T.nc"))[0]
				av[name]['Data']['File'] 	= WOAFolder+'woa13_decav_s00_01v2.nc'	
			else:	
				av[name]['Data']['File'] 	= WOAFolder+'salinity_monthly_1deg.nc'	
				av[name]['NEMO']['File'] 	= NEMOFolder+jobID+"_"+year+'_SAL.nc'	

			av[name]['NEMO']['grid'] 		= modelGrid
			av[name]['depthLevels'] 		= ['Surface','Transect','PTransect']	 
			av[name]['plottingSlices'] 	= AllStandard
			
			av[name]['Data']['coords'] 	= woaCoords
			av[name]['NEMO']['coords']	= medusaCoords
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['NEMO']['source']	= 'NEMO'			
	
			av[name]['NEMO']['details']	= {'name': name, 'vars':['vosaline',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['s_an',], 'convert': ukp.NoChange,}	# no units?
			
					
		if doTemp:
			name = 'temperature'
			if annual:
				av[name]['NEMO']['File'] 	= sorted(glob(NEMOFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_grid_T.nc"))[0]				
				av[name]['Data']['File'] 	= WOAFolder+'woa13_decav_t00_01v2.nc'	
			else:	
				av[name]['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'	
				av[name]['NEMO']['File'] 	= NEMOFolder+jobID+"_"+year+'_TEMP.nc'	

			av[name]['NEMO']['grid'] 	= modelGrid	
			av[name]['depthLevels'] 	= ['Surface','Transect','PTransect']	
			av[name]['plottingSlices'] 	= AllStandard

			av[name]['Data']['coords'] 	= woaCoords
			av[name]['NEMO']['coords']	= medusaCoords
			av[name]['Data']['source'] 	= 'WOA'
			av[name]['NEMO']['source']	= 'NEMO'			
	
			av[name]['NEMO']['details']	= {'name': name, 'vars':['votemper',], 'convert': ukp.NoChange,}			
			av[name]['Data']['details']	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,}	# no units?
			
						   
		if doMLD:
			name = 'mld'		
			if annual:	
				av[name]['NEMO']['File'] 	= sorted(glob(NEMOFolder_pref+jobID+"/"+jobID+"o_1y_*1201_"+year+"1130_grid_T.nc"))[0]							
				av[name]['Data']['File'] 		= MLDFolder+"mld_DT02_c1m_reg2.0-annual.nc"
			else:	
				av[name]['Data']['File'] 		= MLDFolder+"mld_DT02_c1m_reg2.0.nc"
				av[name]['NEMO']['File'] 		= NEMOFolder+jobID+"_"+year+'_MLD.nc'	
					
			av[name]['NEMO']['grid'] 		= modelGrid
			av[name]['depthLevels'] 		= ['',]
			av[name]['plottingSlices'] 		= AllStandard

			av[name]['Data']['coords'] 	= {'t':'index_t', 'z':'index_z','lat':'lat','lon':'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
			av[name]['NEMO']['coords']	= medusaCoords
			av[name]['Data']['source'] 	= 'IFREMER'
			av[name]['NEMO']['source']	= 'NEMO'			
	
			av[name]['NEMO']['details']	= {'name': name, 'vars':['somxl010',], 'convert': ukp.NoChange,'units':'m'}			
			av[name]['Data']['details']	= {'name': name, 'vars':['mld','mask',], 'convert': ukp.applymask,'units':'m'}	# no units?
			
		
		for model in models:
			workingDir 	= ukp.folder(workDir+model+'-'+jobID+'-'+year)
			imageFolder 	= ukp.folder(imgDir+'/'+jobID)

			shelvesAV.extend(
		    		testsuite_p2p(
					model 		= model,
					jobID 		= jobID,
					year  		= year,
					av 		= av,
					plottingSlices	= [],	# set this so that testsuite_p2p reads the slice list from the av.
					workingDir 	= workingDir,
					imageFolder	= imageFolder,
					noPlots		= noPlots,	# turns off plot making to save space and compute time.
					gridFile	= orcaGridfn,	# enforces custom gridfile.
					annual		= annual,
			 	)
			)
	
		
	#BGCvsPhysics(shelvesAV, jobID, modelGrid )
	#if len(years)>1: InterAnnualPatterns(shelvesAV, jobID, years,modelGrid)	# plots interannual comparison and time series.
#	def outPutForJASMIN(shelvesAV):
#	outdict = shelveToDictionary(shelvesAV)
	
		
if __name__=="__main__":
	#analysis_jasmin()	
	analysis_jasmin(models	= ['NEMO','MEDUSA',],
		jobID 	= 'u-ab749',
		years 	= ['2007'], #'2075','2076',
		modelGrid = 'eORCA1',
		annual 	= True,
		noPlots = False,
		analysisSuite='annual',)
	
			

