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

from calendar import month_name
from UKESMpython import AutoVivification,AutoVivToYaml,folder,YamlToDict
#from itertools import product
from os.path import exists
import numpy as np



#####
#	

regions 	= ['Surface','200m','100m','500m','1000m','Transect','All','',]

MaredatTypes 	= ['chl','diatoms','bac','mesozoo','picophyto','microzoo']

Ocean_names	= ['SouthPacificOcean',  'ArcticOcean',  'AntarcticOcean',
			'NorthAtlanticOcean','SouthAtlanticOcean', 
			'NorthPacificOcean','IndianOcean',
			'EquatorialPacificOcean','EquatorialAtlanticOcean',]

IFREMERTypes 	= ['mld','mld_DT02','mld_DR003','mld_DReqDTm02', ]

WOATypes 	= ['silicate','nitrate','phosphate','salinity','temperature','oxygen']

CMIP5models = [ 'MEDUSA','ERSEM','BNU-ESM', 'IPSL-CM5A-LR', 'CESM1-BGC', 'IPSL-CM5A-MR', 
		'CMCC-CESM', 'IPSL-CM5B-LR', 'CNRM-CM5', 'MPI-ESM-LR', 
		'GFDL-ESM2G', 'MPI-ESM-MR', 'GFDL-ESM2M', 'MRI-ESM1', 
		'HadGEM2-CC', 'NorESM1-ME', 'HadGEM2-ES',]
			

TAKAHASHITypes 	= ['pCO2',]

GEOTRACESTypes 	= ['iron',]

BGCmodels 	= ['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10',]

Seasons		= ['JFM','AMJ','JAS','OND'] 

Hemispheres	=['NorthHemisphere','SouthHemisphere',]

months = [m for m in month_name if m]	# Because months starts at 1, and 0 is empty.
OceanMonth_names = [o+m for o in Ocean_names for m in months]
OceanSeason_names = [o+s for o in Ocean_names for s in Seasons]
HemispheresMonths = [h+m for h in Hemispheres for m in months] 	
SouthHemispheresMonths = [h+m for h in ['SouthHemisphere',] for m in months] 	
NorthHemispheresMonths = [h+m for h in ['NorthHemisphere',] for m in months] 	


#####
# Get Match Type: - _obsolete now
#	
def getmt_obsolete(loadYaml=False): 
	"""
	getmt(): Get Match Type. 
		returns a Nested Dictionary with the following structure options:
		
	Typical usage:
		mt[ModelName or Data source][variable name] = [list of variable names as they appears in the netcdf]
		
		where the mt dict always includes the coordinate names as their appear in the netcdf file too. ie:
			mt['GEOTRACES']['lat'] 		= 'Latitude'
	
		This assumes that all netcdf files from the same source kept the same coordinate names.
		Similarly, dimensions are treated differently than data.
		ie data are a list, dimensions are a string.
		
	Alternative usage:	
		It is possible to plot a combination of variables, such the sum of a series of values or the ratio of two.
		In some cases you may want to divide a value by 1000. 
		These cases work like this:
			mt[ModelName or Data source][variable name]['name'] = variable name 
			mt[ModelName or Data source][variable name]['convert'] = a function defined as defined here	
			mt[ModelName or Data source][variable name][vars] = [list of variables to combine in the convert function]
			mt[ModelName or Data source][variable name][new units] = 'mg m^-3' or whatever. (preferable from fancyUnits, below)
	"""
	assert False
	
	#####
	# Try to open a yaml file. 
	# By detault, this is not used.
	yamlFile = folder('yaml')+'matchMetadata.yaml'
	if exists(yamlFile) and loadYaml:
		print 'getmt:\tLoading mt file from ',yamlFile
		mt = YamlToDict(yamlFile,)
		return mt
		print 'getmt:\tCreating mt file:',yamlFile


	####
	# Some functions for maniulating data:
	def NoChange(nc,keys):	return nc.variables[keys[0]][:]
	def N2Biomass(nc,keys):	return nc.variables[keys[0]][:]* 79.573
	def mul1000(nc,keys):	return nc.variables[keys[0]][:]* 1000.
	def div1000(nc,keys):	return nc.variables[keys[0]][:]/ 1000.	
	def applymask(nc,keys):	return np.ma.masked_where(nc.variables[keys[1]][:]==0.,nc.variables[keys[0]][:])
	def sums(nc,keys):	
		a = nc.variables[keys[0]][:]
		for k in keys[1:]:a += nc.variables[k][:]
		return a 
	def oxconvert(nc,keys): return nc.variables[keys[0]][:] *44.661
	# 1 ml/l = 103/22.391 = 44.661 umol/l
	# http://ocean.ices.dk/Tools/UnitConversion.aspx
		
	tdicts = {	'ZeroToZero': {i  :i     for i in xrange(12)},		
			'OneToOne':   {i+1:i+1   for i in xrange(12)},
			'OneToZero':  {i+1:i     for i in xrange(12)},
			'ZeroToOne':  {i  :i+1   for i in xrange(12)},			
		}			
	#####
	# Models:
	mt = AutoVivification() # match type
	mt['ERSEM']['bac'] 			= ['B1c',]
	mt['ERSEM']['mesozoo'] 			= ['Z4c',]
	mt['ERSEM']['diatoms'] 			= ['P1c',]
	mt['ERSEM']['picophyto'] 		= ['P3c',]
	mt['ERSEM']['microzoo'] 		= ['Z5c',]			 	
	mt['ERSEM']['pCO2'] 			= ['pCO2w',]#'fAirSeaC',]
	mt['ERSEM']['chl'] 			= ['chl',]
    	mt['ERSEM']['phosphate']  		= ['N1p',]	
    	mt['ERSEM']['nitrate']  		= ['N3n',]	    	
    	mt['ERSEM']['silicate']  		= ['N5s',]	
    	mt['ERSEM']['oxygen']	  		= ['O2o',]	    	
	mt['ERSEM']['iron'] 			= ['N7f',]
	mt['ERSEM']['intpp'] 			= ['intpp',]	
	mt['ERSEM']['t']			= 'index_t'	
	mt['ERSEM']['z'] 			= 'deptht'
	mt['ERSEM']['lat'] 			= 'nav_lat'
	mt['ERSEM']['lon'] 			= 'nav_lon'
	mt['ERSEM']['cal'] 			= '365_day'
	
	    		    	
	mt['NEMO']['temperature'] 		= ['votemper',]	
 	mt['NEMO']['salinity'] 			= ['vosaline',]				
	mt['NEMO']['mld'] 			= ['somxl010',]
	mt['NEMO']['mld'] 			= ['somxl010',]	
	mt['NEMO']['mld_DT02'] 			= ['somxl010',]
	mt['NEMO']['mld_DR003'] 		= ['somxl010',]
	mt['NEMO']['mld_DReqDTm02'] 		= ['somxl010',]
	mt['NEMO']['U'] 			= ['vozocrtx',]
	mt['NEMO']['V'] 			= ['vomecrty',]				
	mt['NEMO']['t'] 			= 'index_t'	
	mt['NEMO']['z'] 			= 'deptht'
	mt['NEMO']['lat'] 			= 'nav_lat'
	mt['NEMO']['lon'] 			= 'nav_lon'
	mt['NEMO']['cal']			= '365_day'


	#mt['MEDUSA']['chl']['name'] 		=  'CHL'
	#mt['MEDUSA']['chl']['vars'] 		=  ['CHN','CHD']		
	#mt['MEDUSA']['chl']['convert'] 		=  sums
	#mt['MEDUSA']['chl']['units'] 		=  'mg C/m^3'
	#mt['MEDUSA']['chl'] 			=  ['CHL',]	
	
	mt['MEDUSA']['chl']['name'] 		=  'CHL'
	mt['MEDUSA']['chl']['vars'] 		=  ['CHL',]		
	mt['MEDUSA']['chl']['convert'] 		=  NoChange
	mt['MEDUSA']['chl']['units'] 		=  'mg C/m^3'
	#mt['MEDUSA']['chl'] 			=  ['CHL',]	
	
	mt['MEDUSA']['diatoms']['name'] 	=  'PHD'
	mt['MEDUSA']['diatoms']['vars'] 	=  ['PHD',]		
	mt['MEDUSA']['diatoms']['convert'] 	=  N2Biomass
	mt['MEDUSA']['diatoms']['units'] 	=  'mg C/m^3'
	
	mt['MEDUSA']['iron']['name']		=  'FER'
	mt['MEDUSA']['iron']['vars']		=  ['FER',]
	mt['MEDUSA']['iron']['convert']		=  mul1000
	mt['MEDUSA']['iron']['units']		=  'umol F/m^3'		

	mt['MEDUSA']['mesozoo']['name'] 	=  'ZME'	
	mt['MEDUSA']['mesozoo']['vars'] 	=  ['ZME',]				
	mt['MEDUSA']['mesozoo']['convert'] 	=  N2Biomass	
	mt['MEDUSA']['mesozoo']['units'] 	=  'mg C/m^3'		
	
	mt['MEDUSA']['microzoo']['name'] 	=  'ZMI'
	mt['MEDUSA']['microzoo']['vars'] 	=  ['ZMI',]					
	mt['MEDUSA']['microzoo']['convert'] 	=  N2Biomass			
	mt['MEDUSA']['microzoo']['units'] 	=  'mg C/m^3'


	dmsd= {'dms_and':'anderson','dms_ara':'aranamit','dms_hal':'halloran','dms_sim':'simodach'} 
	for dms in ['dms_and','dms_ara','dms_hal','dms_sim']:	
		mt['MEDUSA'][dms]['name'] 	=  dmsd[dms]
		mt['MEDUSA'][dms]['vars'] 	=  [dmsd[dms],]					
		mt['MEDUSA'][dms]['convert'] 	=  NoChange			
		mt['MEDUSA'][dms]['units'] 	=  'umol / m3'
	dmsd= {'dms_p_and':'anderson','dms_p_ara':'aranamit','dms_p_hal':'halloran','dms_p_sim':'simodach'} 
	for d in ['dms_p_sim','dms_p_and','dms_p_ara','dms_p_hal',]:	
	  for i in ['','1','2']:
	  	dms = d+i
		mt['MEDUSA'][dms]['name'] 	=  dmsd[d]+i
		mt['MEDUSA'][dms]['vars'] 	=  [dmsd[d]+i,]					
		mt['MEDUSA'][dms]['convert'] 	=  NoChange			
		mt['MEDUSA'][dms]['units'] 	=  'umol / m3'
		

#	mt['MEDUSA']['dms'][

    	mt['MEDUSA']['nitrate'] 	 	= ['DIN',]
    	mt['MEDUSA']['pCO2']	  		= ['OCN_PCO2',]		    	
    	mt['MEDUSA']['silicate']	  	= ['SIL',]			
    	mt['MEDUSA']['oxygen']	  		= ['OXY',]			    	
	mt['MEDUSA']['t'] 			= 'index_t'	
	mt['MEDUSA']['z'] 			= 'deptht'
	mt['MEDUSA']['lat'] 			= 'nav_lat'
	mt['MEDUSA']['lon'] 			= 'nav_lon'
	mt['MEDUSA']['cal'] 			= '365_day'
	#mt['Medusa']				= mt['MEDUSA']



    	mt['IMARNET']['chl']['vars']  		= ['chl',]
    	mt['IMARNET']['chl']['name']  		= 'chl'
	mt['IMARNET']['chl']['units'] 		= 'mg Chl/m3'
	mt['IMARNET']['chl']['convert'] 	=  NoChange		    	
    	mt['IMARNET']['nitrate']['vars']  	= ['no3',]
    	mt['IMARNET']['nitrate']['name']  	= 'Nitrate'
	mt['IMARNET']['nitrate']['units'] 	= 'mmol/m^3'
	mt['IMARNET']['nitrate']['convert'] 	=  NoChange				

    	mt['IMARNET']['phosphate']['vars']  	= ['po4',]
    	mt['IMARNET']['phosphate']['name']  	= 'Phosphate'
	mt['IMARNET']['phosphate']['units'] 	= 'mmol/m^3'
	mt['IMARNET']['phosphate']['convert'] 	=  NoChange	
		
    	mt['IMARNET']['pCO2']['vars']  		= ['spco2',]
    	mt['IMARNET']['pCO2']['name']  		= 'Partial pressure pCO2'
	mt['IMARNET']['pCO2']['units'] 		= 'uatm'
	mt['IMARNET']['pCO2']['convert'] 	=  NoChange	

    	mt['IMARNET']['intpp']['vars']  	= ['intpp',]
    	mt['IMARNET']['intpp']['name']  	= 'Integrated Primary Production'
	mt['IMARNET']['intpp']['units'] 	= 'gC/m2/d'
	mt['IMARNET']['intpp']['convert'] 	=  mul1000	

    	mt['IMARNET']['oxygen']['vars']  	= ['o2',]
    	mt['IMARNET']['oxygen']['name']  	= 'Dissolved Oxygen'
	mt['IMARNET']['oxygen']['units'] 	= 'mmol/m^3'
	mt['IMARNET']['oxygen']['convert'] 	=  NoChange	
		

    	mt['IMARNET']['silicate']['vars']  	= ['si',]
    	mt['IMARNET']['silicate']['name']  	= 'Silicate'
	mt['IMARNET']['silicate']['units'] 	= 'mmol/m^3'
	mt['IMARNET']['silicate']['convert'] 	=  NoChange
				
	mt['IMARNET']['t'] 			= 'index_t'	
	mt['IMARNET']['z'] 			= 'index_z' 
	mt['IMARNET']['lat'] 			= 'nav_lat'
	mt['IMARNET']['lon'] 			= 'nav_lon'
	mt['IMARNET']['cal'] 			= '365_day'
	for model in BGCmodels:
		mt['IMARNET_'+model] = mt['IMARNET']
		mt['iMarNet_'+model] = mt['IMARNET']

	mt['CMIP5']['t'] 			= 'time'	
	mt['CMIP5']['z'] 			= 'lev' 
	mt['CMIP5']['lat'] 			= 'lat'
	mt['CMIP5']['lon'] 			= 'lon'
	mt['CMIP5']['cal'] 			= '365_day'
    	mt['CMIP5']['oxygen']['vars']  	= ['o2',]
    	mt['CMIP5']['oxygen']['name']  	= 'Dissolved Oxygen'
	mt['CMIP5']['oxygen']['units'] 	= 'mmol/m^3'
	mt['CMIP5']['oxygen']['convert'] =  NoChange	

    	mt['CMIP5']['U']['vars']  	= ['uo',]
    	mt['CMIP5']['U']['name']  	= 'Zonal Velocity'
	mt['CMIP5']['U']['units'] 	= 'm/s'
	mt['CMIP5']['U']['convert'] 	=  NoChange
	
    	mt['CMIP5']['V']['vars']  	= ['vo',]
    	mt['CMIP5']['V']['name']  	= 'Meridional Velocity'
	mt['CMIP5']['V']['units'] 	= 'm/s'
	mt['CMIP5']['V']['convert'] 	=  NoChange
					
	CMIP5models = [ 'MEDUSA','ERSEM','BNU-ESM', 'IPSL-CM5A-LR', 'CESM1-BGC', 'IPSL-CM5A-MR', 
			'CMCC-CESM', 'IPSL-CM5B-LR', 'CNRM-CM5', 'MPI-ESM-LR', 
			'GFDL-ESM2G', 'MPI-ESM-MR', 'GFDL-ESM2M', 'MRI-ESM1', 
			'HadGEM2-CC', 'NorESM1-ME', 'HadGEM2-ES',]
	for model in CMIP5models:
		mt['CMIP5_'+model] = mt['CMIP5']
		mt['cmip5_'+model] = mt['CMIP5']		
	#models= ['MEDUSA','NEMO']
	
	
			
	#####
	# Data:
	mt['MAREDAT']['bac'] 			= ['BIOMASS',]
	mt['MAREDAT']['mesozoo'] 		= ['BIOMASS',]
	mt['MAREDAT']['diatoms'] 		= ['BIOMASS',]
	mt['MAREDAT']['picophyto'] 		= ['BIOMASS',]
	mt['MAREDAT']['microzoo'] 		= ['BIOMASS',]
	mt['MAREDAT']['PP'] 			= ['PP',]
	mt['MAREDAT']['intpp'] 			= ['PPint',]	
	mt['MAREDAT']['chl']['name']		= 'Chlorophylla'
	mt['MAREDAT']['chl']['vars']		= ['Chlorophylla',]
	mt['MAREDAT']['chl']['convert']		= div1000	
	mt['MAREDAT']['chl']['units']		= ['ug/L',]
	mt['MAREDAT']['t'] 			= 'index_t'
	mt['MAREDAT']['z'] 			= 'DEPTH'
	mt['MAREDAT']['lat'] 			= 'LATITUDE'
	mt['MAREDAT']['lon'] 			= 'LONGITUDE'
	mt['MAREDAT']['cal'] 			= 'standard'
	mt['MAREDAT']['tdict']			= tdicts['ZeroToZero']
	
	#mt['Maredat'] 				= mt['MAREDAT']
			
	mt['WOA']['temperature'] 		= ['t_an',]#'t_mn',
  	mt['WOA']['salinity'] 			= ['s_an',]#'s_mn',
  	mt['WOA']['nitrate'] 			= ['n_an',]#'s_mn',  	
	mt['WOA']['silicate'] 			= ['i_an',]#'i_mn',
	mt['WOA']['phosphate'] 			= ['p_an',]#'p_mn',	    	  		
	#mt['WOA']['oxygen'] 			= ['o_an',]#'p_mn',	
    	mt['WOA']['oxygen']['vars']  		= ['o_an',]
    	mt['WOA']['oxygen']['name']  		= 'o_an'
	mt['WOA']['oxygen']['units'] 		= 'mmol/m^3'
	mt['WOA']['oxygen']['convert'] 		=  oxconvert
		    	  			
	mt['WOA']['t'] 				= 'index_t'
	mt['WOA']['z'] 				= 'depth'
	mt['WOA']['lat'] 			= 'lat'
	mt['WOA']['lon'] 			= 'lon'
	mt['WOA']['cal'] 			= 'standard'    
	mt['WOA']['tdict']			= tdicts['ZeroToZero']
	

	mt['TAKAHASHI']['pCO2'] 		= ['PCO2_SW',]#'DELTA_PCO2',]	'TFLUXSW06',
	mt['TAKAHASHI']['t'] 			= 'TIME'
	mt['TAKAHASHI']['z'] 			= 'index_z'
	mt['TAKAHASHI']['lat'] 			= 'LAT'
	mt['TAKAHASHI']['lon'] 			= 'LON'
	mt['TAKAHASHI']['cal'] 			= 'standard'
	mt['TAKAHASHI']['tdict']		= tdicts['ZeroToZero']
	#mt['intPP']['intPP']			= ['PPint',]	
								
	mt['GEOTRACES']['iron']			= ['Fe_D_CONC_BOTTLE',]#'Fe_D_CONC_BOTTLE_FIA','Fe_S_CONC_BOTTLE',]
	mt['GEOTRACES']['t']			= 'MONTH'
	mt['GEOTRACES']['z'] 			= 'DEPTH'
	mt['GEOTRACES']['lat'] 			= 'Latitude'
	mt['GEOTRACES']['lon'] 			= 'Longitude'
	mt['GEOTRACES']['cal'] 			= 'standard'
	mt['GEOTRACES']['tdict']		= tdicts['OneToZero']
		
	mt['IFREMER']['mld']['name']		= 'mld'
	mt['IFREMER']['mld']['vars']		= ['mld','mask']
	mt['IFREMER']['mld']['convert']		= applymask	
	mt['IFREMER']['mld']['units']		= ['m',]
	#mt['IFREMER']['mld_DT02']		= ['mld','mask']
	#mt['IFREMER']['mld_DR003']		= ['mld','mask']
	#mt['IFREMER']['mld_DReqDTm02']		= ['mld','mask']
	mt['IFREMER']['t'] 			= 'index_t'
	mt['IFREMER']['z'] 			= 'index_z'
	mt['IFREMER']['lat'] 			= 'lat'
	mt['IFREMER']['lon'] 			= 'lon'
	mt['IFREMER']['cal'] 			= 'standard'	
	mt['IFREMER']['tdict']			= tdicts['ZeroToZero']
	
	
#	mt['DRIFTERS']['U']			= ['U',]
#	mt['DRIFTERS']['V']			= ['V',]
#	mt['DRIFTERS']['t']			= 'Time'
#	mt['DRIFTERS']['z'] 			= 'Depth'
#	mt['DRIFTERS']['lat'] 			= 'Lat'
#	mt['DRIFTERS']['lon'] 			= 'Lon'
#	mt['DRIFTERS']['cal'] 			= 'standard'
#	mt['DRIFTERS']['tdict']			= tdicts['OneToZero']	

	mt['GODAS']['U']			= ['uo',]
	mt['GODAS']['V']			= ['vo',]
	mt['GODAS']['t']			= 'time'
	mt['GODAS']['z'] 			= 'lev'
	mt['GODAS']['lat'] 			= 'lat'
	mt['GODAS']['lon'] 			= 'lon'
	mt['GODAS']['cal'] 			= 'standard'
	mt['GODAS']['tdict']			= tdicts['OneToZero']
	
		
	for dms in ['dms_and','dms_ara','dms_hal','dms_sim']:	
		mt['LANA'][dms]['name'] 	=  'lanaetal'
		mt['LANA'][dms]['vars'] 	=  ['lanaetal',]					
		mt['LANA'][dms]['convert'] 	=  NoChange			
		mt['LANA'][dms]['units'] 	=  'umol / m3'
	mt['LANA']['t'] 		= 'index_t'
	mt['LANA']['z'] 		= 'deptht'
	mt['LANA']['lat'] 		= 'nav_lat'
	mt['LANA']['lon'] 		= 'nav_lon'
	mt['LANA']['tdict']		= tdicts['ZeroToZero']	

	for d in ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']:
	  for i in ['','1','2']:
	  	dms = d+i
		mt['LANA_p'][dms]['name'] 	=  'DMS'
		mt['LANA_p'][dms]['vars'] 	=  ['DMS',]					
		mt['LANA_p'][dms]['convert'] 	=  NoChange			
		mt['LANA_p'][dms]['units'] 	=  'nM'
		
	mt['LANA_p']['t'] 		= 'index_t'
	mt['LANA_p']['z'] 		= 'deptht'
	mt['LANA_p']['lat'] 		= 'nav_lat'
	mt['LANA_p']['lon'] 		= 'nav_lon'
	mt['LANA_p']['tdict']		= tdicts['ZeroToZero']	
	#mt['PP']['PP'] 		= ['PP',]
	#'AutoVivToYaml(mt,yamlFile)
	return mt	
	
#def getShelveName(workingDir,model,jobID,year,name,newslice,xkey='*',ykey='*'):
	#workingDir  = "/data/euryale7/scratch/ledm/ukesm_postProcessed/"
	#Model = "ERSEM"
	#jobID = "xhonp"
	#year="1998"
	#name="pCO2"
#	shelvename = workingDir + model+'-'+jobID+'-'+year+'/'+name+'_'+newslice+'_'xkey+'vs'+ykey+'.shelve'
	#return shelve

def getLongName(text,debug=False):
	if debug: print "Getting long name:",text
	if type(text) in [type(['a','b',]),type(('a','b',))]:
		return ' '.join([getLongName(t) for t in text])
		#out = ''
		#for t in text:out+=getLongName(t)+' '
		#return out
	
	#mt = getmt()
	
	noChange = ['Surface',]
	#noChange.extend(mt.keys())
	
	if text in noChange:return text
	
	firstLetterCaps = ['temperature', "salinity", "nitrate", "phosphate" ,'silicate','iron',]
	
	if text.lower() in firstLetterCaps:	return text.title()


  	if text == 'Transect':		return "Transect"	
  	if text == 'PTransect':		return "Pacific Transect"
#  	if text == 'Surface':		return "Surface"
  	if text == '100m':		return "100m deep"  	
  	if text == '200m':		return "200m deep"  	  	  	
  	if text == '500m':		return "500m deep"  	
  	if text == '1000m':		return "1000m deep"  	

#  	if text == 'temperature':	return "Temperature" 
# 	if text == 'salinity':		return "salinity" 
#  	if text == 'nitrate':		return "nitrate" 
#  	if text == 'phosphate':		return "phosphate" 
#  	if text == 'silicate':		return "silicate"   	  	  	
	  		  	  	
  	if text == 'temperaturePTransect':	return "Pacific Transect Temperature"
  	if text == 'temperatureSOTransect':	return "Southern Ocean Transect Temperature"  	
  	if text == 'temperatureSurface':	return "Surface Temperature"
  	if text == 'temperatureAll':		return "Temperature"  	
  	if text == 'temperature100m':		return "Temperature (100m deep)"  	
  	if text == 'temperature200m':		return "Temperature (200m deep)"  
  	if text == 'temperature500m':		return "Temperature (500m deep)"  	
  	if text == 'temperature1000m':		return "Temperature (1000m deep)"  	  		  	  	
  	if text == 'TemperaturePTransect':	return "Pacific Transect Temperature"
  	if text == 'TemperatureSOTransect':	return "Southern Ocean Transect Temperature"  	
  	if text == 'TemperatureSurface':	return "Surface Temperature"
  	if text == 'TemperatureAll':		return "Temperature"  	
  	if text == 'Temperature100m':		return "Temperature (100m deep)"  	
  	if text == 'Temperature200m':		return "Temperature (200m deep)"  
  	if text == 'Temperature500m':		return "Temperature (500m deep)"  	
  	if text == 'Temperature1000m':		return "Temperature (1000m deep)"  
  	  	
  	if text == 'salinityPTransect':		return "Pacific Transect Salinity"
  	if text == 'salinitySOTransect':		return "Southern Ocean Transect Salinity"  	
  	if text == 'salinitySurface':		return "Surface Salinity"
  	if text == 'salinityAll':		return "Salinity"  	
  	if text == 'salinity100m':		return "Salinity (100m deep)"  	
  	if text == 'salinity200m':		return "Salinity (200m deep)"  	  	  	
  	if text == 'salinity500m':		return "Salinity (500m deep)"  	
  	if text == 'salinity1000m':		return "Salinity (1000m deep)"  	 	  	
        if text == 'SalinityTransect':          return "Atlantic Transect Salinity"
  	if text == 'SalinityPTransect':		return "Pacific Transect Salinity"
  	if text == 'SalinitySOTransect':	return "Southern Ocean Transect Salinity"  	
  	if text == 'SalinitySurface':		return "Surface Salinity"
  	if text == 'SalinityAll':		return "Salinity"  	
  	if text == 'Salinity100m':		return "Salinity (100m deep)"  	
  	if text == 'Salinity200m':		return "Salinity (200m deep)"  	  	  	
  	if text == 'Salinity500m':		return "Salinity (500m deep)"  	
  	if text == 'Salinity1000m':		return "Salinity (1000m deep)"  	
  	  	
  	if text == 'nitratePTransect':	return "Pacific Transect Nitrate (WOA14)"
  	if text == 'nitrateSOTransect':	return "Southern Ocean Transect Nitrate (WOA14)"  	
  	if text == 'nitrateSurface':	return "Surface Nitrate (WOA14)"
  	if text == 'nitrateAll':	return "Nitrate (WOA14)"  	
  	if text == 'nitrate100m':	return "Nitrate (100m deep)"  	  	
  	if text == 'nitrate200m':	return "Nitrate (200m deep)"  
  	if text == 'nitrate500m':	return "Nitrate (500m deep)"  	  	
  	if text == 'NitrateTransect':	return "Atlantic Transect Nitrate (WOA14)"  	
  	if text == 'NitratePTransect':	return "Pacific Transect Nitrate (WOA14)"  	
  	if text == 'NitrateSOTransect':	return "Southern Ocean Transect Nitrate (WOA14)"
  	if text == 'NitrateSurface':	return "Surface Nitrate (WOA14)"
  	if text == 'NitrateAll':	return "Nitrate (WOA14)"  	
  	if text == 'Nitrate100m':	return "Nitrate (100m deep)"  	  	
  	if text == 'Nitrate200m':	return "Nitrate (200m deep)"  
  	if text == 'Nitrate500m':	return "Nitrate (500m deep)"  	
  	 
  	  		  	  	
  	if text == 'phosphatePTransect':return "Pacific Transect Phosphate"
  	if text == 'phosphateSOTransect':return "Southern Ocean Transect Phosphate"  	  	
  	if text == 'phosphateSurface':	return "Surface Phosphate"
  	if text == 'phosphateAll':	return "Phosphate"  	
  	if text == 'phosphate100m':	return "Phosphate (100m deep)"  	  	
  	if text == 'phosphate200m':	return "Phosphate (200m deep)"    	
  	if text == 'phosphate500m':	return "Phosphate (500m deep)"    	  	
  	if text == 'PhosphatePTransect':return "Pacific Transect Phosphate"
  	if text == 'PhosphateSOTransect':return "Southern Ocean Transect Phosphate"  	
  	if text == 'PhosphateSurface':	return "Surface Phosphate"
  	if text == 'PhosphateAll':	return "Phosphate"  	
  	if text == 'Phosphate100m':	return "Phosphate (100m deep)"  	  	
  	if text == 'Phosphate200m':	return "Phosphate (200m deep)"    	
  	if text == 'Phosphate500m':	return "Phosphate (500m deep)"   
  	 	  	
  	if text == 'silicatePTransect':	return "Pacific Transect Silicate"
  	if text == 'silicateSurface':	return "Surface Silicate"
  	if text == 'silicateAll':	return "Silicate"  	
  	if text == 'silicate100m':	return "Silicate (100m deep)"  	  	
  	if text == 'silicate200m':	return "Silicate (200m deep)"    	
  	if text == 'silicate500m':	return "Silicate (500m deep)"   
  	if text == 'SilicateTransect':	return "Atlantic Transect Silicate"  	
  	if text == 'SilicatePTransect':	return "Pacific Transect Silicate"
  	if text == 'SilicateSurface':	return "Surface Silicate"
  	if text == 'SilicateAll':	return "Silicate"  	
  	if text == 'Silicate100m':	return "Silicate (100m deep)"  	  	
  	if text == 'Silicate200m':	return "Silicate (200m deep)"    	
  	if text == 'Silicate500m':	return "Silicate (500m deep)"    	

	if text == 'OxygenSurface':	return "Surface Oxygen"  
	if text == 'OxygenTransect':	return "Atlantic Transect Oxygen"  
	if text == 'OxygenPTransect':	return "Pacific Transect Oxygen"  
	if text == 'OxygenSOTransect':	return "Southern Ocean Transect Oxygen"  	
	if text == 'Oxygen100m':		return "Oxygen (100m deep) "  
	if text == 'Oxygen200m':		return "Oxygen (200m deep)"
	if text == 'Oxygen500m':		return "Oxygen (500m deep)"  
	if text == 'Oxygen1000m':		return "Oxygen (1000m deep)"  

	if text == 'Alkalinity':	return "Alkalinity"  
  	if text == 'alkalinityPTransect':return "Pacific Transect Alkalinity"
  	if text == 'alkalinitySurface':	return "Surface Alkalinity"
  	if text == 'alkalinityAll':	return "Alkalinity"  	
  	if text == 'alkalinity100m':	return "Alkalinity (100m deep)"  	  	
  	if text == 'alkalinity200m':	return "Alkalinity (200m deep)"    	
  	if text == 'alkalinity500m':	return "Alkalinity (500m deep)"   
  	if text == 'AlkalinityTransect':return "Atlantic Transect Alkalinity"
  	if text == 'AlkalinityPTransect':return "Pacific Transect Alkalinity"
  	if text == 'AlkalinitySOTransect':return "South Ocean Transect Alkalinity"  	
  	if text == 'AlkalinitySurface':	return "Surface Alkalinity"
  	if text == 'AlkalinityAll':	return "Alkalinity"  	
  	if text == 'Alkalinity100m':	return "Alkalinity (100m deep)"  	  	
  	if text == 'Alkalinity200m':	return "Alkalinity (200m deep)"    	
  	if text == 'Alkalinity500m':	return "Alkalinity (500m deep)"    	

  	if text == 'DICTransect':	return "Atlantic Transect DIC"
  	if text == 'DICSOTransect':	return "Southern Ocean Transect DIC"
  	if text == 'DICPTransect':	return "Pacific Transect DIC"
  	if text == 'DICSurface':	return "Surface DIC"
  	if text == 'DICAll':	return "DIC"  	
  	if text == 'DIC100m':	return "DIC (100m deep)"  	  	
  	if text == 'DIC200m':	return "DIC (200m deep)"    	
  	if text == 'DIC500m':	return "DIC (500m deep)"   
  	
  	if text == 'AMOC_26N':	return "AMOC 26N"     	
  	if text == 'AMOC_32S':	return "AMOC 23S"     	  	
  		
	if text == 'percentiles':	return "Time series"  
	if text == 'Chlorophyll':	return "Chlorophyll"  
	if text == 'Chlorophyll_cci':	return "Chlorophyll (CCI)"  
	if text == 'Chlorophyll_pig':	return "Chlorophyll (Pigments)"  		
	if text == 'cci':		return "CCI"  
	if text == 'pig':		return "pigments"  
	if text == 'exportRatio':	return "Export Ratio"  
	if text == 'ExportRatio':	return "Export Ratio"  	
	if text == 'LocalExportRatio':	return "Export Ratio"  		
	
	if text == 'IntegratedPrimaryProduction':	return "Integrated Primary Production"  
	if text == 'IntegratedPrimaryProduction_OSU':	return "Integrated Primary Production (OSU)"  
	if text == 'IntegratedPrimaryProduction':	return "Integrated Primary Production"
	if text == 'TotalIntegratedPrimaryProduction':	return "Total Integrated Primary Production"	
	  		
	if text == 'OSU':		return "(OSU)"  
	if text == '1x1':		return "(iMarNet)"  
	if text == 'Oxygen':		return "Oxygen"  
	if text == 'AirSeaFluxCO2':	return "Air Sea CO2 Flux"  
	if text == 'TotalAirSeaFluxCO2':return "Total Air Sea CO2 Flux"  
	if text == 'DIC':		return "DIC"  
	if text == 'DICSurface':	return "Surface DIC"  
	if text == 'DICTransect':	return "Atlantic Transect DIC"  
	if text == 'DICPTransect':	return "Pacific Transect DIC"  
	if text == 'DICSOTransect':	return "Southern Ocean Transect DIC"  	
	if text == 'DIC100m':		return "DIC (100m deep) "  
	if text == 'DIC200m':		return "DIC (200m deep)"
	if text == 'DIC500m':		return "DIC (500m deep)"  
	if text == 'DIC1000m':		return "DIC (1000m deep)"  	
  	if text == 'OMZ':	 	return 'Oxygen minimum zone' 
  	if text == 'TotalOMZVolume':	return 'Total Oxygen minimum zone volume (<20 mmol O2/m^3)'   	
  	if text == 'TotalOMZVolume50':	return 'Total Oxygen minimum zone volume (<50 mmol O2/m^3)'   	  	
  	if text == 'OMZThickness':	return 'Oxygen minimum zone thickness (<20 mmol O2/m^3)'  
  	
  	if text == 'DrakePassageTransport':	return 'Drake Passage Transport'    	
  	if text == 'DPT':			return 'Drake Passage Transport'    	  		
	
	if text == 'hist':		return "Histogram"  
	if text == 'scatter':		return "scatter diagram"  	
	if text == 'robinquad':		return "Maps"  
	if text == 'robinquad-cartopy':	return "Interpolated Map"
	if text == 'hov':		return "Hovmoeller"  			
	if text == '10-90pc':		return ""  				
	if text == 'Target':		return "Target Diagram"  	
	if text == 'Taylor':		return "Taylor Diagram"  	
	if text == 'SummaryTargets':	return "Summary Diagrams"  		
	if text == 'RobustTarget':	return "Robust Statisitcs Target Diagram" 			


  	if text == 'AtlanticTransect':	return "Atlantic Transect"    
  	if text == 'PacificTransect':	return "Pacific Transect"    	
  	if text == 'SouthernTransect':	return "Southern Transect"    	
  	if text == 'SOTransect':	return "Southern Transect"    
  	if text == 'Transect':		return "Atlantic Transect"    
  	if text == 'PTransect':		return "Pacific Transect"      	  		  	
  	
  	if text == '10N':		return "10 degree North Transect"    	
  	if text == '10S':		return "10 degree South Transect"    	  	  	  	  	
  	  	
  	#if text == 'nitTransect': return "Pacific Transect Nitrate (WOA40)"
  	#if text == 'nitSurface': return "Surface Nitrate (WOA40)"
  	
    	if text == 'GEOTRACES':		return "GEOTRACES"  		
  	if text == 'iron':		return "Iron"  	
  	if text == 'Fe_D_CONC_BOTTLE':	return "Iron (Dissolved)"  		

  	if text ==  'WOA': 	return 'WOA'  	  	
  	if text ==  't_mn': 	return 'Mean Temperature'  	
  	if text ==  't_an': 	return 'Temperature'  
  	if text ==  's_mn': 	return 'Mean Salinity'  	
  	if text ==  's_an': 	return 'Salinity'  
  	if text ==  'n_mn': 	return 'Mean Nitrate'  	
  	if text ==  'n_an': 	return 'Nitrate'  
  	if text ==  'p_mn': 	return 'Mean Phosphate'  	
  	if text ==  'p_an': 	return 'Phosphate'  
  	if text ==  'p_an': 	return 'Oxygen'    	
  	if text ==  'i_mn': 	return 'Mean Silicate'  	
  	if text ==  'i_an': 	return 'Silicate'  
  	  	
  	if text == 'NorthernTotalIceExtent': 	return 'Northern Hemisphere Ice Extent'
  	if text == 'SouthernTotalIceExtent': 	return 'Southern Hemisphere Ice Extent'
  	if text == 'TotalIceExtent': 		return 'Total Ice Extent'  	  	

  	if text == 'NorthernTotalIceArea': 	return 'Northern Hemisphere Ice Area'
  	if text == 'SouthernTotalIceArea': 	return 'Southern Hemisphere Ice Area'
  	if text == 'TotalIceArea': 		return 'Total Ice Area'  
  	  		
	if text == 'Seasons':   return 'Seasons'
	if text in ['JFM','AMJ','JAS','OND',]:return text   	
  	
     	if text == 'IFREMER':	return "IFREMER"
  	if text == 'mld':		return 'Mixed Layer Depth'     	
  	if text == 'MLD':		return 'Mixed Layer Depth'     	  	
  	if text == 'mld_DT02':		return 'MLD: Fixed Threshold Temperature '
  	if text == 'mld_DR003':		return 'MLD: Fixed Threshold Density'
  	if text == 'mld_DReqDTm02':	return 'MLD: Variable Threshold Density'  	  	

 	if text == 'anderson':		return 'Anderson et al.'
 	if text == 'dms_and':		return 'Anderson et al.'
 	if text == 'dms_andSurface':	return 'Anderson et al.'
 	if text == 'dms_p_and':		return 'Anderson et al.'
 	if text == 'dms_p_andSurface':	return 'Anderson et al.'
 	if text == 'dms_p_and1':	return 'DMS (Anderson - all CHL)'
 	if text == 'dms_p_and1Surface':	return 'DMS (Anderson - all CHL)'
 	if text == 'dms_p_and2':	return 'DMS (Anderson - CHN only)'
 	if text == 'dms_p_and2Surface':	return 'DMS (Anderson - CHN only)' 	

 	if text == 'aranamit':		return 'Aranami et al.'
 	if text == 'dms_ara':		return 'Aranami et al.'
 	if text == 'dms_araSurface':	return 'Aranami et al.'
 	if text == 'dms_p_ara':		return 'Aranami et al.'
 	if text == 'dms_p_araSurface':	return 'Aranami et al.'
 	if text == 'dms_p_ara1':	return 'DMS (Aranamit - all CHL)'
 	if text == 'dms_p_ara1Surface':	return 'DMS (Aranamit - all CHL)'
 	if text == 'dms_p_ara2':	return 'DMS (Aranamit - CHN only)'
 	if text == 'dms_p_ara2Surface':	return 'DMS (Aranamit - CHN only)' 	
 	 	
 	if text == 'halloran':		return 'Halloran et al.' 	
 	if text == 'dms_hal':		return 'Halloran et al.' 	
 	if text == 'dms_halSurface':	return 'Halloran et al.' 	
 	if text == 'dms_p_hal':		return 'Halloran et al.' 	
 	if text == 'dms_p_halSurface':	return 'Halloran et al.' 	
 	if text == 'dms_p_hal1':	return 'DMS (Halloran - all CHL)'
 	if text == 'dms_p_hal1Surface':	return 'DMS (Halloran - all CHL)'
 	if text == 'dms_p_hal2':	return 'DMS (Halloran - CHN only)'
 	if text == 'dms_p_hal2Surface':	return 'DMS (Halloran - CHN only)'
 	 	 	
 	if text == 'simodach':		return 'Simo & Dach'
 	if text == 'dms_sim':		return 'Simo & Dach'
 	if text == 'dms_simSurface':	return 'Simo & Dach'
 	if text == 'dms_p_sim':		return 'Simo & Dach'
 	if text == 'dms_p_simSurface':	return 'Simo & Dach'
 	if text == 'dms_p_sim1':	return 'DMS (Simodach - all CHL)'
 	if text == 'dms_p_sim1Surface':	return 'DMS (Simodach - all CHL)'
 	if text == 'dms_p_sim2':	return 'DMS (Simodach - CHN only)'
 	if text == 'dms_p_sim2Surface':	return 'DMS (Simodach - CHN only)'

 	if text in ['LANA',]:		return 'Lana et al. (extrapolated)'
 	if text == 'LANA_p':		return 'Lana et al. (pixels)'
 	 	 	 	 	
 	if text == 'lanaetal':		return 'DMS extrapolated (Lana et al. 2011)'
 	if text == 'DMS':		return 'DMS pixels (Lana et al. 2011)'
	if text == 'DMS_p':		return 'DMS (pixels)'
	if text == 'DMS_e':		return 'DMS (extrapolated)'	
 
  	if text == 'All':	return 'Global'
  	if text == 'Best':	return 'Best'  	
  	if text == 'Standard':	return 'Standard'  	  	
  	if text =='SalArtifact':return 'Salinity Artifact (<15psu)'
  	if text =='NitArtifact':return 'Nitrogen Artifact'  	
  	if text =='OffAxis':	return 'Off Axis'  	
  	if text =='Depth':	return 'Depth >200m'  	
  	if text =='Shallow':	return 'Depth <200m'
  	
  	if text == 'BIOMASS':	return 'Biomass'
  	
  	if text =='Depth_0-10m':	return 'Depth <10m'
  	if text =='Depth_10-20m':	return '10m <= Depth < 20m'
  	if text =='Depth_20-50m':	return '20m <= Depth < 50m'
  	if text =='Depth_50-100m':	return '50m <= Depth < 100m'
  	if text =='Depth_100-500m':	return '100m <= Depth < 500m'
  	if text =='Depth_500m':		return 'Depth > 500m'  	  	
  	if text =='Depth_1000m':	return 'Depth > 1000m'  	  	
  	
  	if text =='Depth_0-50m':	return 'Depth <50m'
  	if text =='Depth_50-100m':	return '50m <= Depth < 100m'
  	if text =='Depth_100-200m':	return '100m <= Depth < 200m'
  	if text =='Depth_200-500m':	return '200m <= Depth < 500m'
  	if text =='Depth_500-1000m':	return '500m <= Depth < 1000m'
  	if text =='Depth_1000-2000m':	return '1000m <= Depth < 2000m'  	
  	if text =='Depth_2000m':	return 'Depth > 2000m'
  	
	
	  	
   	if text =='maskBelowBathy':	return 'Masked Below Bathymetery' 
   	if text =='OnShelf':		return 'On Shelf'
   	if text =='OffShelf':		return 'Off Shelf'    	   	  	
  	  	
   	if text =='nonZero':	return 'Non zero' 
   	if text =='aboveZero':	return ''#> zero'    	
   	if text =='1-99pc':	return '1-99 percentiles' 
   	if text =='5-95pc':	return '5-95 percentiles' 
   	if text =='0-99pc':	return 'up to 99th percentile' 
   
   	if text =='TypicalIron':	return 'Iron < 4 umol m^-3' 
   	   	
   	if text =='0-1pc':	return '0-1 percentiles' 
   	if text =='1-5pc':	return '1-5 percentiles' 
   	if text =='5-25pc':	return '5-25 percentiles' 
   	if text =='25-40pc':	return '25-40 percentiles' 
   	if text =='40-60pc':	return '40-60 percentiles' 
   	if text =='60-75pc':	return '60-75 percentiles' 
   	if text =='75-95pc':	return '75-95 percentiles' 
   	if text =='95-99pc':	return '95-99 percentiles' 
   	if text =='99-100pc':	return '99-100 percentiles' 
   	if text =='metricless':	return ''    	

	if text =='Overestimate_2sig': return "Overestimate 2 sigma"
	if text =='Overestimate_3sig': return "Overestimate 3 sigma"
	if text =='Underestimate_2sig': return "Overestimate 2 sigma"
	if text =='Underestimate_3sig': return "Overestimate 3 sigma"			
	   	
  	if text in ['Tropics','Temperate','Arctic','Surface','Depth', 'SalArtifact','OffAxis','pCO2',]:return text
	if text in ['Overestimate','Underestimate','Matched', 'Arctic','Tropics','Temperate', 'Surface',]:return text
	
  	if text == 'SurfaceNoArtics':	return ""


  	if text == 'Global':			return "Global"	
  	if text == 'Equator10':			return "Equator (+/-10)"	
  	if text == 'Remainder':			return "Oligotrophic Gyres"
  	if text == 'ArcticOcean':		return "Arctic Ocean"	
  	if text == 'NorthernSubpolarAtlantic':	return "Northern Subpolar Atlantic"	
  	if text == 'NorthernSubpolarPacific':	return "Northern Subpolar Pacific"	
  	if text == 'ignoreInlandSeas':		return "No Inland Seas"	
  	if text == 'SouthernOcean':		return "Southern Ocean"	  	  	  	  	  	  		
  	if text == 'regionless':		return ""		
  	if text == 'layerless':			return ""			
	
  	if text == 'NorthTemperate':	return "North Temperate"
  	if text == 'SouthTemperate':	return "South Temperate"  	
  	if text == 'NorthTropics':	return "North Tropics"
  	if text == 'SouthTropics':	return "South Tropics"  	
  	if text == 'NorthArctic':	return "North Arctic"
  	if text == 'Antarctic':		return "Antarctic"  	
  	if text == 'Equatorial':	return "Equatorial"  	
  	if text == 'AMT':		return "AMT"  	  	
  	if text == 'AMTTop40m':		return "AMT (Top 40m)"  	  	
   	if text == 'AMTTop200m':	return "AMT (Top 200m)"  	  	 	
   	
  	if text == 'BlackSea':		return "Black Sea"  
  	if text == 'RedSea':		return "Red Sea"  
  	if text == 'BalticSea':		return "Baltic Sea"  
   	if text == 'PersianGulf':	return "Persian Gulf"   	
  	  	  	 		
  	if text == 'ignoreBlackSea':	return "No Black Sea"  	  	
  	if text == 'ignoreRedSea':	return "No Red Sea"  	  	
  	if text == 'ignoreBalticSea':	return "No Baltic Sea"  
   	if text == 'ignorePersianGulf':	return "No Persian Gulf"   	  		  	
  	if text == 'ignoreInlandSeas':	return "No Inland Seas"
  	if text == 'ignoreMediteranean':return "No Mediteranean"  	

  	if text in ['Oceans',]:		return 'Oceans'	 	  	  	  	
  	if text in ['OceansMonths',]:		return 'Oceans Months'	 	  	  	  	  	
  	if text == 'ArcticOcean':	return "Arctic Ocean"	  		  	
  	if text == 'AntarcticOcean':	return "Antarctic Ocean" 		  	
  	if text == 'NorthAtlanticOcean':return "North Atlantic Ocean"
  	if text == 'SouthAtlanticOcean':return "South Atlantic Ocean"	  		  	
  	if text == 'NorthPacificOcean':	return "North Pacific Ocean"		  	
  	if text == 'SouthPacificOcean':	return "South Pacific Ocean"	  		  	
  	if text == 'EquatorialAtlanticOcean':	return "Equatorial Atlantic Ocean"	  		  	  	
  	if text == 'EquatorialPacificOcean':	return "Equatorial Pacific Ocean"	  		  	  	  	  	
  	if text == 'IndianOcean':	return "Indian Ocean"
  	if text == 'ignoreExtraArtics':	return "No Arctic Oceans (50 degrees)"  	
  	if text == 'ignoreMoreArtics':	return "No Arctic Oceans (60 degrees)"
  	if text == 'ignoreMidArtics':	return "No Arctic Oceans (65 degrees)"  
  	if text == 'ignoreArtics':	return "No Arctic Oceans (70 degrees)" 
  	if text == 'NorthHemisphere':	return "North Hemisphere"
  	if text == 'SouthHemisphere':	return "South Hemisphere"  
  	
  	if text == 'Top40m':	return "Top 40m"
  	if text == 'Top200m':	return "Top 200m"
  	if text == 'Top40mNoArtics':	return "Top 40m (No Arctics)"
  	if text == 'Top200mNoArtics':	return "Top 200m  (No Arctics)"

  	if text == 'Transect':	return "Pacifc Transect"  	
  	if text == 'AtlanticTransect':	return "Atlantic Transect"  
  		  	
  	if text == 'NoShelf':	return "No Shelf"  
  	if text == 'NoShelfTop40':	return "No Shelf (Top 40m)"    	
  	if text == 'NoShelfSurface':	return "No Shelf (Surface)"    	  	

  	  		  	
  	if text == 'picophyto':	return 'Picophytoplankton'
  	if text == 'microzoo':	return 'Microzooplankton'
  	if text == 'mesozoo':	return 'Mesozooplankton'
  	if text == 'diatoms':	return 'Diatoms'
  	if text ==  'bac': 	return 'Bacteria'
  	if text in  ['chl','Chlorophylla',]: 
  		return 'Chlorophyll'  	
  	if text in  ['chlSurface','ChlorophyllaSurface',]:  return 'Surface Chlorophyll'  	
  	  		
	if text in ['PCO2_SW', 'pCO2']:	return 'pCO2'

  	if text ==  'NEMO': 	return 'NEMO'  	  	
  	if text ==  'Nemo': 	return 'NEMO'
  	if text ==  'CICE': 	return 'CICE'  	  	  	
  	
  	if text ==  'ERSEM': 	return 'ERSEM'  	  	
  	if text ==  'ERSEM': 	return 'ERSEM'  	  	
  	if text ==  'ERSEM-1891': return 'ERSEM (1891)'    	  	
  	if text ==  'ERSEM-1893': return 'ERSEM (1893)'  	
  	if text ==  'ERSEM-1894': return 'ERSEM (1894)'
  	if text ==  'ERSEM-1895': return 'ERSEM (1895)'  	
  	if text ==  'ERSEM-1899': return 'ERSEM (1899)'  	  	
  	if text ==  'ERSEM-1909': return 'ERSEM (1909)'
  	if text ==  'ERSEM-1948': return 'ERSEM (1948)'  
  	if text ==  'ERSEM-1982': return 'ERSEM (1982)'    	
  	if text ==  'ERSEM-2001': return 'ERSEM (2001)'    		  	
  	if text ==  'ERSEM-2006': return 'ERSEM (2006)'    		
  	if text ==  'ERSEM-clim': return 'ERSEM'  	
  	if text ==  'ERSEM-clim_97-07': return 'ERSEM (\'97-\'07)'
  	if text ==  'ERSEM-2001': return 'ERSEM (2001)'    	  	
  	if text ==  'ERSEM-HighResp': return 'ERSEM (High Respiration)'  	  	   	  	 
  	if text ==  'Maredat': 	return 'Maredat'  	  	
  	if text ==  'MAREDAT': 	return 'Maredat'  	  	  	

  	if text ==  'Takahashi': 	return 'Takahashi 2009'  	  	  	
  	if text ==  'Seawifs': 		return 'Seawifs'  	  	
  	if text ==  'Seawifs-micro': 	return 'Seawifs Microphyto. chl.'  	  	
  	if text ==  'Seawifs-nano': 	return 'Seawifs Nanophyto. chl.'  	  	
  	if text ==  'Seawifs-pico': 	return 'Seawifs Picophyto. chl.'  
  	if text ==  'SeawifsBM-micro': 	return 'Seawifs Microphyto. Biomass'  	  	
  	if text ==  'SeawifsBM-nano': 	return 'Seawifs Nanophyto. Biomass'  	  	
  	if text ==  'SeawifsBM-pico': 	return 'Seawifs Picophyto. Biomass'    		  	
  	if text ==  'Seawifs-biomass': 	return 'Phytoplankton Biomass'
  	if text.lower() ==  'intpp': 	return 'Integrated Primary Production'  	  	  	
  	if text ==  'intppSurface': 	return 'Integrated PP'  	
  	if text ==  'PPint': 	return 'Integrated PP'  	
  	if text ==  'ppint': 	return 'Integrated PP'  	
   	if text ==  'PP': 	return 'MareDat PP'

  	if text ==  'medusa_1998': 	return 'MEDUSA (1998)'  
  	if text ==  'MEDUSA': 		return 'MEDUSA'    
  	if text == 'ZMI':		return 'Microzooplankton'
  	if text == 'ZME':		return 'Mesozooplankton'
  	if text == 'PHD':		return 'Diatoms'
  	if text in  ['CHL',]: 		return 'Chlorophyll'  

  	if text ==  'InitialConditions': return 'Initial Conditions'    
  	  	
  	if text in ['Months','months']:return 'Months'	 	  	  	  	
  	
  	if text == 'HighLatWinter':	return 'High Latitude Winter'
  	
  	if text in month_name: return text
  	for m in month_name:
  		t =  text.find(m)
  		if t>0: return getLongName(text[:t]) +' '+m
  	
  	if len(text) ==4:
  		try:	
  			year = int(text)
  			return text
  		except: pass
  	#if text in ['picophyto','microzoo','mesozoo','diatoms', 'bac', ]:
  	#	print "need to add ",text,"to get longname"
  	print "getLongName:\tERROR:\tCould not find Longname for ",text
	return text
	#assert False


	
		
	
	
def fancyUnits(units,debug=False):
	"""	Converts ascii units string into latex style formatting.
	"""
	units = units.replace('[','').replace(']','')
		
  	#if units in ['mg C/m^3','mg C/m^2',]:		return 'mg C m'+r'$^{-3}$'
  	if units in ['umol/l, uM, mo/l, ug/l, ',]:	return 'mg m'+r'$^{-3}$' # silly nitrates multi units
  	if units in ['mg C/m^3',]:			return 'mg C m'+r'$^{-3}$'
  	if units in ['mg Chl/m3','ng/L','mgCh/m3',]:		return 'mg Chl m'+r'$^{-3}$'  	
  	if units in ['mg C/m^3/d',]:			return 'mg C m'+r'$^{-3}$/day'
  	if units in ['mg N/m^3',]:			return 'mg N m'+r'$^{-3}$'  
  	if units in ['mg P/m^3',]:			return 'mg P m'+r'$^{-3}$'
  	if units in ['mmol N/m^3', 'mmol-N/m3' ]: 	return 'mmol N m'+r'$^{-3}$'
  	if units in ['mmol P/m^3', ]: 			return 'mmol P m'+r'$^{-3}$'
  	if units in ['mmol C/m^3', ]: 			return 'mmol C m'+r'$^{-3}$'
  	if units in ['umol F/m^3',]:			return r'$\mu$'+'mol m'+r'$^{-3}$'
  	if units in ['umol /m^3','umol / m3',]:		return r'$\mu$'+'mol m'+r'$^{-3}$' 
  	if units in ['mmol S/m^3', ]: 			return 'mmol S m'+r'$^{-3}$'  	
  	if units in ['mmolSi/m3', 'mmol Si/m^3', ]: 	return 'mmol Si m'+r'$^{-3}$'  	  	
  	if units in ['mmolFe/m3',]:			return 'mmol Fe m'+r'$^{-3}$'  	
  	
	if units in ['ug/l','mg/m^3','ug/L',]:  	return 'mg m'+r'$^{-3}$'
	if units in ['10^12 g Carbon year^-1',]:	return r'$10^{12}$'+' g Carbon/year'
	if units in ['mol C/m^',]:			return 'mol C/m'+r'$^{2}$'
  	if units in ['mmmol/m^3', 'mmol/m^3','umol/l','micromoles/l','mmolO2/m3']:
  							return 'mmol m'+r'$^{-3}$'
	if units in ['mmol/m^2']:			return 'mmol m'+r'$^{-2}$' 
	#if units in ['mmol/m^3']:			return 'mmol m'+r'$^{-3}$' 	
	if units in ['degrees Celsius', 'degreesC', 'C', 'degC', 'degrees_celsius',]:
							return r'$\,^{\circ}\mathrm{C}$'
	if units in ['psu','PSU',]:			return 'psu'
	#if units in ['umol/l',]:			return r'$\mu$'+'mol/l'
	if units in ['m','meters','meter',]:		return 'm'	
	if units in ['1/m',]:				return r'$\mathrm{m}^{-1}$'
	if units in ['m/s',]:				return r'$\mathrm{ms}^{-1}$'	
	#if units in ['ug/l']:			#	return 'mg m'+r'$^{-3}$'
	if units in ['W/m^2']:				return 'W m'+r'$^{-2}$'
	if units in ['umol/kg',]:			return r'$\mu$'+'mol kg'+r'$^{-1}$'
	if units in ['nmol/kg',]:			return 'nmol kg'+r'$^{-1}$'
	if units in ['tons C/d',]:			return 'tons C/day'
	if units in ['ug/L/d','ug                  ']:	return 'mg m'+r'$^{-3}$'+'/day'	#yes, there are lots of spaces
	if units.replace(' ','') in ['ug',]:		return r'$\mu$'+'g' #r'$\mu$'+	
	if units in ['1',]:			
		print 'fancyUnits:\tWarning:\tStrange units:',units
		return ''
	if units in ['uatm',]:				return r'$\mu$'+'atm'
	if units in ['ppmv',]:				return 'ppm'	
	if units in ['milliliters_per_liter',]:		return 'ml/l'
	print 'fancyUnits:\tERROR:\t',units,' not found in fancyUnits.'
	if debug:
		assert False
	return units 		
		
	
