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

from calendar import month_name
from UKESMpython import AutoVivification,AutoVivToYaml,folder,YamlToDict
from itertools import product
from os.path import exists
import numpy as np



#####
#	

regions 	= ['Surface','200m','100m','500m','1000m','Transect','All','',]

MaredatTypes 	= ['chl','diatoms','bac','mesozoo','picophyto','microzoo']

Ocean_names	= ['SouthPacificOcean',  'ArcticOcean',  'AntarcticOcean','NorthAtlanticOcean','SouthAtlanticOcean', 'NorthPacificOcean','IndianOcean',]

IFREMERTypes 	= ['mld','mld_DT02','mld_DR003','mld_DReqDTm02', ]

WOATypes 	= ['silicate','nitrate','phosphate','salinity','temperature',]

TAKAHASHITypes 	= ['pCO2',]

GEOTRACESTypes 	= ['iron',]


#####
# Get Match Type:
#	
def getmt(loadYaml=False): 
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
	def N2Biomass(nc,keys):	return nc.variables[keys[0]][:]* 79.573
	def mul1000(nc,keys):	return nc.variables[keys[0]][:]* 1000.
	def div1000(nc,keys):	return nc.variables[keys[0]][:]/ 1000.	
	def applymask(nc,keys):	return np.ma.masked_where(nc.variables[keys[1]][:]==0.,nc.variables[keys[0]][:])
	 #np.ma.masked_where(nc.variables[keys[1]][:],nc.variables[keys[0]][:])
					
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
	mt['ERSEM']['iron'] 			= ['N7f',]
	mt['ERSEM']['t']			= 'index_t'	
	mt['ERSEM']['z'] 			= 'deptht'
	mt['ERSEM']['lat'] 			= 'nav_lat'
	mt['ERSEM']['lon'] 			= 'nav_lon'
	mt['ERSEM']['cal'] 			= '365_day'
	    		    	
	mt['NEMO']['temperature'] 		= ['votemper',]	
 	mt['NEMO']['salinity'] 			= ['vosaline',]				
	mt['NEMO']['mld'] 			= ['somxl010',]
	mt['NEMO']['mld_DT02'] 			= ['somxl010',]
	mt['NEMO']['mld_DR003'] 		= ['somxl010',]
	mt['NEMO']['mld_DReqDTm02'] 		= ['somxl010',]		
	mt['NEMO']['t'] 			= 'index_t'	
	mt['NEMO']['z'] 			= 'deptht'
	mt['NEMO']['lat'] 			= 'nav_lat'
	mt['NEMO']['lon'] 			= 'nav_lon'
	mt['NEMO']['cal']			= '365_day'


	
	mt['MEDUSA']['chl'] 			=  ['CHL',]	
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

    	mt['MEDUSA']['nitrate'] 	 	= ['DIN',]
    	mt['MEDUSA']['pCO2']	  		= ['OCN_PCO2',]			    	
    	mt['MEDUSA']['silicate']	  	= ['SIL',]			
	mt['MEDUSA']['t'] 			= 'index_t'	
	mt['MEDUSA']['z'] 			= 'deptht'
	mt['MEDUSA']['lat'] 			= 'nav_lat'
	mt['MEDUSA']['lon'] 			= 'nav_lon'
	mt['MEDUSA']['cal'] 			= '365_day'
	#mt['Medusa']				= mt['MEDUSA']
	
	#####
	# Data:
	mt['MAREDAT']['bac'] 			= ['BIOMASS',]
	mt['MAREDAT']['mesozoo'] 		= ['BIOMASS',]
	mt['MAREDAT']['diatoms'] 		= ['BIOMASS',]
	mt['MAREDAT']['picophyto'] 		= ['BIOMASS',]
	mt['MAREDAT']['microzoo'] 		= ['BIOMASS',]
	mt['MAREDAT']['PP'] 			= ['PP',]
	mt['MAREDAT']['chl']['name']		= 'Chlorophylla'
	mt['MAREDAT']['chl']['vars']		= ['Chlorophylla',]
	mt['MAREDAT']['chl']['convert']		= div1000	
	mt['MAREDAT']['chl']['units']		= ['ug/L',]
	mt['MAREDAT']['t'] 			= 'index_t'
	mt['MAREDAT']['z'] 			= 'DEPTH'
	mt['MAREDAT']['lat'] 			= 'LATITUDE'
	mt['MAREDAT']['lon'] 			= 'LONGITUDE'
	mt['MAREDAT']['cal'] 			= 'standard'	
	#mt['Maredat'] 				= mt['MAREDAT']
			
	mt['WOA']['temperature'] 		= ['t_an',]#'t_mn',
  	mt['WOA']['salinity'] 			= ['s_an',]#'s_mn',
  	mt['WOA']['nitrate'] 			= ['n_an',]#'s_mn',  	
	mt['WOA']['silicate'] 			= ['i_an',]#'i_mn',
	mt['WOA']['phosphate'] 			= ['p_an',]#'p_mn',	    	  		
	mt['WOA']['t'] 				= 'index_t'
	mt['WOA']['z'] 				= 'depth'
	mt['WOA']['lat'] 			= 'lat'
	mt['WOA']['lon'] 			= 'lon'
	mt['WOA']['cal'] 			= 'standard'    
			 	
	mt['TAKAHASHI']['pCO2'] 		= ['PCO2_SW',]#'DELTA_PCO2',]	'TFLUXSW06',
	mt['TAKAHASHI']['t'] 			= 'TIME'
	mt['TAKAHASHI']['z'] 			= 'index_z'
	mt['TAKAHASHI']['lat'] 			= 'LAT'
	mt['TAKAHASHI']['lon'] 			= 'LON'
	mt['TAKAHASHI']['cal'] 			= 'standard'
	#mt['intPP']['intPP']			= ['PPint',]	
								
	mt['GEOTRACES']['iron']			= ['Fe_D_CONC_BOTTLE',]#'Fe_D_CONC_BOTTLE_FIA','Fe_S_CONC_BOTTLE',]
	mt['GEOTRACES']['t']			= 'MONTH'
	mt['GEOTRACES']['z'] 			= 'DEPTH'
	mt['GEOTRACES']['lat'] 			= 'Latitude'
	mt['GEOTRACES']['lon'] 			= 'Longitude'
	mt['GEOTRACES']['cal'] 			= 'standard'
		
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

def getLongName(text):
	print "Getting long name:",text
	if type(text) in [type(['a','b',]),type(('a','b',))]:
		return ' '.join([getLongName(t) for t in text])
		#out = ''
		#for t in text:out+=getLongName(t)+' '
		#return out
	
	mt = getmt()
	
	noChange = ['Surface',]
	noChange.extend(mt.keys())
	
	if text in noChange:return text
	
	firstLetterCaps = ['temperature', "salinity", "nitrate", "phosphate" ,'silicate','iron',]
	
	if text.lower() in firstLetterCaps:	return text.title()


	
  	if text == 'Transect':		return "Pacific Transect"
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
	  		  	  	
  	if text == 'temperatureTransect':	return "Pacific Transect Temperature"
  	if text == 'temperatureSurface':	return "Surface Temperature"
  	if text == 'temperatureAll':		return "Temperature"  	
  	if text == 'temperature100m':		return "Temperature (100m deep)"  	
  	if text == 'temperature200m':		return "Temperature (200m deep)"  
  	if text == 'temperature500m':		return "Temperature (500m deep)"  	
  	if text == 'temperature1000m':		return "Temperature (1000m deep)"  	  		  	  	
  	
  	if text == 'salinityTransect':		return "Pacific Transect Salinity"
  	if text == 'salinitySurface':		return "Surface Salinity"
  	if text == 'salinityAll':		return "Salinity"  	
  	if text == 'salinity100m':		return "Salinity (100m deep)"  	
  	if text == 'salinity200m':		return "Salinity (200m deep)"  	  	  	
  	if text == 'salinity500m':		return "Salinity (500m deep)"  	
  	if text == 'salinity1000m':		return "Salinity (1000m deep)"  	 	  	
  	
  	if text == 'nitrateTransect':	return "Pacific Transect Nitrate (WOA14)"
  	if text == 'nitrateSurface':	return "Surface Nitrate (WOA14)"
  	if text == 'nitrateAll':	return "Nitrate (WOA14)"  	
  	if text == 'nitrate100m':	return "Nitrate (100m deep)"  	  	
  	if text == 'nitrate200m':	return "Nitrate (200m deep)"  
  	if text == 'nitrate500m':	return "Nitrate (500m deep)"  	  	
  	  		  	  	
  	if text == 'phosphateTransect':	return "Pacific Transect Phosphate"
  	if text == 'phosphateSurface':	return "Surface Phosphate"
  	if text == 'phosphateAll':	return "Phosphate"  	
  	if text == 'phosphate100m':	return "Phosphate (100m deep)"  	  	
  	if text == 'phosphate200m':	return "Phosphate (200m deep)"    	
  	if text == 'phosphate500m':	return "Phosphate (500m deep)"    	  	
  	if text == 'silicateTransect':	return "Pacific Transect Silicate"
  	if text == 'silicateSurface':	return "Surface Silicate"
  	if text == 'silicateAll':	return "Silicate"  	
  	if text == 'silicate100m':	return "Silicate (100m deep)"  	  	
  	if text == 'silicate200m':	return "Silicate (200m deep)"    	
  	if text == 'silicate500m':	return "Silicate (500m deep)"    	
  	  	
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
  	if text ==  'i_mn': 	return 'Mean Silicate'  	
  	if text ==  'i_an': 	return 'Silicate'  
  	  	  	  	  	
  	
     	if text == 'IFREMER':	return "IFREMER"
  	if text == 'mld':		return 'Mixed Layer Depth'     	
  	if text == 'mld_DT02':		return 'MLD: Fixed Threshold Temperature '
  	if text == 'mld_DR003':		return 'MLD: Fixed Threshold Density'
  	if text == 'mld_DReqDTm02':	return 'MLD: Variable Threshold Density'  	  	
 
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

	if text =='Overestimate_2sig': return "Overestimate 2 sigma"
	if text =='Overestimate_3sig': return "Overestimate 3 sigma"
	if text =='Underestimate_2sig': return "Overestimate 2 sigma"
	if text =='Underestimate_3sig': return "Overestimate 3 sigma"			
	   	
  	if text in ['Tropics','Temperate','Arctic','Surface','Depth', 'SalArtifact','OffAxis','pCO2',]:return text
	if text in ['Overestimate','Underestimate','Matched', 'Arctic','Tropics','Temperate', 'Surface',]:return text
	
  	if text == 'SurfaceNoArtics':	return ""
	
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


  	if text == 'ArcticOcean':	return "Arctic Ocean"	  		  	
  	if text == 'AntarcticOcean':	return "Antarctic Ocean" 		  	
  	if text == 'NorthAtlanticOcean':return "North Atlantic Ocean"
  	if text == 'SouthAtlanticOcean':return "South Atlantic Ocean"	  		  	
  	if text == 'NorthPacificOcean':	return "North Pacific Ocean"		  	
  	if text == 'SouthPacificOcean':	return "South Pacific Ocean"	  		  	
  	if text == 'IndianOcean':	return "Indian Ocean"
  	if text == 'ignoreExtraArtics':	return "No Arctic Oceans (50 degrees)"  	
  	if text == 'ignoreMoreArtics':	return "No Arctic Oceans (60 degrees)"
  	if text == 'ignoreMidArtics':	return "No Arctic Oceans (65 degrees)"  
  	if text == 'ignoreArtics':	return "No Arctic Oceans (70 degrees)" 
  	
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
  		
	if text in ['PCO2_SW', 'pCO2']:	return 'pCO2'

  	if text ==  'NEMO': 	return 'NEMO'  	  	
  	if text ==  'Nemo': 	return 'NEMO'
  	
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

  	if text ==  'Takahashi': 	return 'Takahashi 2009'  	  	  	
  	if text ==  'Seawifs': 		return 'Seawifs'  	  	
  	if text ==  'Seawifs-micro': 	return 'Seawifs Microphyto. chl.'  	  	
  	if text ==  'Seawifs-nano': 	return 'Seawifs Nanophyto. chl.'  	  	
  	if text ==  'Seawifs-pico': 	return 'Seawifs Picophyto. chl.'  
  	if text ==  'SeawifsBM-micro': 	return 'Seawifs Microphyto. Biomass'  	  	
  	if text ==  'SeawifsBM-nano': 	return 'Seawifs Nanophyto. Biomass'  	  	
  	if text ==  'SeawifsBM-pico': 	return 'Seawifs Picophyto. Biomass'    		  	
  	if text ==  'Seawifs-biomass': 	return 'Phytoplankton Biomass'
  	if text ==  'intPP': 	return 'WOA Integrated PP'  	  	  	
   	if text ==  'PP': 	return 'MareDat PP'

  	if text ==  'medusa_1998': 	return 'MEDUSA (1998)'  
  	if text ==  'MEDUSA': 		return 'MEDUSA'    
  	if text == 'ZMI':		return 'Microzooplankton'
  	if text == 'ZME':		return 'Mesozooplankton'
  	if text == 'PHD':		return 'Diatoms'
  	if text in  ['CHL',]: 		return 'Chlorophyll'  

  	if text ==  'InitialConditions': return 'Initial Conditions'    
  	  		 	  	  	  	
  	if text in month_name: return text
  	#if text in ['picophyto','microzoo','mesozoo','diatoms', 'bac', ]:
  	#	print "need to add ",text,"to get longname"
  	print "getLongName:\tERROR:\tCould not find Longname for ",text

	assert False


	
		
	
	
def fancyUnits(units,debug=False):
	"""	Converts ascii units string into latex style formatting.
	"""
  	#if units in ['mg C/m^3','mg C/m^2',]:		return 'mg C m'+r'$^{-3}$'
  	if units in ['umol/l, uM, mo/l, ug/l, ',]:	return 'mg m'+r'$^{-3}$' # silly nitrates multi units
  	if units in ['mg C/m^3',]:			return 'mg C m'+r'$^{-3}$'
  	if units in ['mg Chl/m3','ng/L',]:		return 'mg Chl m'+r'$^{-3}$'  	
  	if units in ['mg C/m^3/d',]:			return 'mg C m'+r'$^{-3}$/day'
  	if units in ['mg N/m^3',]:			return 'mg N m'+r'$^{-3}$'  
  	if units in ['mg P/m^3',]:			return 'mg P m'+r'$^{-3}$'
  	if units in ['mmol N/m^3', 'mmol-N/m3' ]: 	return 'mmol N m'+r'$^{-3}$'
  	if units in ['mmol P/m^3', ]: 			return 'mmol P m'+r'$^{-3}$'
  	if units in ['mmol C/m^3', ]: 			return 'mmol C m'+r'$^{-3}$'
  	if units in ['umol F/m^3',]:			return r'$\mu$'+'mol m'+r'$^{-3}$'
  	if units in ['mmol S/m^3', ]: 			return 'mmol S m'+r'$^{-3}$'  	
  	if units in ['mmolSi/m3', 'mmol Si/m^3', ]: 	return 'mmol Si m'+r'$^{-3}$'  	  	
  	if units in ['mmolFe/m3',]:			return 'mmol Fe m'+r'$^{-3}$'  	
  	
	if units in ['ug/l','mg/m^3','ug/L',]:  	return 'mg m'+r'$^{-3}$'
	if units in ['10^12 g Carbon year^-1',]:	return r'$10^{12}$'+' g Carbon/year'
	if units in ['mol C/m^',]:			return 'mol C/m'+r'$^{2}$'
  	if units in ['mmmol/m^3', 'mmol/m^3','umol/l','micromoles/l',]:
  							return 'mmol m'+r'$^{-3}$'
	if units in ['mmol/m^2']:			return 'mmol m'+r'$^{-2}$' 
	#if units in ['mmol/m^3']:			return 'mmol m'+r'$^{-3}$' 	
	if units in ['degrees Celsius', 'degreesC', 'C',]:
							return r'$\,^{\circ}\mathrm{C}$'
	if units in ['psu','PSU',]:			return 'psu'
	#if units in ['umol/l',]:			return r'$\mu$'+'mol/l'
	if units in ['m','meters','meter',]:		return 'm'	
	if units in ['1/m',]:				return r'$\mathrm{m}^{-1}$'
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
	print 'fancyUnits:\tERROR:\t',units,' not found in fancyUnits.'
	if debug:
		assert False
	return units 		
		
	
