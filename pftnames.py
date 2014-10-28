from calendar import month_name
from UKESMpython import AutoVivification,AutoVivToYaml,folder,YamlToDict
from itertoolsmodule import product
"""	This is a dictionary of all the terms that you'd need to pick out.
	pftnames['Model Name']['Functional type']['currency']

"""



            
            
pftnames = AutoVivification()
pftnames['ERSEM']['diatoms']['c'] = 'P1c'
pftnames['ERSEM']['diatoms']['n'] = 'P1n'
pftnames['ERSEM']['diatoms']['p'] = 'P1p'
pftnames['ERSEM']['diatoms']['f'] = 'P1f'
pftnames['ERSEM']['diatoms']['s'] = 'P1s'
pftnames['ERSEM']['diatoms']['chl'] = 'Chl1'

pftnames['ERSEM' ]['total']['chl'] = 'chl'
pftnames['MEDUSA']['total']['chl'] = 'CHL'

pftnames['MEDUSA']['diatoms']['n'] = ['']
pftnames['MEDUSA']['diatoms']['f'] = ['']
pftnames['MEDUSA']['diatoms']['s'] = ['']


# overloading different names with the same dict.
pftnames['E'] = pftnames['ERSEM']
pftnames['e'] = pftnames['ERSEM']
pftnames['ersem'] = pftnames['ERSEM']
pftnames['M'] = pftnames['MEDUSA']
pftnames['m'] = pftnames['MEDUSA']
pftnames['medusa'] = pftnames['MEDUSA']


regions = ['Surface','200m','100m','500m','1000m','Transect','All','',]
MaredatTypes = ['chl','diatoms','bac','mesozoo','picophyto','microzoo']

#WOATypes = [a+b for a,b in product(['silicate','nitrate','phosphate','salinity','temperature',],['Surface','500m','100m','200m','1000m','','Transect',''])]

Ocean_names	=['SouthPacificOcean',  'ArcticOcean',  'AntarcticOcean','NorthAtlanticOcean','SouthAtlanticOcean', 'NorthPacificOcean','IndianOcean',]

MLDTypes = ['mld','mld_DT02','mld_DR003','mld_DReqDTm02', ]
WOATypes = [a+b for a,b in product(['silicate','nitrate','phosphate','salinity','temperature',],regions)]

GEOTRACESTypes = ['iron',]

def getLongName(text):
	if type(text) in [type(['a','b',]),type(('a','b',))]:
		out = ''
		for t in text:out+=getLongName(t)+' '
		return out

			 		
  	if text == 'tempTransect':	return "Pacific Transect Temperature"
  	if text == 'tempSurface':	return "Surface Temperature"
  	if text == 'tempAll':		return "Temperature"  	
  	if text == 'temp100m':		return "Temperature (100m deep)"  	
  	if text == 'temp200m':		return "Temperature (200m deep)"  
  	if text == 'temp500m':		return "Temperature (500m deep)"  	
  	if text == 'temp1000m':		return "Temperature (1000m deep)"  	  		  	  	
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
  	if text == 'salTransect':	return "Pacific Transect Salinity"
  	if text == 'salSurface':	return "Surface Salinity"
  	if text == 'salAll':		return "Salinity"  	
  	if text == 'sal100m':		return "Salinity (100m deep)"  	
  	if text == 'sal200m':		return "Salinity (200m deep)"  	  	  	
  	if text == 'sal500m':		return "Salinity (500m deep)"  	
  	if text == 'sal1000m':		return "Salinity (1000m deep)"  	  	  	
  	
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


	
def getmt(): # Match Type

	"""
	getmt(): Get Match Type. 
	Typical usage:
		mt[ModelName or Data source][varaible name] = [list of variable names as they appears in the netcdf]
		
	mt always includes the coordinate names as their appear in the netcdf file too. ie:
		mt['GEOTRACES']['lat'] 		= 'Latitude'
	This assumes that all netcdf files from the same source kept the same coordinate names.
		
	However, it is possible tha you want to plot a combination of variables, such the sum of a series of values or the ratio of two.
	In some cases you may want to divide a value by 1000. 
	These cases work like this:
		mt[ModelName or Data source][varaible name]['name'] = variable name 
		mt[ModelName or Data source][varaible name][function in "excatractData()"] = [list of varaibles to combine in extractData]
		mt[ModelName or Data source][varaible name][new units] = 'mg m^-3' or whatever. (preferable from fancyUnits, below)
	"""
	
	#####
	# Models:
	mt = AutoVivification() # match type
	mt['ERSEM']['bac'] 		= ['B1c',]
	mt['ERSEM']['mesozoo'] 		= ['Z4c',]
	mt['ERSEM']['diatoms'] 		= ['P1c',]
	mt['ERSEM']['picophyto'] 	= ['P3c',]
	mt['ERSEM']['microzoo'] 	= ['Z5c',]			 	
	mt['ERSEM']['pCO2'] 		= ['pCO2w',]#'fAirSeaC',]
	mt['ERSEM']['chl'] 		= ['chl',]
	mt['ERSEM']['iron'] 		= ['N7f',]
	mt['ERSEM']['t']		= 'index_t'	
	mt['ERSEM']['z'] 		= 'deptht'
	mt['ERSEM']['lat'] 		= 'nav_lat'
	mt['ERSEM']['lon'] 		= 'nav_lon'
	mt['ERSEM']['cal'] 		= '365_day'

	for d in regions:#['Surface','200m','100m','500m','1000m','Transect','All',]:
	    for t in ['temperature','temp',]:
	    	mt['WOA'][t+d] 		= ['t_an',]#'t_mn',
	    	mt['NEMO'][t+d] 	= ['votemper',]

	    for t in ['salinity','sal',]:
	    	mt['WOA' ][t+d] 	= ['s_an',]#'s_mn',
	    	mt['NEMO'][t+d] 	= ['vosaline',]

	    for t in ['nitrate','nit']:
	    	mt['WOA'   ][t+d] 	= ['n_an',]#'n_mn',
	    	mt['MEDUSA'][t+d]  	= ['DIN',]
	    	mt['ERSEM' ][t+d]  	= ['N3n',]	    	
	    for t in ['silicate',]:
	    	mt['WOA'   ][t+d] 	= ['i_an',]#'i_mn',
	    	mt['MEDUSA'][t+d]  	= ['SIL',]
	    	mt['ERSEM' ][t+d]  	= ['N5s',]	
	    for t in ['phosphate',]:
	    	mt['WOA'   ][t+d] 	= ['p_an',]#'p_mn',
	    	mt['ERSEM' ][t+d]  	= ['N1p',]	
	
	mt['WOA']['t'] 			= 'index_t'
	mt['WOA']['z'] 			= 'depth'
	mt['WOA']['lat'] 		= 'lat'
	mt['WOA']['lon'] 		= 'lon'
	mt['WOA']['cal'] 		= 'standard'    		    	
    	
	    		    	
				
	mt['NEMO']['mld'] 		=  ['somxl010',]
	mt['NEMO']['mld_DT02'] 		=  ['somxl010',]
	mt['NEMO']['mld_DR003'] 	=  ['somxl010',]
	mt['NEMO']['mld_DReqDTm02'] 	=  ['somxl010',]		
	mt['NEMO']['t'] 		= 'index_t'	
	mt['NEMO']['z'] 		= 'deptht'
	mt['NEMO']['lat'] 		= 'nav_lat'
	mt['NEMO']['lon'] 		= 'nav_lon'
	mt['NEMO']['cal']		= '365_day'
	
	mt['MEDUSA']['iron']['mul1000']		=  ['FER',]
	mt['MEDUSA']['iron']['name']		=  'FER'
	mt['MEDUSA']['iron']['units']		=  'umol F/m^3'		
	mt['MEDUSA']['chl'] 			=  ['CHL',]	
	mt['MEDUSA']['diatoms']['N2Biomass'] 	=  ['PHD',]		
	mt['MEDUSA']['diatoms']['name'] 	=  'PHD'
	mt['MEDUSA']['diatoms']['units'] 	=  'mg C/m^3'	
	mt['MEDUSA']['mesozoo']['N2Biomass'] 	=  ['ZME',]			
	mt['MEDUSA']['mesozoo']['name'] 	=  'ZME'
	mt['MEDUSA']['mesozoo']['units'] 	=  'mg C/m^3'		
	mt['MEDUSA']['microzoo']['N2Biomass'] 	=  ['ZMI',]				
	mt['MEDUSA']['microzoo']['name'] 	=  'ZMI'
	mt['MEDUSA']['microzoo']['units'] 	=  'mg C/m^3'		
	mt['MEDUSA']['t'] 			= 'index_t'	
	mt['MEDUSA']['z'] 			= 'deptht'
	mt['MEDUSA']['lat'] 			= 'nav_lat'
	mt['MEDUSA']['lon'] 			= 'nav_lon'
	mt['MEDUSA']['cal'] 			= '365_day'
	
	mt['Medusa'] = mt['MEDUSA']
	
	#####
	# Data:
	mt['Maredat']['bac'] 		= ['BIOMASS',]
	mt['Maredat']['mesozoo'] 	= ['BIOMASS',]
	mt['Maredat']['diatoms'] 	= ['BIOMASS',]
	mt['Maredat']['picophyto'] 	= ['BIOMASS',]
	mt['Maredat']['microzoo'] 	= ['BIOMASS',]
	mt['Maredat']['PP'] 		= ['PP',]
	mt['Maredat']['chl']['name']	= 'Chlorophylla'
	mt['Maredat']['chl']['div1000']	= ['Chlorophylla',]
	mt['Maredat']['chl']['units']	= ['ug/L',]
	mt['Maredat']['t'] 		= 'index_t'
	mt['Maredat']['z'] 		= 'DEPTH'
	mt['Maredat']['lat'] 		= 'LATITUDE'
	mt['Maredat']['lon'] 		= 'LONGITUDE'
	mt['Maredat']['cal'] 		= 'standard'	
	mt['MAREDAT'] 			= mt['Maredat']
			
		 
	mt['Takahashi']['pCO2'] 	= ['PCO2_SW',]#'DELTA_PCO2',]	'TFLUXSW06',
	mt['Takahashi']['t'] 		= 'index_t'
	mt['Takahashi']['z'] 		= 'index_z'
	mt['Takahashi']['lat'] 		= 'LAT'
	mt['Takahashi']['lon'] 		= 'LON'
	mt['Takahashi']['cal'] 		= 'standard'	
	
	#mt['intPP']['intPP']		= ['PPint',]								
	mt['GEOTRACES']['iron']		= ['Fe_D_CONC_BOTTLE',]#'Fe_D_CONC_BOTTLE_FIA','Fe_S_CONC_BOTTLE',]
	mt['GEOTRACES']['t']		= 'MONTH'
	mt['GEOTRACES']['z'] 		= 'DEPTH'
	mt['GEOTRACES']['lat'] 		= 'Latitude'
	mt['GEOTRACES']['lon'] 		= 'Longitude'
	mt['GEOTRACES']['cal'] 		= 'standard'
		
	mt['IFREMER']['mld']		= ['mld',]	
	mt['IFREMER']['mld_DT02']	= ['mld',]
	mt['IFREMER']['mld_DR003']	= ['mld',]
	mt['IFREMER']['mld_DReqDTm02']	= ['mld',]		
	mt['IFREMER']['t'] 		= 'index_t'
	mt['IFREMER']['z'] 		= 'index_z'
	mt['IFREMER']['lat'] 		= 'lat'
	mt['IFREMER']['lon'] 		= 'lon'
	mt['IFREMER']['cal'] 		= 'standard'	
		
	#mt['PP']['PP'] 		= ['PP',]
	AutoVivToYaml(mt,folder('yaml')+'matchMetadata.yaml')
	
	mt = 0
	print mt
	mt = YamlToDict(folder('yaml')+'matchMetadata.yaml',)
	return mt	
	
def fancyUnits(units,debug=False):#'mg C/m^2',
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
	#if units in ['umol/l',]:return r'$\mu$'+'mol/l'
	if units in ['m','meters','meter',]:		return 'm'	
	if units in ['1/m',]:				return r'$\mathrm{m}^{-1}$'
	#if units in ['ug/l']:			#	return 'mg m'+r'$^{-3}$'
	if units in ['W/m^2']:				return 'W m'+r'$^{-2}$'
	if units in ['umol/kg',]:			return r'$\mu$'+'mol kg'+r'$^{-1}$'
	if units in ['nmol/kg',]:			return 'nmol kg'+r'$^{-1}$'
	if units in ['tons C/d',]:			return 'tons C/day'
	if units in ['ug/L/d','ug                  ']:	return 'mg m'+r'$^{-3}$'+'/day'#yes, there are lots of spaces
	if units.replace(' ','') in ['ug',]:		return r'$\mu$'+'g' #r'$\mu$'+	
	if units in ['1',]:			
		print 'fancyUnits:\tWarning:\tStrange units:',units
		return ''
	if units in ['uatm',]:				return r'$\mu$'+'atm'
	print 'fancyUnits:\tERROR:\t',units,' not found in fancyUnits.'
	if debug:assert False
	return units 		
		
	
