#!/usr/bin/ipython 
#Standard Python modules:
from sys import argv
from os.path import exists
#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, NestedDict
from p2p import matchDataAndModel,makePlots,makeTargets

from pftnames import getmt

###	Potential problems?
###		Reliance on ORCA1 grid
###		Can't take data from more than one file at a time.


def main():

    if len(argv[1:]): models  = argv[1:]
    else:models = ['MEDUSA','ERSEM','NEMO']
    print 'Models:',models
    for model in models:
	ERSEMjobID = 'xhonc'
	
	if model == 'NEMO':	year='1893'
	if model == 'ERSEM':	year='1893'		
	if model == 'MEDUSA':	year='1998'
	plotallcuts = False
	av = AutoVivification()
	
	mt = getmt()
	
	MAREDATFolder 	= "/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"
	WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"	
	MEDUSAFolder	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
	ERSEMFolder	= "/data/euryale7/scratch/ledm/UKESM/ERSEM/"
	
	workingDir = folder("/data/euryale7/scratch/ledm/ukesm_postProcessed/"+model+'-'+year)
	
	av['chl']['Data'  ]['File'] 		= MAREDATFolder+"MarEDat20121001Pigments.nc"	
	av['chl']['MEDUSA']['File'] 		= MEDUSAFolder+"medusa_bio_"+year+".nc"	
	av['chl']['ERSEM' ]['File'] 		= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMMisc.nc'			
	av['chl']['Data']['Vars'] 		= ['Chlorophylla',]
	av['chl']['MEDUSA']['Vars'] 		= ['CHL',]	
	av['chl']['ERSEM']['Vars'] 		= ['chl',]
	av['chl']['region'] 			= ''
	
	av['diatoms']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20120716Diatoms.nc"	
	av['diatoms']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"	
	av['diatoms']['ERSEM' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMphytoBm.nc'				
	av['diatoms']['Data']['Vars'] 		= ['BIOMASS',]
	av['diatoms']['MEDUSA']['Vars'] 	= ['PHD',]	
	av['diatoms']['ERSEM']['Vars'] 		= ['P1c',]
	av['diatoms']['region'] 		= ''	
	

	av['bac']['Data'  ]['File'] 		= MAREDATFolder+"MarEDat20120214Bacteria.nc"	
	av['bac']['ERSEM' ]['File'] 		= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMbac.nc'			
	av['bac']['Data']['Vars'] 		= ['BIOMASS',]
	av['bac']['ERSEM']['Vars'] 		= ['B1c',]
	av['bac']['region'] 			= ''	
	
	av['picophyto']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20111206Picophytoplankton.nc"	
	av['picophyto']['ERSEM' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMphytoBm.nc'			
	av['picophyto']['Data']['Vars'] 	= ['BIOMASS',]
	av['picophyto']['ERSEM']['Vars'] 	= ['P3c',]
	av['picophyto']['region'] 		= ''	
		

	av['microzoo']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20120424Microzooplankton.nc"	
	av['microzoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"	
	av['microzoo']['ERSEM' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMzoo.nc'			
	av['microzoo']['Data']['Vars'] 		= ['BIOMASS',]
	av['microzoo']['MEDUSA']['Vars'] 	= ['ZMI',]	
	av['microzoo']['ERSEM']['Vars'] 	= ['Z5c',]
	av['microzoo']['region'] 		= ''	
	
	
	av['mesozoo']['Data'  ]['File'] 	= MAREDATFolder+"MarEDat20120705Mesozooplankton.nc"	
	av['mesozoo']['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"	
	av['mesozoo']['ERSEM' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMzoo.nc'			
	av['mesozoo']['Data']['Vars'] 		= ['BIOMASS',]
	av['mesozoo']['MEDUSA']['Vars'] 	= ['ZME',]	
	av['mesozoo']['ERSEM']['Vars'] 		= ['Z4c',]
	av['mesozoo']['region'] 		= ''
							
			

	for woa in ['salinity','temperature',]:
		if model == 'MEDUSA' or model == 'ERSEM' : continue
		if woa == 'salinity':		NEMOVars  	= ['vosaline',]
		if woa == 'temperature':	NEMOVars  	= ['votemper',]
		for s in ['Surface','500m','100m','200m','1000m',]:
			av[woa+s]['Data'  ]['File'] 	= WOAFolder+woa+'_monthly_1deg.nc'	
			av[woa+s]['NEMO' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_NEMO.nc'	
			av[woa+s]['Data']['Vars'] 	= [woa[0]+'_mn',woa[0]+'_an','depth','lat','lon','time'] 
			av[woa+s]['NEMO']['Vars'] 	= NEMOVars
			av[woa+s]['region'] 		= s	    
	
	av['mld']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
	av['mld']['NEMO' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_NEMO.nc'			
	av['mld']['Data']['Vars'] 	= ['mld','mask',]
	av['mld']['NEMO']['Vars'] 	= ['somxl010',]	
	av['mld']['region'] 		= ''
	

	av['mld_DR003']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
	av['mld_DR003']['NEMO' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_NEMO.nc'			
	av['mld_DR003']['Data']['Vars'] 	= ['mld','mask',]
	av['mld_DR003']['NEMO']['Vars'] 	= ['somxl010',]	
	av['mld_DR003']['region'] 		= ''

	av['mld_DReqDTm02']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
	av['mld_DReqDTm02']['NEMO' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_NEMO.nc'			
	av['mld_DReqDTm02']['Data']['Vars'] 	= ['mld','mask',]
	av['mld_DReqDTm02']['NEMO']['Vars'] 	= ['somxl010',]	
	av['mld_DReqDTm02']['region'] 		= ''

	#av['mld_DT02']['Data'  ]['File'] 	= "/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc"
	#av['mld_DT02']['NEMO' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_NEMO.nc'			
	#av['mld_DT02']['Data']['Vars'] 	= ['mld','mask',]
	#av['mld_DT02']['NEMO']['Vars'] 	= ['somxl010',]	
	#av['mld_DT02']['region'] 		= ''
				
	#    mldvars = ['somxl010',]
	#    for datafile in ["/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc",
	#   			"/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DR003_c1m_reg2.0.nc",
	#    			"/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DReqDTm02_c1m_reg2.0.nc", ]:
	    				
	

	for woa in ['silicate','nitrate','phosphate',]:
		if model == 'MEDUSA' and woa == 'phosphate' : continue

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
			MEDVars 	= []
		
		for s in ['Surface','100m','200m','500m',]:#'Transect',]:#'All',
			av[woa+s]['Data'  ]['File'] 	= WOAFolder+woa+'_monthly_1deg.nc'	
			av[woa+s]['ERSEM' ]['File'] 	= ERSEMFolder + ERSEMjobID+'/'+year+'/'+ERSEMjobID+'_'+year+'_ERSEMNuts.nc'	
			av[woa+s]['MEDUSA']['File'] 	= MEDUSAFolder+"medusa_bio_"+year+".nc"									
			av[woa+s]['Data']['Vars'] 	= [l+'_mn',l+'_an','depth','lat','lon','time'] 
			av[woa+s]['ERSEM']['Vars'] 	= ERSEMVars
			av[woa+s]['MEDUSA']['Vars'] 	= MEDVars				
			av[woa+s]['region'] 		= s	    


	shelves = []
	avshelves = AutoVivification()
	
	#for name in sorted(['nitrateSurface','phosphateSurface','silicateSurface',]):
	for name in sorted(av.keys()):
		if not av[name][model]: continue
		
		if not exists(av[name]['Data']['File']):
			print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name]['Data']['File']
			continue
		if not exists(av[name][model]['File']):
			print "testsuite_p2p.py:\tWARNING:\tFile does not exist", av[name][model]['File']
			continue
		print "\n\n\ntestsuite_p2p.py:\tINFO:\tRunning:",name
		b = matchDataAndModel.matchDataAndModel(av[name]['Data']['File'], 
							av[name][model]['File'],
							name,
							DataVars  = av[name]['Data']['Vars'],
							ModelVars = av[name][model]['Vars'],
							jobID=model,
							year=year,
							workingDir = folder(workingDir+name),
							region = av[name]['region'])
		m = makePlots.makePlots(	b.MatchedDataFile, 
						b.MatchedModelFile, 
						name, 
						model, 
						year = year, 
						plotallcuts=plotallcuts, 
						workingDir = folder(workingDir+name),
						compareCoords=True)
		#print 'OutPutShelves:', m.shelves
		
		shelves.extend(m.shelves)
		filename = folder('images/'+model+'/Targets/')+model+'_'+year+'_'+name+'.png'
		
		t = makeTargets.makeTargets(	m.shelves, 
						filename,
						#name=name,
						legendKeys = ['newSlice','ykey',],
						debug=True)
						#imageDir='', 
						
		for ns, sh in zip(m.newSlices,m.shelves):
			region = av[name]['region']
			avshelves[name.replace(region,'')][region][ns] = sh
	
		
	#####
	# All possible targets in  this model
	
	filename = folder('images/'+model+'/Targets/')+model+'_'+year+'_everything.png'		
	makeTargets.makeTargets(	shelves, 
					filename,
					legendKeys = ['name', 'newSlice','ykey',],					
					debug=True)#imageDir='', diagramTypes=['Target',]			
	
	#avshelves[model][name]
	
	
	#print "shelves:",shelves
	print "Working dir:",workingDir
	







	
if __name__=="__main__":
	main() 
	print 'The end.'
	
	
