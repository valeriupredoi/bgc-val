#!/usr/bin/ipython
from netCDF4 import num2date
from datetime import datetime
from sys import argv
#from re import sub
from os.path import exists,split, getmtime, basename
from glob import glob
from shelve import open as shOpen
#from xtraPlots import hovmoeller
from matplotlib.colors import LogNorm
from matplotlib import pyplot, ticker
#from numpy import array,ones
#from numpy.ma import  min as mmin, max as mmax, np.ma.masked_where, array as ma,zeros as mzeros, masked,divide,nonzero,array as marray,masked_outside,masked_inside,clip, power
#from numpy import std, isnan as npNaN, isinf as npInf, array,arange, delete, tile,NaN,mean,zeros,abs, ma,random
from calendar import month_name
from itertoolsmodule import product
from scipy.stats import linregress
from StatsDiagram import StatsDiagram
from scipy.stats.mstats import scoreatpercentile
#from BGCnames import fancyUnits

import numpy as np 

#from deMoraTools import getLogTicks,folder,mnStr
from ncdfView import ncdfView
#from iMarNetPython import noLeapConvert, ukp.shouldIMakeFile,depthKeys,alwaysInclude,getFileList,getJobID,getMareDatatype,stripmask,getLongName,getDataOrigin
#from iMarNetPython import getOrcaIndexCC,getORCAdepth
#from iMarNetPython import makeLatSafe, ukp.makeLonSafe,histPlot, ukp.makeLonSafeArr
#from iMarNetPython import cumulativeSumPlot, scatterPlot, robinPlotPair, robinPlot


import UKESMpython as ukp 
#from UKESMpython import ukp.shouldIMakeFile, ukp.folder,robinPlotPair,histPlot,scatterPlot,getOrcaIndexCC,makeLatSafe
from pftnames import getLongName, AutoVivification, getkd, getmt

#from wavelets import wavelet

    

"""	This code makes matched plots, hexbins, scatter plots, and so on.

"""

BioLogScales 	= ['bac','mesozoo','diatoms','picophyto','microzoo','PP','Seawifs', 'iron'] 	
noXYLogs 	= [ 'pCO2',
		#'nitrateSurface', 	'nitrateAll',	'nitrateTransect',
		#'phosphateSurface',	'phosphateAll',	'phosphateTransect',
		'silicateSurface',	'silicateAll',	'silicateTransect', 'silicate100m','silicate200m','silicate500m',
		'tempSurface',		'tempAll',	'tempTransect',	'temp100m',	'temp200m','temp1000m',	'temp500m',
		'salSurface', 		'salAll',	'salTransect',	'sal100m',	'sal200m','sal1000m',	'sal500m',]

class makePlots:
  def __init__(self,matchedDataFile,matchedModelFile,name, jobID='MEDUSA',year='clim',region='', compareCoords=True,saveShelve=True,plotallcuts=  False): #xfilename,yfilename,

  	
  	#if not self.setUpFromName(name,jobID): return
  	self.xfn =matchedModelFile
  	self.yfn =matchedDataFile  	
    	self.name = name
  	self.xtype = jobID
	if self.name == 'chl': self.ytype = 'Maredat'

	self.compType = self.xtype+self.ytype+region

  	self.saveShelve = saveShelve
	#self.getplotsInfoName()
  	self.kd = getkd()
  	self.BioLogScales = BioLogScales
  	self.xnc = ncdfView(self.xfn,Quiet=True)
  	self.ync = ncdfView(self.yfn,Quiet=True)

	if compareCoords: self.CompareCoords()
	
	self.basicCut	  = 'All'#'Standard'##,
	self.standardCuts = ['5-95pc','ignoreInlandSeas','OffShelf','ignoreExtraArtics','aboveZero','0.1','0.2','0.01']	
	self.newSlices 		=['All','Standard',
				#'aboveZero', 'TypicalIron','maskBelowBathy',
				#'ignoreMoreArtics','ignoreArtics', 'ignoreMidArtics','ignoreExtraArtics'
				#'OffShelf','OnShelf', '5-95pc','maskBelowBathy', 
				#'Shallow','Depth','nonZero','aboveZero'
				#'NitArtifact','SalArtifact','0-99pc',
				]
	self.months = {month_name[i+1]:i for i in xrange(0,12) }
	self.depthRanges	=['OffShelf','maskBelowBathy', 'OnShelf',] # 'Depth_0-10m','Depth_10-20m','Depth_20-50m','Depth_50-100m','Depth_100-500m','Depth_500m',
	self.percentiles	=['0-1pc','1-5pc','5-25pc','25-40pc','40-60pc','60-75pc','75-95pc','95-99pc','99-100pc',]
	self.latregions		=['NorthTemperate','SouthTemperate','NorthTropics','Equatorial',
				'SouthTropics','Antarctic','NorthArctic',]#'Arctic','Tropics','Temperate']
	self.Seas		=['ignoreMediteranean','BlackSea','ignoreBlackSea','RedSea','BalticSea','PersianGulf', 'ignoreInlandSeas',]	# 'ignoreRedSea', 'ignoreBalticSea','ignorePersianGulf',]
	self.Oceans		=['SouthPacificOcean',  'ignoreMoreArtics','ArcticOcean','AntarcticOcean','NorthAtlanticOcean',
				 'SouthAtlanticOcean','NorthPacificOcean','IndianOcean',] #'ignoreExtraArtics','ignoreMidArtics','ignoreArtics',

	self.QualityCuts 	=['Overestimate','Underestimate','Overestimate_2sig','Underestimate_2sig',
			    	  'Overestimate_3sig','Underestimate_3sig', 'Matched','OffAxis','1-99pc','5-95pc','0-99pc',]
	self.OceanMonths	= sorted([i for i in product(self.Oceans,self.months)] )
	self.OceanMonths.extend(sorted([i for i in product(['All',],self.months)]))
	
	if plotallcuts:
	 self.plotMonths	= 0#True
	 self.plotdepthRanges	=0
	 self.plotpercentiles	=0#True	
	 self.plotLatRegions	=0
	 self.plotQualityCuts	=0#True	
	 self.plotSeas		=0#True		 
	 self.plotOceans	= 0#True	
	 self.plotOceanMonths   = True	
	else: 	
	 self.plotMonths	=0#True
	 self.plotdepthRanges	=0#True	
	 self.plotpercentiles	=0#True	
	 self.plotLatRegions	=0#True
	 self.plotQualityCuts	=0#True
	 self.plotSeas		=0#True		 
	 self.plotOceans	=0#True	
	 self.plotOceanMonths   = 0	 	 	 

	
	if self.plotMonths: 	 self.newSlices.extend(self.months.keys())
	if self.plotdepthRanges: self.newSlices.extend(self.depthRanges)
	if self.plotpercentiles: self.newSlices.extend(self.percentiles)
	if self.plotLatRegions:	 self.newSlices.extend(self.latregions)	
	if self.plotQualityCuts: self.newSlices.extend(self.QualityCuts)		
	if self.plotSeas: 	 self.newSlices.extend(self.Seas)			
	if self.plotOceans: 	 self.newSlices.extend(self.Oceans)
	if self.plotOceanMonths: self.newSlices.extend(self.OceanMonths)
	print "New SLICES:", 	 self.newSlices
	self.plotWithSlices()
	    	
  	self.xnc.close()
  	self.ync.close()  	
  	
  def setUpFromName(self,name,jobID):
  	print "setUpFromName:\t",name,"\t",jobID
    	self.name = name
    	self.jobID = jobID
    	if jobID.find('Restart')>-1:
	    	self.compType = 'MaredatMatched-'+jobID#.replace('_','-')
    		self.model = 'InitialConditions'	    	
	elif jobID.lower().find('medusa')>-1:
	    	self.compType = 'MaredatMatched-'+jobID.replace('_','-')
    		self.model = 'MEDUSA'
	else:
	    	self.compType = 'MaredatMatched-'+jobID.replace('_','-')
    		self.model = 'Ersem'
    		
	 	 
	
	jobID = jobID.replace('_','*' )	
	#Matched_NetCDF_Folder = "/users/modellers/ledm/workspace/iMarNetp2p/outNetCDF/"+self.compType+"/"+self.name+"/"
	Matched_NetCDF_Folder = "/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/"+self.compType+"/"+self.name+"/"
	if self.name in ['salSurface','tempSurface',
			 'salTransect','tempTransect',
			 'sal100m','temp100m',
			 'sal200m','temp200m',			 			 
			 'sal500m','temp500m',			 			 
			 'sal1000m','temp1000m',			 			 			 
			 'salAll','tempAll',]:
	   	mfiles =Matched_NetCDF_Folder +self.name[:3]+"*_monthly_1deg*_2.nc"
	   	efiles =Matched_NetCDF_Folder +jobID+"_SalTemp*_"+self.name[:3]+"*_1D.nc"	   
	elif self.name in ['nitrateSurface', 'phosphateSurface', 'silicateSurface',
			   'nitrateAll',  'phosphateAll','silicateAll',
			   'nitrateTransect', 'phosphateTransect',  'silicateTransect',
			   'nitrate100m',  'phosphate100m','silicate100m',
			   'nitrate200m',  'phosphate200m','silicate200m',
			   'nitrate500m',  'phosphate500m','silicate500m',
			   ]:
	   	mfiles =Matched_NetCDF_Folder+"*_monthly_1deg*_2.nc"
	   	if self.model == 'InitialConditions': 
	   		efiles =Matched_NetCDF_Folder +"*"+self.name+"_1D.nc"	   
	   	else:	efiles =Matched_NetCDF_Folder +jobID+"_*"+self.name+"*_1D.nc"	   	   	
	elif self.name in ['pCO2',]:
	   	mfiles =Matched_NetCDF_Folder+"takahashi2009_month_flux_pCO2_2006c_noHead_2.nc"
	   	efiles =Matched_NetCDF_Folder +jobID+"_Diag*_"+self.name[:3]+"*_1D.nc"  
	elif self.name in ['iron',]:
	   	mfiles =Matched_NetCDF_Folder+"Iron_GEOTRACES*_2.nc"
	   	if self.model == 'InitialConditions': 	   	
	   		efiles =Matched_NetCDF_Folder +"*iron_1D.nc"  
		else:   efiles =Matched_NetCDF_Folder +jobID+"_*_1D.nc"  
	elif self.name in ['chl','intPP','nitSurface','nitTransect','PP', 'mld_DT02','mld_DR003','mld_DReqDTm02']:
	   	mfiles =Matched_NetCDF_Folder + "*_2.nc"
	   	efiles =Matched_NetCDF_Folder + jobID+"_*_1D.nc"
	elif self.name in ['Seawifs',]:
	   	mfiles =Matched_NetCDF_Folder + "*_2.nc"
	   	efiles =Matched_NetCDF_Folder + jobID+"*Diag*_1D.nc"
	elif self.name in ['Seawifs-micro', 'Seawifs-nano','Seawifs-pico',
				'SeawifsBM-micro', 'SeawifsBM-nano','SeawifsBM-pico',
				'Seawifs-biomass']:
		mfiles = "/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/"+self.compType+"/Seawifs/*_2.nc"
	   	efiles = "/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/"+self.compType+"/Seawifs/*ERSEM*_1D.nc"
		#mfiles = "/users/modellers/ledm/workspace/iMarNetp2p/outNetCDF/"+self.compType+"/Seawifs/*_2.nc"
	   	#efiles = "/users/modellers/ledm/workspace/iMarNetp2p/outNetCDF/"+self.compType+"/Seawifs/*ERSEM*_1D.nc"
	   		   	
	else:
	   	mfiles =Matched_NetCDF_Folder + "/Mar*_2.nc"
	   	efiles =Matched_NetCDF_Folder + jobID+"_*"+self.name+"_1D.nc"

	marefiles 	= glob(mfiles)
	ERSEMfiles 	= glob(efiles)

	harshKill = False
	if not len(marefiles):
		print "could not find MAREDAT file:",self.name,'in', mfiles
		if harshKill:assert False
		else: return False
		
	if not len(ERSEMfiles):
		print "could not find Model file:",self.name	,'in',efiles		
		if harshKill:assert False
		else: return False
		
	if len(marefiles)>1:
		print "ERROR: Found Too many MAREDAT files:"+self.name+'in'+ mfiles+':\n\t',marefiles
		if harshKill:assert False
		else: return False
	
	if len(ERSEMfiles)>1:
		print "ERROR: Found Too many ERSEM files:"+self.name+'in'+ efiles+':\n\t',ERSEMfiles
		if harshKill:assert False
		else: return False

	if not exists(marefiles[0]): 
		print "maredatFile does not exist:",marefiles[0]
		if harshKill:assert False
		else: return False
		
	if not exists(ERSEMfiles[0]): 
		print "ERSEMfile does not exist:",ERSEMfiles[0]
		if harshKill:assert False
		else: return False
	print "setUpFromName:\txfile:\t", ERSEMfiles[0]			
	print "setUpFromName:\tyfile:\t", marefiles[0]
	
  	self.xfn = ERSEMfiles[0]
  	self.yfn = marefiles[0]

	xfn = basename(self.xfn).lower()
	yfn = basename(self.yfn).lower()	

	self.xtype=getDataOrigin(xfn)
	self.ytype=getDataOrigin(yfn)

		
	try: print 'setUpFromName:\tINFO:','Found:',self.xtype
	except:	self.xtype = self.xfn.split('/')[-2]
	try: print 'setUpFromName:\tINFO:','Found:',self.ytype
	except:	self.ytype = self.yfn.split('/')[-2]
	print 'setUpFromName:\tINFO:','Found:',self.xtype,self.ytype
	
	if self.xtype ==  self.ytype: self.plotType = self.xtype
	if self.xtype !=  self.ytype: self.plotType = self.xtype+self.ytype
		
	print 'setUpFromName:\tINFO:',"found Types:",self.plotType, '\tx:',self.xtype,'\ty:',self.ytype,'\tname:',self.name
	return True
	

  def CompareCoords(self,):
	xkeys,ykeys = ['index_t',] , ['index_t',]
	for k in ['t','lat','lon','z','lon',]:	xkeys.append(self.kd[self.xtype][k])
	for k in ['t','lat','lon','z','lat',]:	ykeys.append(self.kd[self.ytype][k])
		
	  	 	  	
  	for xkey,ykey in zip(xkeys,ykeys):
	    	if xkey not in self.xnc.variables.keys():continue  	    
	    	if ykey not in self.ync.variables.keys():continue
		filename = ukp.folder(['images',self.compType,'CompareCoords',self.name])+self.name+xkey+'vs'+ykey+'.png'	    	
		if not ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):continue
		print "CompareCoords:\tx:",xkey, "\ty:",ykey
		if xkey not in self.xnc.variables.keys():
			print xkey, "not in xnc"
			assert False
		if ykey not in self.ync.variables.keys():
			print ykey, "not in ync"
			assert False		
		
		dxtmp = np.ma.array(self.xnc(xkey)[:])
		dytmp = np.ma.array(self.ync(ykey)[:])
		print "CompareCoords:\t",xkey,':', len(dxtmp),"\tx min:", dxtmp.min(),'x max:',dxtmp.max()
		print "CompareCoords:\t",ykey,':', len(dytmp),'\ty min:', dytmp.min(),'y max:',dytmp.max()
		dx = np.ma.masked_where(dxtmp.mask + dytmp.mask, dxtmp).compressed()
		dy = np.ma.masked_where(dxtmp.mask + dytmp.mask, dytmp).compressed()
		print "CompareCoords:\t",xkey,':', len(dx),"\t:",ykey,':',len(dy), dx.min(),dx.max(), dy.min(), dy.max()		
		del dxtmp
		del dytmp


		fig = pyplot.figure()
		fig.set_size_inches(8, 12)		
		ax = pyplot.subplot(411)

		rects1 = pyplot.hist((dx,dy),label=[xkey,ykey],histtype='bar',bins=72/2)#,alpha=0.2)
		pyplot.legend()

		ax.set_title(xkey + ' and '+ykey)		
		ax = pyplot.subplot(412)		
		rects3 = pyplot.hist(dx - dy,bins=72,label=[xkey + ' - '+ykey])
		pyplot.legend()

		ax = pyplot.subplot(212)
		pyplot.hexbin(dx, dy, bins='log',gridsize = 72, cmap=pyplot.get_cmap('gist_yarg'),mincnt=0)#extent=plotrange,
		cb = pyplot.colorbar()				
		mmax = max(dx.max(),dy.max())
		mmin = min(dx.min(),dy.min())
			
		fx = np.arange(mmin, mmax, (mmax-mmin)/20.)
		pyplot.plot(fx,fx, 'k--')
		ax.set_aspect("equal")

		pyplot.xlabel(self.xtype+' '+xkey)
		pyplot.ylabel(self.ytype+' '+ykey)

		print "\tSaving: " + filename
		pyplot.savefig(filename,dpi=100,)
		pyplot.close()	  	

  	
 		
  def getFileName(self,newSlice,xkey,ykey):
	  file_prefix = ukp.folder(['images',self.compType,'Slices-'+self.basicCut,self.name,])

	  file_suffix = '_'+self.xtype+'.png'
	  if newSlice in self.months.keys():
		filename = ukp.folder([file_prefix,'months'])+self.name+'_'+mnStr(self.months[newSlice])+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	  elif newSlice in self.depthRanges:
		filename = ukp.folder([file_prefix,'DepthRanges'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix	  
	  elif newSlice in self.percentiles:
		filename = ukp.folder([file_prefix,'Percentiles'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix	  
	  elif newSlice in self.latregions:
		filename = ukp.folder([file_prefix,'LatRegions'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix	 
	  elif newSlice in self.QualityCuts:
		filename = ukp.folder([file_prefix,'QualityCuts'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	  elif newSlice in self.Seas:
		filename = ukp.folder([file_prefix,'Seas'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	  elif newSlice in self.Oceans:
		filename = ukp.folder([file_prefix,'Oceans'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	  elif newSlice in self.OceanMonths:
		if type(newSlice) in [type(['a','b',]),type(('a','b',))]:
		  	print 'getFileName:', newSlice,
		  	newSlice = ''.join(newSlice)
		  	print '-->',newSlice	  
		filename = ukp.folder([file_prefix,'OceanMonths'])+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	  else:
		filename = file_prefix+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	  return filename	    



  def makeMask(self, newSlice, xt,xz,xy,xx,xd):	  
	print "makePlots:\tmakeMask:\tinitialise:\t",self.name, '\t',newSlice
		
	if newSlice in ['OffAxis', 'Overestimate','Underestimate','Matched', 
			'Overestimate_2sig','Underestimate_2sig', 
			'Overestimate_3sig','Underestimate_3sig',]:
		print "makeMask:\tSlice", newSlice, "requires both datasets, and you should never see this"
		assert False

  	if newSlice == 'All': return np.zeros(len(xd))		
			  	
	nmask = np.ones(len(xd))	# everything masked	
	nmask = np.zeros(len(xd))	# nothing masked	

 	if newSlice in ['maskBelowBathy', 'OnShelf','OffShelf',]:
		bathync = ncdfView("/users/modellers/ledm/workspace/iMarNet/outNetCDF/ORCA1bathy.nc",Quiet=True)
		bathy = abs(bathync("bathymetry")[:])
		latcc, loncc =  bathync("lat")[:], bathync("lon")[:]	
		bathync.close()
		shelfDepth=500.
		shelveFn = ukp.folder("shelves/MaredatErsemMatch-latlon/"+self.name)+"diag_"+self.name+".shelve"
		try:
			s = shOpen(shelveFn)		
			lldict  = s['lldict']
			s.close()
		except:	lldict={}	
	
		print "Bathy mask: before mask:", newSlice, nmask.sum(), 'of', len(nmask)			  	
		for i,z in enumerate(xz):
			try:
				la,lo = lldict[(xy[i],xx[i])]
			except:
				la,lo = ukp.getOrcaIndexCC(xy[i],xx[i],latcc,loncc,debug=False)
				lldict[(xy[i],xx[i])] = la,lo
			if la==lo==-1:
				print "Corner case:", la,lo,bathy[la,lo] 
				nmask[i]=1
			if newSlice == "maskBelowBathy":
				if (bathy[la,lo]-10.) > abs(z): nmask[i]=1	
			if newSlice == "OnShelf":
				if  bathy[la,lo] >= shelfDepth: nmask[i]=1	
			if newSlice == "OffShelf":
				if  bathy[la,lo] < shelfDepth:  nmask[i]=1
		
			if i%100000==0 or i==(len(xz)+1):
				s = shOpen(shelveFn)		
				s['lldict'] = lldict 
				s.close()			
		print "Bathy mask:", newSlice, nmask.sum(), 'of', len(nmask)
		return nmask
	

	if newSlice in ['1-99pc','5-95pc','0-99pc'] or newSlice in self.percentiles:	  		
		if newSlice in self.percentiles:
		  	tmp  = newSlice.replace('pc','').split('-')
		  	pcmin,pcmax = float(tmp[0]),float(tmp[1])
		  	print newSlice, pcmin, pcmax
		  	if pcmin == 0:	ymin = yd.min()
			else:		ymin = scoreatpercentile(yd,pcmin)
		  	if pcmax == 100:ymax = yd.max()
			else:		ymax = scoreatpercentile(yd,pcmax) 	
		if newSlice in ['1-99pc',]:
			ymin = scoreatpercentile(xd,1)
			ymax = scoreatpercentile(xd,99)
		if newSlice in ['5-95pc',]:
			ymin = scoreatpercentile(xd,5)
			ymax = scoreatpercentile(xd,95)
		
		if newSlice in ['0-99pc',]:
			ymin = xd.min()
			ymax = scoreatpercentile(xd,99)	
		print  "makeMask:\t",newSlice,ymin,ymax
		return  np.ma.masked_outside(xd,ymin, ymax).mask		

	if newSlice in self.months.keys():
		print "masking a month:",newSlice,xt[0], xt[-1]
		return np.ma.masked_where( xt != self.months[newSlice],nmask).mask 
	
	if newSlice == "0.1":	return np.ma.masked_where( xd==0.1, xd).mask
	if newSlice == "0.2":	return np.ma.masked_where( xd==0.2, xd).mask
	if newSlice == "0.01":	return np.ma.masked_where( xd==0.01, xd).mask	
	
	
	if newSlice == 'Shallow':	return np.ma.masked_where( xz > 200.,nmask).mask
	if newSlice == 'Depth':		return np.ma.masked_where( xz < 200.,nmask).mask
	if newSlice == 'Zoom':		return np.ma.masked_where( xd > 10.,nmask).mask 
	if newSlice == 'Zoom5':		return np.ma.masked_where( xd > 5., nmask).mask 
	if newSlice == 'Zoom2':		return np.ma.masked_where( xd > 2., nmask).mask 
	if newSlice == 'nonZero':	return np.ma.masked_where( xd == 0.,nmask).mask 			
	if newSlice == 'aboveZero':	return np.ma.masked_where( xd <= 0.,nmask).mask
	if newSlice == 'Tropics':	return np.ma.masked_where( abs(xy) >23.,nmask).mask 			
	if newSlice == 'Equatorial':	return np.ma.masked_where( abs(xy) >7.,nmask).mask 
	if newSlice == 'Temperate':	return np.ma.masked_where( (abs(xy) <23.)+(abs(xy) >60.),nmask).mask 	
	if newSlice == 'NorthTropics':	return np.ma.masked_where( (xy >23.)+(xy < 7.),nmask).mask 			
	if newSlice == 'SouthTropics':	return np.ma.masked_where( (xy <-23.)+(xy > -7.),nmask).mask 				
	if newSlice == 'NorthTemperate':return np.ma.masked_where( (xy <23.)+(xy >60.),nmask).mask 			
	if newSlice == 'SouthTemperate':return np.ma.masked_where( (xy >-23.)+(xy <-60.),nmask).mask 	
	if newSlice == 'Arctic':	return np.ma.masked_where( abs(xy) < 60.,nmask).mask
	if newSlice == 'Antarctic':	return np.ma.masked_where( xy > -60.,nmask).mask 			
	if newSlice == 'NorthArctic':	return np.ma.masked_where( xy < 60.,nmask).mask 														
	if newSlice == 'SalArtifact': 	return np.ma.masked_where( (xd > 15.)+(xd < 10.),nmask).mask 
	if newSlice == 'NitArtifact':	return np.ma.masked_where( (xd > 6.)+(xd < 2.),  nmask).mask 
	if newSlice == 'Depth_0-10m': 	return np.ma.masked_where( abs(xz) > 10.,nmask).mask 
	if newSlice == 'Depth_10-20m': 	return np.ma.masked_where( (abs(xz) < 10.)+(abs(xz) > 20.),nmask).mask 
	if newSlice == 'Depth_20-50m': 	return np.ma.masked_where( (abs(xz) > 20.)+(abs(xz) > 50.),nmask).mask 
	if newSlice == 'Depth_50-100m': return np.ma.masked_where( (abs(xz) < 50.)+(abs(xz) > 100.),nmask).mask
	if newSlice == 'Depth_100-500m':return np.ma.masked_where( (abs(xz) < 100.)+(abs(xz) > 500.),nmask).mask
	if newSlice == 'Depth_500m': 	return np.ma.masked_where(  abs(xz) < 500.,nmask).mask	

	if newSlice == 'TypicalIron': 	return np.ma.masked_where( (xd<=0.) *(xd<=4.),nmask).mask 
	
	if newSlice == 'BlackSea': 	
		mx = np.ma.masked_outside(xx, 25.9,41.7).mask
		my = np.ma.masked_outside(xy, 39.8,48.1).mask				
		return np.ma.masked_where( mx+my,nmask).mask 
		
	if newSlice == 'ignoreBlackSea':
		mx = np.ma.masked_inside(xx, 25.9,41.7).mask
		my = np.ma.masked_inside(xy, 39.8,48.1).mask				
		return np.ma.masked_where( mx*my,nmask).mask 	
		
	if newSlice == 'BalticSea': 	
		mx = np.ma.masked_outside(xx, 12.5,30.7).mask
		my = np.ma.masked_outside(xy, 53.0,66.4).mask				
		return np.ma.masked_where( mx+my,nmask).mask 
		
	if newSlice == 'ignoreBalticSea':
		mx = np.ma.masked_inside(xx, 12.5,30.7).mask
		my = np.ma.masked_inside(xy, 53.0,66.4).mask				
		return np.ma.masked_where( mx*my,nmask).mask 	

	if newSlice == 'RedSea': 	
		mx = np.ma.masked_outside(xx, 30.0,43.0).mask
		my = np.ma.masked_outside(xy, 12.4,30.4).mask				
		return np.ma.masked_where( mx+my,nmask).mask 
		
	if newSlice == 'ignoreRedSea':
		mx = np.ma.masked_inside(xx, 30.0,43.0).mask
		my = np.ma.masked_inside(xy, 12.4,30.4).mask				
		return np.ma.masked_where( mx*my,nmask).mask 	
		
	if newSlice == 'PersianGulf': 	
		mx = np.ma.masked_outside(xx, 47.5, 56.8).mask
		my = np.ma.masked_outside(xy, 22.3, 32.1).mask				
		return np.ma.masked_where( mx+my,nmask).mask 
		
	if newSlice == 'ignorePersianGulf':
		mx = np.ma.masked_inside(xx, 47.5, 56.8).mask
		my = np.ma.masked_inside(xy, 22.3, 32.1).mask				
		return np.ma.masked_where( mx*my,nmask).mask 										
		
	if newSlice == 'ignoreMediteranean':
		mx  = np.ma.masked_inside(xx, -5.8, 42.5).mask #E
		my  = np.ma.masked_inside(xy, 30., 43.).mask	#N			
		mx2 = np.ma.masked_inside(xx, 0., 20.).mask #E
		my2 = np.ma.masked_inside(xy, 32., 47.).mask #N
		m = mx*my+ mx2*my2
		return np.ma.masked_where( m,nmask).mask 		

	if newSlice in ['ignoreInlandSeas', 'IndianOcean']: 
		mx = np.ma.masked_inside(xx, 47.5,  56.8).mask * np.ma.masked_inside(xy, 22.3, 32.1).mask	
		mx += np.ma.masked_inside(xx, 30.0, 43.0).mask * np.ma.masked_inside(xy, 12.4,30.4).mask	
		mx += np.ma.masked_inside(xx, 12.5, 30.7).mask * np.ma.masked_inside(xy, 53.0,66.4).mask
		mx += np.ma.masked_inside(xx, 25.9, 41.7).mask * np.ma.masked_inside(xy, 39.8,48.1).mask		
		mx += np.ma.masked_inside(xx, -5.8, 42.5).mask * np.ma.masked_inside(xy, 30., 43.).mask
		mx += np.ma.masked_inside(xx, 0.0,  20.0).mask * np.ma.masked_inside(xy, 32., 47.).mask 
		if newSlice == 'ignoreInlandSeas':return np.ma.masked_where( mx,nmask).mask 		
		mx += np.ma.masked_outside(xx, 25.,100.).mask
		my = np.ma.masked_outside(xy, -50.,30.).mask
		if newSlice == 'IndianOcean':return np.ma.masked_where( mx+my,nmask).mask 

	if newSlice == 'AntarcticOcean': 	return np.ma.masked_where(  xy >-50.,nmask).mask 
	if newSlice == 'ArcticOcean': 		return np.ma.masked_where(  xy < 60.,nmask).mask 
	if newSlice == 'ignoreArtics':		return np.ma.masked_outside(xy,-70., 70.).mask
	if newSlice == 'ignoreMidArtics':	return np.ma.masked_outside(xy,-65., 65.).mask
	if newSlice == 'ignoreMoreArtics':	return np.ma.masked_outside(xy,-60., 60.).mask
	if newSlice == 'ignoreExtraArtics':	return np.ma.masked_outside(xy,-50., 50.).mask 
	if newSlice == 'NorthAtlanticOcean': 	return np.ma.masked_outside(ukp.makeLonSafeArr(xx), -80.,0.).mask + np.ma.masked_outside(xy, 10.,60.).mask
	if newSlice == 'SouthAtlanticOcean':	return np.ma.masked_outside(ukp.makeLonSafeArr(xx), -65.,20.).mask + np.ma.masked_outside(xy, -50.,-10.).mask

	if newSlice == 'NorthPacificOcean':
		mx = np.ma.masked_inside(xx,-100., 120. ).mask
		mx += np.ma.masked_inside(xx,260., 365. ).mask		
		mx += np.ma.masked_outside(xy,10., 60. ).mask
		return mx
	
	if newSlice == 'SouthPacificOcean': 	
		mx = np.ma.masked_inside(xx,-70., 140. ).mask
		mx += np.ma.masked_inside(xx,290., 365. ).mask		
		my = np.ma.masked_outside(xy,-10., -50. ).mask
		return np.ma.masked_where( mx+my,nmask).mask 
	print "Mask region not accepted:",newSlice
	assert False

  def extractData(self, nc, mt,key = ['',]):
  	""" The idea here is that the data in mt[type][name] can be a list, where the first value is the operation that you want to run.
  	    Also accepts the keywords: 'product', 'sum', 'mean', 'quotient'
  	"""
  	try: 
  		a = mt.keys()
  		print "extractData: MT is a dict", a

  	except:
  		print "extractData: mt Not a dict:", mt, key
  		return np.ma.array(nc(key)[:])
  	if mt[0] == {}:	print "warning:MT may be empty"

  	
  	print "extractData: Extracting data:\tinit:",mt
  		
  	#for name in mt[1:]:
  	#	if  name not  in nc.variables.keys(): 	
  	#		print "ExtractData:\tWarning:\t",name , ' not in nc (', mt,')'
  	#		assert False
  	
	if mt['sum']:
		print "Extracting data:\tSumming:", mt['sum']
  	  	xd = nc(mt['sum'][0])[:]
		for  name in mt['sum'][1:]:	xd +=nc(name)[:]
  	
  	if mt['product']:
  		print "Extracting data:\tmultiplying:", mt['product']
  	  	xd = nc(mt['product'][0])[:]
		for  name in mt['product'][1:]:	xd *= nc(name)[:]

  	if mt['productPC']:
  		print "Extracting data:\tmultiplying:", mt['productPC']
  	  	xd = nc(mt['productPC'][0])[:]*nc(mt['productPC'][1])[:]/100.   	  	

  	if mt['N2Biomass']:
  		
  		print "Extracting data:\tmultiplying:", mt['N2Biomass'] ,'by 79.573'
  	  	xd = nc(mt['N2Biomass'][0])[:] * 79.573
  	  		
	if mt['div1000']:
  		print "Extracting data:\tDividing by 1000. ", mt['div1000']
   	  	xd = nc(mt['div1000'][0])[:]/1000.   	  	 				

  	if mt['SWtoBmass']:
  		print "Extracting data:\tconverting seawifsPFT% into Biomass:", mt['SWtoBmass']
  	  	xd = nc(mt['SWtoBmass'][0])[:]*nc(mt['SWtoBmass'][1])[:]/100. 
		xd = 79. * power(xd, 0.65)
		#fit From http://www.int-res.com/articles/meps_oa/m383p073.pdf
		# doi: 10.3354/meps07998
	if mt['Chl2BM']:
  		print "Extracting data:\tconverting Chl to Biomass:", mt['Chl2BM']
  		xd = 79. * power(nc(mt['Chl2BM'][0])[:], 0.65)
  				
  	if mt['mean']:
  		print "Extracting data:\tmeaning:", mt['mean'] 
  	  	xd = nc(mt['mean'][0])[:]
		for  name in mt['mean'][1:]:	xd +=nc(name)[:]
		xd  = xd/float(len(mt['mean']))

	if mt['divide']:
  		print "Extracting data:\tdividing", mt['divide'] 	
		xd = nc(mt['divide'][0])[:]/nc(mt['divide'][1])[:]
		
  	return np.ma.array(xd)


  		


  def plotWithSlices(self):  
	print "plotWithSlices:\txtype:",self.xtype,"\tytype:",self.ytype,"\tname:",self.name
  	
	mt = getmt()

	#time and depth
	xt = self.xnc(self.kd[self.xtype]['t'])[:]
	yt = self.ync(self.kd[self.ytype]['t'])[:]
	xz = self.xnc(self.kd[self.xtype]['z'])[:]
	yz = self.ync(self.kd[self.ytype]['z'])[:]

	#lat and lon
	xy = self.xnc(self.kd[self.xtype]['lat'])[:]
	yy = self.ync(self.kd[self.ytype]['lat'])[:]
	xx = self.xnc(self.kd[self.xtype]['lon'])[:]
	yx = self.ync(self.kd[self.ytype]['lon'])[:]	

	#data
	try:	xd = self.xnc(mt[self.xtype][self.name][0])[:]
	except:	xd = self.extractData(self.xnc,mt[self.xtype][self.name])	
	try:	yd = self.ync(mt[self.ytype][self.name][0])[:]
	except:	yd = self.extractData(self.ync,mt[self.ytype][self.name])	
	print "plotWithSlices:\tx",xd.min(),xd.mean(),xd.max()
	print "plotWithSlices:\ty",yd.min(),yd.mean(),yd.max()
	
	if mt[self.ytype][self.name] == ['Chlorophylla',]:	yd = yd/1000.

	basicCutMask = self.makeMask('All',xt,xz,xy,xx,xd).astype(int)
        if self.name in ['mld_DT02','mld_DR003','mld_DReqDTm02']:
        	mldMask = self.ync('mask')[:]
        	mldMask = np.ma.masked_where(mldMask==0.,mldMask).mask
        	#basicCutMask += np.ma.masked_where(yd > 10.E6,yd).mask
        	        
	if self.basicCut in ['Standard',]:	  	  
	    print "Preparing Basic Cuts is Standard (5-95, no inland seas, off shelf)"	    
	    for ns in self.standardCuts: 
		if self.name in ['tempSurface','tempTransect', 'tempAll'] and ns in ['aboveZero',]:continue # no negative or zero values allowed.
		if self.name not in ['nitrateSurface','nitrateAll','nitrateTransect',] 		and ns in ['0.1','0.2']:continue	    			
		if self.name not in ['phosphateSurface','phosphateAll','phosphateTransect',] 	and ns in ['0.01',]:continue
		basicCutMask += self.makeMask(ns,xt,xz,xy,xx,xd).astype(int)
	   	basicCutMask += self.makeMask(ns,yt,yz,yy,yx,yd).astype(int)
	   
	
	for newSlice in self.newSlices:
	  xmask = self.makeMask('All',xt,xz,xy,xx,xd).astype(int)
	  ymask = self.makeMask('All',yt,yz,yy,yx,yd).astype(int)
	  	  	
	  if newSlice == "maskBelowBathy" and (self.name.lower().find('surface')>-1 or self.name in ['Seawifs',]):continue
	  if newSlice == 'Depth' and self.name.lower().find('surface')>-1: continue
	  if newSlice in  self.depthRanges and self.name.lower().find('surface')>-1: continue	 
	  if newSlice == 'SalArtifact' and self.name not in ['salTransect', 'salSurface']:continue
	  if newSlice == 'NitArtifact' and self.name not in ['nitTransect', 'nitSurface']:continue
	  

	  
	  xkeys = []
	  ykeys = []
	  plots = [] 	  
	  try: 		xkeys.append(mt[self.xtype][self.name]['name'])
	  except: 	xkeys = mt[self.xtype][self.name]
	  try: 		ykeys.append(mt[self.ytype][self.name]['name'])
	  except: 	ykeys = mt[self.ytype][self.name]  
	  print xkeys, ykeys
	  
	  for xk,yk in product(xkeys,ykeys):
	  	print 'listing plots:', xk, yk
		plots.append((xk,yk))
	  # this is a test to check if any files exist in the series
	  filename = self.getFileName(newSlice,xkeys[0],ykeys[0])
	  #if not ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):continue	  

	  if type(newSlice) in [type(['a','b',]),type(('a','b',))]:
	  	for n in newSlice:
	  	  xmask += self.makeMask(n,xt,xz,xy,xx,xd).astype(int)	  
	  	  ymask += self.makeMask(n,yt,yz,yy,yx,yd).astype(int)
		  
	  elif newSlice in  ['OffAxis', 'Overestimate','Underestimate','Matched', 'Overestimate_2sig','Underestimate_2sig', 'Overestimate_3sig','Underestimate_3sig',]:
	  	# These ones require both!
		  if newSlice in ['OffAxis', 'Overestimate','Underestimate','Matched']: distance = (std(xd)+std(yd))/2.
		  if newSlice in ['Overestimate_2sig','Underestimate_2sig',]:  	  	distance = std(xd)+std(yd)
		  if newSlice in ['Overestimate_3sig','Underestimate_3sig',]:  	  	distance = 3.*(std(xd)+std(yd))/2.
		  
	  	  for i in xrange(len(self.xnc('index')[:])):		  
	  	  	if newSlice in ['Overestimate','Overestimate_2sig','Overestimate_3sig'] :
				if abs(xd[i] - yd[i]) < distance: xmask[i] = 1
				if xd[i] < yd[i]: xmask[i] = 1
		
			if newSlice in ['Underestimate','Underestimate_2sig','Underestimate_3sig',]:
				if abs(xd[i] - yd[i]) < distance: xmask[i] = 1
				if xd[i] > yd[i]: xmask[i] = 1
	
			if newSlice == 'Matched':
				if abs(xd[i] - yd[i]) > distance: xmask[i] = 1
			
	  elif newSlice == 'Standard':
	  	if self.basicCut in ['Standard',]: 
		  	xmask = basicCutMask
		  	ymask = basicCutMask
		else:  
	 	    xmask = self.makeMask('All',xt,xz,xy,xx,xd).astype(int)
		    ymask = self.makeMask('All',yt,yz,yy,yx,yd).astype(int)	  	
	  	    for ns in self.standardCuts: #['5-95pc','ignoreInlandSeas','OffShelf',]:
			if self.name in ['tempSurface','tempTransect', 'tempAll'] and ns in ['aboveZero',]:continue # no negative or zero values allowed.	    	  	    
			if self.name not in ['nitrateSurface','nitrateAll','nitrateTransect',] 		and ns in ['0.1','0.2']:continue	    			
			if self.name not in ['phosphateSurface','phosphateAll','phosphateTransect',] 	and ns in ['0.01',]:continue	    						
	  		xmask += self.makeMask(ns,xt,xz,xy,xx,xd).astype(int)
	  	 	ymask += self.makeMask(ns,yt,yz,yy,yx,yd).astype(int)
	  	 	
	  elif newSlice in ['All','aboveZero',]:  	
	 	xmask = self.makeMask(newSlice,xt,xz,xy,xx,xd).astype(int)
		ymask = self.makeMask(newSlice,yt,yz,yy,yx,yd).astype(int)	  	
	  else:
	  	xmask += self.makeMask(newSlice,xt,xz,xy,xx,xd).astype(int)
	  	ymask += self.makeMask(newSlice,yt,yz,yy,yx,yd).astype(int)
	  	print xmask.sum(), ymask.sum() 

	  print self.name, newSlice
	  
          if self.name in ['mld_DT02','mld_DR003','mld_DReqDTm02']:
        	ymask += mldMask

	  keepMaskAsIs = ['All','TypicalIron',]
	  keepMaskAsIs.extend(self.QualityCuts)
	  keepMaskAsIs.extend(self.percentiles)	  
	  
	  nmask = (xmask + ymask).astype(bool)
	  	
	  if self.basicCut in ['Standard',] and newSlice not in keepMaskAsIs:	
		  nmask = (xmask + ymask + basicCutMask).astype(bool)  	# Apply cuts across thge board.
	  else:
	  	  nmask = (xmask + ymask).astype(bool)
	  	
				
	  print "New Mask,",newSlice,", covers ",nmask.sum(),' of ', len(self.xnc('index')[:])
	
	  if nmask.sum() == len(self.xnc('index')[:]):
		print "New Mask,",newSlice,", covers entire dataset.",nmask.sum(), len(self.xnc('index')[:])
		continue


	  
	  for xkey,ykey in plots:
	  	filename = self.getFileName(newSlice,xkey,ykey)
		print "investigating:", filename
		if not ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):continue
		
		x = self.extractData(self.xnc,mt[self.xtype][self.name], key=xkey)
		y = self.extractData(self.ync,mt[self.ytype][self.name], key=ykey)
		#try:	
		#	x = np.ma.array(self.xnc(xkey)[:])
		#	print "No convension here"
		#except:	x = np.ma.array(self.extractData(self.xnc,mt[self.xtype][self.name]))
		#try:	y = np.ma.array(self.ync(ykey)[:])
		#except:	y = np.ma.array(self.extractData(self.ync,mt[self.ytype][self.name]))

		#data
		#try:	x = self.xnc(mt[self.xtype][self.name][0])[:]
		#except:	x = self.extractData(self.xnc,mt[self.xtype][self.name])	
		#try:	y = self.ync(mt[self.ytype][self.name][0])[:]
		#except:	y = self.extractData(self.ync,mt[self.ytype][self.name])
			
		#xd = self.extractData(self.xnc,mt['Ersem'][self.name])
		#yd = self.extractData(self.ync,mt[self.ytype][self.name])			
		#x = ma(self.xnc(xkey))
		#y = ma(self.ync(ykey))
		
		fullmask = nmask + x.mask + y.mask + np.ma.masked_invalid(x).mask + np.ma.masked_invalid(y).mask 
		
		nmxx = np.ma.masked_where(fullmask, xx).compressed()
		nmxy = np.ma.masked_where(fullmask, xy).compressed()
		nmxz = np.ma.masked_where(fullmask, xz).compressed()
		nmxt = np.ma.masked_where(fullmask, xt).compressed()	
		nmyx = np.ma.masked_where(fullmask, yx).compressed()
		nmyy = np.ma.masked_where(fullmask, yy).compressed()
		nmyz = np.ma.masked_where(fullmask, yz).compressed()
		nmyt = np.ma.masked_where(fullmask, yt).compressed()
		datax = np.ma.masked_where(fullmask, x).compressed()
		datay = np.ma.masked_where(fullmask, y).compressed()
		
		if mt[self.ytype][self.name] == ['Chlorophylla',]:	datay = datay/1000.		
		
		if 0 in [len(datax),len(datay),len(nmxx),len(nmxy),len(nmxz),len(nmyx),len(nmyy),len(nmyz)]:
			print 'slice:',newSlice,'There is a zero in :',	 [len(datax),len(datay),len(nmxx),len(nmxy),len(nmxz),len(nmyx),len(nmyy),len(nmyz)]
			continue	

		if newSlice == 'SalArtifact':
			print 'SalArtifact: x data:',xkey,datax.min(),datax.mean(),datax.max()
			print 'SalArtifact: y data:',ykey,datay.min(),datay.mean(),datay.max()

		try:xunits = fancyUnits(self.xnc.variables[xkey].units,debug=True)
		except:xunits = ''
		try:yunits = fancyUnits(self.ync.variables[ykey].units,debug=True)	
		except:yunits = ''
				
		#try:xunits = fancyUnits(self.xnc.variables[xkey].units)
		#except:xunits=''
		#try:yunits = fancyUnits(self.ync.variables[ykey].units)	
		#except:yunits=''
							
		if xkey == self.kd[self.xtype]['lon']:
		  labelx = "Longitude"
		  labely = "Latitude"
		else:
		  #labelx = self.xtype+' '+getLongName(self.name)+', '+ xunits
		  #labely = self.ytype+' '+getLongName(self.name)+', '+ yunits
		  labelx = self.xtype+' '+self.name+', '+ xunits
		  labely = self.ytype+' '+self.name+', '+ yunits		  
		try: title = getLongName(newSlice)+' '+getLongName(self.name)
		except:title = newSlice+' '+xkey+' vs '+ykey
		

		gs = 50
		robfnx = filename[:-4]+'_xrobin.png'
		robfny = filename[:-4]+'_yrobin.png'
		robfnxy = filename[:-4]+'_xyrobin.png'		

		histfnxy = filename[:-4]+'_hist.png'
		cufn = 	filename.replace('.png','_cusum.png')	
		scat = True
		dmin = min([datax.min(),datay.min()])
		dmax = max([datax.max(),datay.max()])

		if dmin == dmax: 
			print "dmin == dmax"
			continue
		if not datax.size or not datay.size:
			print "No data for this cut"
			continue
					
		documusum = False
		if documusum:			
			cumulativeSumPlot(datax, datay,  cufn, Title=title, labelx='Count',labely='Cumulative Sum',)
			print "Made", cufn

		
		if ukp.shouldIMakeFile([self.xfn,self.yfn],robfnxy,debug=False) or True:
			ti1 = getLongName(self.xtype)+' ' +getLongName(newSlice)+' '+getLongName(self.name)
			ti2 =  getLongName(self.ytype)+' ' +getLongName(newSlice)+' '+getLongName(self.name)	
			if self.name in noXYLogs or dmin*dmax <=0.:
				print "ROBINNOT DOING DOLOG:",[ti1,ti2],False,dmin,dmax
				ukp.robinPlotPair(nmxx, nmxy, datax,datay,robfnxy,titles=[ti1,ti2], vmin=dmin,vmax=dmax, doLog=False)
				
			else:	
				print "ROBIN DOING DOLOG:",[ti1,ti2],False,dmin,dmax
			#	data1 = np.ma.log10(np.ma.masked_less_equal(datax, 0.))
			#	data2 = np.ma.log10(np.ma.masked_less_equal(datay, 0.))
			#	min1 = np.ma.log10(min([data1.min(),data2.min(),]))
			#	max1 = np.ma.log10(max([data1.max(),data2.max(),]))	
			#	ukp.robinPlotPair(nmxx, nmxy, data1, data2,
			#			robfnxy,titles=[ti1,ti2], 
			#			vmin=min1,vmax=max1,
			#			cbarlabel='log$_{10}$('+xunits+')',
			#			doLog=False)
										
				ukp.robinPlotPair(nmxx, nmxy, np.ma.log10(datax),np.ma.log10(datay),
						robfnxy,titles=[ti1,ti2], 
						vmin=np.ma.log10(dmin),vmax=np.ma.log10(dmax),
						cbarlabel='log$_{10}$('+xunits+')',
						doLog=False)
			
		if ukp.shouldIMakeFile([self.xfn,self.yfn],histfnxy,debug=False):
			xaxislabel= getLongName(self.name)+', '+ xunits
			#if self.name in self.BioLogScales: 
			if self.name in noXYLogs:				
				ukp.histPlot(datax, datay,  histfnxy, Title=title, labelx=self.xtype,labely=self.ytype,xaxislabel =xaxislabel)	
			else:	ukp.histPlot(datax, datay,  histfnxy, Title=title, labelx=self.xtype,labely=self.ytype,xaxislabel =xaxislabel, logx = True, )
				
				

										
		
		#if ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):
		if self.name in noXYLogs:
			ukp.scatterPlot(datax, datay,  filename, Title=title, labelx=labelx,labely=labely, bestfitLine=True,gridsize=gs)
		else:	ukp.scatterPlot(datax, datay,  filename, Title=title, labelx=labelx,labely=labely, bestfitLine=True,gridsize=gs,logx = True, logy=True,)					
		#else:
		#			
	
		
		#filename2 = filename[:-4]+'_01-99.png'
		#ukp.scatterPlot(datax, datay,  filename2, Title=title +' 1%-99%', labelx=labelx,labely=labely, bestfitLine=True,gridsize=gs,percentileRange=[1,99])
		
		
		if self.saveShelve:
			if type(newSlice) in [type(['a','b',]),type(('a','b',))]:	
				ns = ''.join(newSlice)
			else: ns = newSlice
			self.shelveName = ukp.folder(['shelves',self.compType, 'Slices'+self.basicCut,self.name])+self.name+'_'+ns+'_'+xkey+'vs'+ykey+'.shelve'
			s = shOpen(self.shelveName)
			print "Saving ",self.shelveName	

			b1, b0, rValue, pValue, stdErr = linregress(datax, datay)
			print "linear regression: \n\tb1:",b1, "\n\tb0:", b0, "\n\trValue:",rValue, "\n\tpValue:",pValue, "\n\tstdErr:",stdErr
			s['b1'] =  b1
			s['b0'] =  b0
			s['rValue'] =  rValue
			s['pValue'] =  pValue
			s['stdErr'] =  stdErr						
			s['N'] 	    =  len(datax)
					
	  		mtaylor = StatsDiagram(datax,datay)
			s['Taylor.E0']=  	mtaylor.E0
			s['Taylor.E']=  	mtaylor.E
			s['Taylor.R']=   	mtaylor.R
			s['Taylor.p']=   	mtaylor.p							
			s['Taylor.gamma']=   	mtaylor.gamma
			
			s['datax'] = datax
			s['datay'] = datay

			s['x_lon'] = nmxx
			s['x_lat'] = nmxy
			s['x_depth'] = nmxz
			s['x_time'] = nmxt			
			s['y_lon'] = nmyx
			s['y_lat'] = nmyy
			s['y_depth'] = nmyz
			s['y_time'] = nmyt					
			
			s['title'] = title
			s['labelx'] = labelx
			s['labely'] = labely
			s['name'] =  self.name
			s['xtype'] =  self.xtype
			s['ytype'] =  self.ytype
			s['xfn'] =  self.xfn
			s['yfn'] =  self.yfn
			s['slice']= newSlice
			s['compType'] = self.compType
			s['basicCut'] = self.basicCut			
			s['newSlice'] = ns
			s['xkey'] = xkey			
			s['ykey'] = ykey
			s.close()
		#waveletA = False
		#if waveletA:
		#	w = wavelet(self.shelveName, plotDecomps= False)




	
def main():
	assert False
	maredatFiles 	= [] 
	ERSEMfiles 	= []

	pfts =	[
		 'iron', 
		'pCO2', 'PP','intPP',

		'phosphateSurface','phosphate100m','phosphate200m','phosphate500m',	
		'silicateSurface','silicate100m','silicate200m','silicate500m',
		'nitrateSurface','nitrate100m','nitrate200m','nitrate500m',	

		 	
		'picophyto',
		'chl',	


	#	]
	#tmp =[
		
		'diatoms',	

		'microzoo',
		'mesozoo',  		
	
		#'iron', 
		#'silicate500m',	'nitrate500m',	
		#'silicate100m','silicate200m',
		#'nitrate100m','nitrate200m',		
		#'silicateSurface','nitrateSurface',

		#'tempSurface',	'salSurface', 			
		#'temp500m',		'temp1000m',	'sal1000m','sal500m',
		#'sal100m','sal200m', 'temp100m','temp200m',			
			
		'bac',

		 
		#'microzoo','mesozoo',  

				  
		#'nitrateAll',	'nitrateTransect',
		#'phosphateSurface',	'phosphateAll',	'phosphateTransect',
		#'silicateSurface',	'silicateAll',	'silicateTransect',
		#		'tempAll',	'tempTransect',	
		#'salSurface', 		'salAll',	'salTransect',
		 #'mld_DT02','mld_DR003','mld_DReqDTm02'		 				
		#'Seawifs',
#		# 'Seawifs-biomass',
#		 'Seawifs-biomass',		 
#			
#		 'Seawifs-nano','Seawifs-pico','Seawifs-micro',
#		 'SeawifsBM-nano','SeawifsBM-pico','SeawifsBM-micro',
		 ]
	plotallcuts = 0#True
	
	
	try: 	
		jobID = argv[1]
		key   = argv[2]
		print "Using command line arguments:", argv[1], argv[2]
	except:
		#jobID= 'xhonp'
		#key = "clim"
		print "Not using command line arguments,Defaults sweep"#,jobID
		for jobID in [  'xhont', 'xhonu', 'xhonv', 'xhonp',]:
		  for yr in ['1894','1899']:
		    for name in pfts:
			m = makePlots(name, jobID+'_'+yr,plotallcuts=plotallcuts, compareCoords=True)		
		
	jobkey = jobID+'_'+key
	
	for name in pfts:
		m = makePlots(name, jobkey,plotallcuts=plotallcuts, compareCoords=True)
		#return
	#	m = makePlots(ERSEMfile,maredatFile,plotallcuts=plotallcuts)

		
if __name__=="__main__":
	main() 
	print 'The end.'
	
