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
  def __init__(self,matchedDataFile,matchedModelFile, name, jobID='MEDUSA',year='clim',region='', compareCoords=True,saveShelve=True,plotallcuts=  False): #xfilename,yfilename,

  
  	self.xfn =matchedModelFile
  	self.yfn =matchedDataFile  	
    	self.name = name
  	self.xtype = jobID
	if self.name == 'chl': self.ytype = 'Maredat'

	self.compType = self.xtype+self.ytype+region

  	self.saveShelve = saveShelve
  	self.kd = getkd()
  	self.BioLogScales = BioLogScales
  	self.xnc = ncdfView(self.xfn,Quiet=True)
  	self.ync = ncdfView(self.yfn,Quiet=True)


	
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

	if compareCoords: self.CompareCoords()	
	self.plotWithSlices()
	    	
  	self.xnc.close()
  	self.ync.close()  	
 
	

  def CompareCoords(self,):
	#xkeys,ykeys = ['index_t',] , ['index_t',]
	#for k in ['t','lat','lon','z','lon',]:	xkeys.append(self.kd[self.xtype][k])
	#for k in ['t','lat','lon','z','lat',]:	ykeys.append(self.kd[self.ytype][k])
	
	xkeys = [self.kd[self.xtype][k] for k in ['index_t','t','lat','lon','z','lon',]]
	ykeys = [self.kd[self.ytype][k] for k in ['index_t','t','lat','lon','z','lat',]]	
	  	 	  	
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
		

		mask = np.ma.array(self.xnc(xkey)[:]).mask + np.ma.array(self.ync(ykey)[:]).mask
		dx = np.ma.masked_where(mask, np.ma.array(self.xnc(xkey)[:])).compressed()
		dy = np.ma.masked_where(mask, np.ma.array(self.ync(ykey)[:])).compressed()
				
		print "CompareCoords:\t",xkey,':', len(dx),"\t:",ykey,':',len(dy), dx.min(),dx.max(), dy.min(), dy.max()		
	


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
		filename = ukp.folder([file_prefix,'months'])+self.name+'_'+ukp.mnStr(self.months[newSlice])+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
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
	except:	xd = extractData(self.xnc,mt[self.xtype][self.name])	
	try:	yd = self.ync(mt[self.ytype][self.name][0])[:]
	except:	yd = extractData(self.ync,mt[self.ytype][self.name])	
	print "plotWithSlices:\tx",xd.min(),xd.mean(),xd.max()
	print "plotWithSlices:\ty",yd.min(),yd.mean(),yd.max()
	
	if mt[self.ytype][self.name] == ['Chlorophylla',]:	yd = yd/1000.

	basicCutMask = ukp.makeMask(self.name,'All',xt,xz,xy,xx,xd).astype(int)
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
		basicCutMask += ukp.makeMask(self.name,ns,xt,xz,xy,xx,xd).astype(int)
	   	basicCutMask += ukp.makeMask(self.name,ns,yt,yz,yy,yx,yd).astype(int)
	   
	
	for newSlice in self.newSlices:
	  xmask = ukp.makeMask(self.name,'All',xt,xz,xy,xx,xd).astype(int)
	  ymask = ukp.makeMask(self.name,'All',yt,yz,yy,yx,yd).astype(int)
	  	  	
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
	  	  xmask += ukp.makeMask(self.name,n,xt,xz,xy,xx,xd).astype(int)	  
	  	  ymask += ukp.makeMask(self.name,n,yt,yz,yy,yx,yd).astype(int)
		  
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
	 	    xmask = ukp.makeMask(self.name,'All',xt,xz,xy,xx,xd).astype(int)
		    ymask = ukp.makeMask(self.name,'All',yt,yz,yy,yx,yd).astype(int)	  	
	  	    for ns in self.standardCuts: #['5-95pc','ignoreInlandSeas','OffShelf',]:
			if self.name in ['tempSurface','tempTransect', 'tempAll'] and ns in ['aboveZero',]:continue # no negative or zero values allowed.	    	  	    
			if self.name not in ['nitrateSurface','nitrateAll','nitrateTransect',] 		and ns in ['0.1','0.2']:continue	    			
			if self.name not in ['phosphateSurface','phosphateAll','phosphateTransect',] 	and ns in ['0.01',]:continue	    						
	  		xmask += ukp.makeMask(self.name,ns,xt,xz,xy,xx,xd).astype(int)
	  	 	ymask += ukp.makeMask(self.name,ns,yt,yz,yy,yx,yd).astype(int)
	  	 	
	  elif newSlice in ['All','aboveZero',]:  	
	 	xmask = ukp.makeMask(self.name,newSlice,xt,xz,xy,xx,xd).astype(int)
		ymask = ukp.makeMask(self.name,newSlice,yt,yz,yy,yx,yd).astype(int)	  	
	  else:
	  	xmask += ukp.makeMask(self.name,newSlice,xt,xz,xy,xx,xd).astype(int)
	  	ymask += ukp.makeMask(self.name,newSlice,yt,yz,yy,yx,yd).astype(int)
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
		
		x = extractData(self.xnc,mt[self.xtype][self.name], key=xkey)
		y = extractData(self.ync,mt[self.ytype][self.name], key=ykey)
		
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

def extractData(nc, mt,key = ['',]):
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
	
