#!/usr/bin/ipython
from netCDF4 import num2date
from datetime import datetime
from sys import argv
from os.path import exists,split, getmtime, basename
from glob import glob
from shelve import open as shOpen
from matplotlib.colors import LogNorm
from matplotlib import pyplot, ticker
from calendar import month_name
from itertoolsmodule import product
from scipy.stats import linregress
from scipy.stats.mstats import scoreatpercentile
import numpy as np 

# imports from other gitlab repositories:
from ncdfView import ncdfView
from StatsDiagram import StatsDiagram

#local imports
import UKESMpython as ukp 
from pftnames import getLongName, AutoVivification, getmt,fancyUnits,MaredatTypes,MLDTypes,WOATypes,GEOTRACESTypes


    

"""	This code makes matched plots, hexbins, scatter plots, and so on.

"""

#BioLogScales 	= ['bac','mesozoo','diatoms','picophyto','microzoo','PP','Seawifs', 'iron'] 
	
noXYLogs 	= [ 'pCO2',
		#'nitrateSurface', 	'nitrateAll',	'nitrateTransect',
		#'phosphateSurface',	'phosphateAll',	'phosphateTransect',
		'silicateSurface',	'silicateAll',	'silicateTransect', 'silicate100m','silicate200m','silicate500m',
		'tempSurface',		'tempAll',	'tempTransect',	'temp100m',	'temp200m','temp1000m',	'temp500m',
		'salSurface', 		'salAll',	'salTransect',	'sal100m',	'sal200m','sal1000m',	'sal500m',]
		


class makePlots:
  def __init__(self,matchedDataFile,matchedModelFile, name, jobID='MEDUSA',year='clim',region='', compareCoords=True,workingDir='',imageDir='',plotallcuts=  True): #xfilename,yfilename,saveShelve=True,

  
  	self.xfn =matchedModelFile
  	self.yfn =matchedDataFile  	
    	self.name = name
  	self.xtype = jobID
  	self.year = year
  	self.plotallcuts = plotallcuts
  	self.compareCoords = compareCoords
  	
	if self.name in MaredatTypes:  	self.ytype = 'Maredat'
	if self.name in WOATypes:  	self.ytype = 'WOA'
	if self.name in MLDTypes:  	self.ytype = 'IFREMER'	
	if self.name in GEOTRACESTypes: self.ytype = 'GEOTRACES'		
	
  	self.shelvedir = workingDir
  	if self.shelvedir == '':self.shelvedir = ukp.folder(['shelves',self.xtype,self.year,self.ytype, 'Slices',self.name])
  	else:			self.shelvedir = ukp.folder(self.shelvedir)		

	if imageDir=='':	self.imageDir = ukp.folder(['images',self.xtype,'P2P_plots',self.year,self.name])
	else: 			self.imageDir = ukp.folder(imageDir)

	self.run()
	
	
	
	
  def run(self,):

  	self.xnc = ncdfView(self.xfn,Quiet=True)
  	self.ync = ncdfView(self.yfn,Quiet=True)
	self.mt = getmt()		    	
	if self.compareCoords: self.CompareCoords()	
	self.defineSlices(self.plotallcuts)
	
	
	#for nslice in self.newSlices:
	self.plotWithSlices()

	    	
  	self.xnc.close()
  	self.ync.close()  	


  		


  def plotWithSlices(self):#,newSlice,):  
	print "plotWithSlices:\txtype:",self.xtype,"\tytype:",self.ytype,"\tname:",self.name
  	
	#####
	# Test if any of the plots exist.
	  	
	xkeys = []
	ykeys = []
	plotpairs = [] 

	
	nx = self.mt[self.xtype][self.name]
	if type(nx) == type(['a',]):	xkeys = self.mt[self.xtype][self.name]
	else:				xkeys.append(self.mt[self.xtype][self.name]['name'])
	ny = self.mt[self.ytype][self.name]
	if type(ny) == type(['a',]):	ykeys = self.mt[self.ytype][self.name]
	else:				ykeys.append(self.mt[self.ytype][self.name]['name'])	

	print "plotWithSlices:\txkeys:", xkeys,'\tykeys:', ykeys
	
	
	
	#####
	# This section of code is a bit of a time saver.
	# It checks to see if the image and the output shelve exist.
	# If they both exist and and are older than the input netcdfs, the rest of this function is skipped.
	# If one is missing, or the input files are newer than the old image, the function runs as normal.
	# Caveat: if some image can not be made, ie the data makes the mask cover 100% of the data, then the code will run as normal (no skipping).  
	self.shelvesAV = AutoVivification()
	
	plotsToMake=0
	for newSlice in self.newSlices:	
	    for xk,yk in product(xkeys,ykeys):
	  	print 'plotWithSlices:\tlisting plotpairs:\tX', xk,': self.mt[',self.xtype,'][',self.name,']'
	  	print 'plotWithSlices:\tlisting plotpairs:\tY', yk,': self.mt[',self.ytype,'][',self.name,']'	 
		plotpairs.append((xk,yk))
		print xk,yk,self.xtype,self.ytype,self.name
		try:fn = newSlice+'_'+xk+'vs'+yk
	  	except:
	  		print "ERROR:\tcan\'t add ",newSlice,xk,yk, 'together as strings. the problem is probably in your mt dictionary in pftnames.'
			assert False
			
		#####
		# Does the image exist?	
		filename = self.getFileName(newSlice,xk,yk)
		if ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):
			plotsToMake+=1
		
		#####
		#Does the shelve file exist?
		if type(newSlice) in [type(['a','b',]),type(('a','b',))]:	
			ns = ''.join(newSlice)
		else: 	ns = newSlice			
		shelveName = self.shelvedir +self.name+'_'+ns+'_'+xk+'vs'+yk+'.shelve'
		if ukp.shouldIMakeFile([self.xfn,self.yfn],shelveName,debug=False):
			plotsToMake+=1
			
		#####
		# Make a list of shelve files for the target plots.
		self.shelvesAV[newSlice][xk][yk] = shelveName				
		try:	self.shelves.append(shelveName)
		except:	self.shelves = [shelveName,]

		
	if plotsToMake == 0: 
	  	print 'plotWithSlices:\tAll plots and shelve files already made',self.name, newSlice, xkeys,ykeys
		return
	
	#####
	# Load data
	
	#time and depth
	xt = self.xnc(self.mt[self.xtype]['t'])[:]
	yt = self.ync(self.mt[self.ytype]['t'])[:]
	xz = self.xnc(self.mt[self.xtype]['z'])[:]
	yz = self.ync(self.mt[self.ytype]['z'])[:]

	#lat and lon
	xy = self.xnc(self.mt[self.xtype]['lat'])[:]
	yy = self.ync(self.mt[self.ytype]['lat'])[:]
	xx = self.xnc(self.mt[self.xtype]['lon'])[:]
	yx = self.ync(self.mt[self.ytype]['lon'])[:]	

	self.xx = xx
	self.xy = xy
	self.xz = xz
	self.xt = xt	
	self.yx = yx
	self.yy = yy
	self.yz = yz
	self.yt = yt	
	
	
	for newSlice in self.newSlices:	
	    for xkey,ykey in product(xkeys,ykeys):
	    	print "plotWithSlices:\t", newSlice, xkey,ykey
		self.plotsFromKeys(newSlice,xkey,ykey)



  def plotsFromKeys(self,newSlice,xkey,ykey):	     

	#####
	# check that the plot and shelve should be made 
	if type(newSlice) in [type(['a','b',]),type(('a','b',))]:	
		ns = ''.join(newSlice)
	else: ns = newSlice
	self.shelveName = self.shelvedir +self.name+'_'+ns+'_'+xkey+'vs'+ykey+'.shelve'		
 	filename = self.getFileName(newSlice,xkey,ykey)
 	
	print "plotWithSlices:\tINFO:\tinvestigating:", filename
	if not ukp.shouldIMakeFile([self.xfn,self.yfn],self.shelveName,debug=False) \
		and not ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):return
	
	
	#####
	# Extract remaining data (already know lat,lon,time,depth)
	xd = extractData(self.xnc,self.mt[self.xtype][self.name],key = xkey)	
	yd = extractData(self.ync,self.mt[self.ytype][self.name],key = ykey)
 
	
	#####
	# Build mask
	fullmask = xd.mask + yd.mask + np.ma.masked_invalid(xd).mask + np.ma.masked_invalid(yd).mask
	
	if type(newSlice) in [type(['a',]),type(('a',))]:    	# newSlice is actaully a list of multiple slices.
	   	for n in newSlice:
	  		fullmask += ukp.makeMask(self.name,n,self.xt,self.xz,self.xy,self.xx,xd).astype(int)	  
		  	fullmask += ukp.makeMask(self.name,n,self.yt,self.yz,self.yy,self.yx,yd).astype(int)	  
		  	
	elif newSlice == 'Standard':				# Standard is a shorthand for my favourite cuts.
	  	for stanSlice in self.standardCuts: 
			if self.name in ['tempSurface','tempTransect', 'tempAll'] and stanSlice in ['aboveZero',]:continue 
				    						
	  		fullmask += ukp.makeMask(self.name,stanSlice,self.xt,self.xz,self.xy,self.xx,xd).astype(int)
	  	 	fullmask += ukp.makeMask(self.name,stanSlice,self.yt,self.yz,self.yy,self.yx,yd).astype(int)	
	  	 	
	else:  	# newSlice is a simple slice.
	  	fullmask += ukp.makeMask(self.name,newSlice,self.xt,self.xz,self.xy,self.xx,xd).astype(int)
	  	fullmask += ukp.makeMask(self.name,newSlice,self.yt,self.yz,self.yy,self.yx,yd).astype(int)
	  	print 'plotWithSlices:\t',fullmask.sum()

	  
        if self.name in ['mld','mld_DT02','mld_DR003','mld_DReqDTm02']:
        	mldMask = self.ync('mask')[:]
        	fullmask += np.ma.masked_where(mldMask==0.,mldMask).mask        
	  
	
	N = len(self.xt)			


	if fullmask.sum() >= N:
		print "plotWithSlices:\tNew Mask,",newSlice,", covers entire dataset.",fullmask.sum(), N
		try:	self.shelves[newSlice][xk][yk] = ''
		except:	pass			
		return
	print "plotWithSlices:\tNew Mask,",newSlice,", covers ",fullmask.sum(),' of ', N
		
	#####
	# Apply mask to all data.	
	nmxx = np.ma.masked_where(fullmask, self.xx).compressed()
	nmxy = np.ma.masked_where(fullmask, self.xy).compressed()
	nmxz = np.ma.masked_where(fullmask, self.xz).compressed()
	nmxt = np.ma.masked_where(fullmask, self.xt).compressed()	
	nmyx = np.ma.masked_where(fullmask, self.yx).compressed()
	nmyy = np.ma.masked_where(fullmask, self.yy).compressed()
	nmyz = np.ma.masked_where(fullmask, self.yz).compressed()
	nmyt = np.ma.masked_where(fullmask, self.yt).compressed()
	datax = np.ma.masked_where(fullmask, xd).compressed()
	datay = np.ma.masked_where(fullmask, yd).compressed()
		
	if 0 in [len(datax),len(datay),len(nmxx),len(nmxy),len(nmxz),len(nmyx),len(nmyy),len(nmyz)]:
		print 'plotWithSlices:\tWARNING:\tslice:',newSlice,'There is a zero in one of the fields.' 
		try:	self.shelvesAV[newSlice][xk][yk] = ''			
		except:	pass			
		return	
						
	dmin = min([datax.min(),datay.min()])
	dmax = max([datax.max(),datay.max()])
	if dmin == dmax: 
		print "plotWithSlices:\tWARNING:\tminimum == maximum,\t (",dmin,' == ',dmax,')'
		try:	self.shelvesAV[newSlice][xk][yk] = ''
		except:	pass			
		return
			

	#####
	# Prepare units, axis labels and titles.
	try:    xunits = fancyUnits(self.mt[self.xtype][self.name]['units'])
	except: xunits = fancyUnits(self.xnc.variables[xkey].units,debug=True)

	try:   yunits = fancyUnits(self.mt[self.ytype][self.name]['units'])
	except:yunits = fancyUnits(self.ync.variables[ykey].units,debug=True)	

	labelx = getLongName(self.xtype)+' '+getLongName(self.name)+', '+ xunits
	labely = getLongName(self.ytype)+' '+getLongName(self.name)+', '+ yunits	
		
	try: title = getLongName(newSlice)+' '+getLongName(self.name)
	except:title = newSlice+' '+xkey+' vs '+ykey


	robfnxy  = filename.replace('.png','_xyrobin.png')
	histfnxy = filename.replace('.png','_hist.png')
				
	#####
	# Robinson projection plot		
	if ukp.shouldIMakeFile([self.xfn,self.yfn],robfnxy,debug=False) or True:
		ti1 = getLongName(self.xtype)+' ' +getLongName(newSlice)+' '+getLongName(self.name)
		ti2 =  getLongName(self.ytype)+' ' +getLongName(newSlice)+' '+getLongName(self.name)	
		if self.name in noXYLogs or dmin*dmax <=0.:
			print "plotWithSlices:\tROBIN NOT DOING DOLOG:",[ti1,ti2],False,dmin,dmax
			ukp.robinPlotPair(nmxx, nmxy, datax,datay,robfnxy,titles=[ti1,ti2], vmin=dmin,vmax=dmax, doLog=False)
				
		else:	
			print "plotWithSlices:\tROBIN DOING DOLOG:",[ti1,ti2],False,dmin,dmax					
			ukp.robinPlotPair(nmxx, nmxy, np.ma.log10(datax),np.ma.log10(datay),
						robfnxy,titles=[ti1,ti2], 
						vmin=np.ma.log10(dmin),vmax=np.ma.log10(dmax),
						cbarlabel='log$_{10}$('+xunits+')',
						doLog=False)
		
	#####
	# Simultaneous histograms plot	
	if ukp.shouldIMakeFile([self.xfn,self.yfn],histfnxy,debug=False):
		xaxislabel= getLongName(self.name)+', '+ xunits
		if self.name in noXYLogs or dmin*dmax <=0.:				
			ukp.histPlot(datax, datay,  histfnxy, Title=title, labelx=self.xtype,labely=self.ytype,xaxislabel =xaxislabel)	
		else:	ukp.histPlot(datax, datay,  histfnxy, Title=title, labelx=self.xtype,labely=self.ytype,xaxislabel =xaxislabel, logx = True, )
				
	#####
	# Scatter  (hexbin) plot
	if ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):		
		gs = 50					
		if self.name in noXYLogs or dmin*dmax <=0.:
			ukp.scatterPlot(datax, datay,  filename, Title=title, labelx=labelx,labely=labely, bestfitLine=True,gridsize=gs)
		else:	ukp.scatterPlot(datax, datay,  filename, Title=title, labelx=labelx,labely=labely, bestfitLine=True,gridsize=gs,logx = True, logy=True,)

	#####
	# Save fit in a shelve file.		
	s = shOpen(self.shelveName)
	print "plotWithSlices:\tSaving ",self.shelveName	
	b1, b0, rValue, pValue, stdErr = linregress(datax, datay)
	print "plotWithSlices:\tlinear regression: \n\tb1:",b1, "\n\tb0:", b0, "\n\trValue:",rValue, "\n\tpValue:",pValue, "\n\tstdErr:",stdErr
	s['b1'] 	=  b1
	s['b0'] 	=  b0
	s['rValue'] 	=  rValue
	s['pValue'] 	=  pValue
	s['stdErr'] 	=  stdErr						
	s['N'] 	    	=  len(datax)
					
  	mtaylor = StatsDiagram(datax,datay)
	s['Taylor.E0'] 	= mtaylor.E0
	s['Taylor.E']	= mtaylor.E
	s['Taylor.R']	= mtaylor.R
	s['Taylor.p']	= mtaylor.p							
	s['Taylor.gamma']=mtaylor.gamma
			
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
			
	s['title'] = 	title
	s['labelx'] = 	labelx
	s['labely'] = 	labely
	s['name'] =   	self.name
	s['year'] =   	self.year 
	s['xtype'] =  	self.xtype
	s['ytype'] =  	self.ytype
	s['xfn'] =  	self.xfn
	s['yfn'] =  	self.yfn
	s['slice']= 	newSlice
	s['newSlice'] = ns
	s['xkey'] = 	xkey			
	s['ykey'] = 	ykey
	s.close()

 
  	

  def CompareCoords(self,):
	"""	This routine plots the coordinates of the data against the coordinates of the model.
		This should produce a straight line plot, ensuring that the matching has been performed correctly.
	"""
	
	xcoords = [self.mt[self.xtype][k] for k in ['index_t','t','lat','lon','z','lon',]]
	ycoords = [self.mt[self.ytype][k] for k in ['index_t','t','lat','lon','z','lat',]]	
	  	 	  	
  	for xkey,ykey in zip(xcoords,ycoords):
	    	if xkey not in self.xnc.variables.keys():continue  	    
	    	if ykey not in self.ync.variables.keys():continue
		filename = self.imageDir+'CompareCoords'+self.name+xkey+'vs'+ykey+'.png'	    	
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


  def defineSlices(self,plotallcuts):	
	self.newSlices 		=['All','Standard',]

	self.standardCuts = ['5-95pc','ignoreInlandSeas','OffShelf','ignoreExtraArtics','aboveZero',]	

	self.months = {month_name[i+1]:i for i in xrange(0,12) }
	
	self.depthRanges	=['OffShelf','maskBelowBathy', 'OnShelf',] 
				 # 'Depth_0-10m','Depth_10-20m','Depth_20-50m','Depth_50-100m','Depth_100-500m','Depth_500m',
	self.percentiles	=['0-1pc','1-5pc','5-25pc',
				  '25-40pc','40-60pc','60-75pc',
				  '75-95pc','95-99pc','99-100pc',]
	self.latregions		=['NorthTemperate','SouthTemperate','NorthTropics',
				  'Equatorial',  'SouthTropics','Antarctic',
				  'NorthArctic',]
				  #'Arctic','Tropics','Temperate']
	self.Seas		=['ignoreMediteranean','BlackSea','ignoreBlackSea',
				  'RedSea','BalticSea','PersianGulf',
				  'ignoreInlandSeas',]	
				  # 'ignoreRedSea', 'ignoreBalticSea','ignorePersianGulf',]
	self.Oceans		=['SouthPacificOcean',  'ArcticOcean',
				  'AntarcticOcean','NorthAtlanticOcean','SouthAtlanticOcean',
				 'NorthPacificOcean','IndianOcean',] 
				 #'ignoreExtraArtics','ignoreMidArtics','ignoreArtics','ignoreMoreArtics',

	self.QualityCuts 	=['Overestimate','Underestimate','Overestimate_2sig',
				  'Underestimate_2sig','Overestimate_3sig','Underestimate_3sig', 
				  'Matched','OffAxis','1-99pc',
				  '5-95pc','0-99pc',]
	self.OceanMonths	= sorted([i for i in product(self.Oceans,self.months)] )
	self.OceanMonths.extend(sorted([i for i in product(['All',],self.months)]))
	
	if plotallcuts:
		 self.plotMonths	= True
		 self.plotdepthRanges	=0
		 self.plotpercentiles	=0#True	
		 self.plotLatRegions	=0
		 self.plotQualityCuts	=0#True	
		 self.plotSeas		=0#True		 
		 self.plotOceans	= True	
		 self.plotOceanMonths   = 0#True	
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
	#print "defineSlices:\tSLICES:", 	 self.newSlices
	
 		
  def getFileName(self,newSlice,xkey,ykey):
	file_prefix = ukp.folder(['images',self.xtype,'P2P_plots',self.year,self.name,])

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
	  	print 'getFileName:', newSlice,	''.join(newSlice)  ,xkey,ykey
	  	try:fn = newSlice+'_'+xkey+'vs'+ykey
	  	except:
	  		print "ERROR:\tcan't add ",newSlice,xkey,ykey, 'together as strings. It breaks in getFileName, but the problem is probably in your mt dictionary in pftnames'
		filename = file_prefix+self.name+'_'+newSlice+'_'+xkey+'vs'+ykey+file_suffix
	return filename



		

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

	if mt['mul1000']:
  		print "Extracting data:\tMultipling by 1000. ", mt['mul1000']
   	  	xd = nc(mt['mul1000'][0])[:]*1000.   	



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



		
if __name__=="__main__":
	print "makePlots isn't written to be run as a __main__"
	print "Look at testsuite_p2p.py for examples on how to run this."
	print 'The end.'
	
