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

import numpy as np
from shelve import open as shOpen
from netCDF4 import Dataset,num2date
import os

#Specific local code:
import UKESMpython as ukp

import timeseriesTools as tst 
import timeseriesPlots as tsp 
#getTimes, loadData




class timeseriesAnalysis:
  def __init__(self,
  		modelFiles, 
		dataFile,
		dataType	= '',
		modelcoords 	= '',
		modeldetails 	= '',
		datacoords 	= '',
		datadetails 	= '',								
		datasource	= '',
		model 		= '',
		jobID		= '',
		layers	 	= '',
		regions	 	= '',			
		metrics	 	= '',
		workingDir	= '',
		imageDir	= '',						
		grid		= '',
		gridFile	= '',
		clean		= True,
		debug		= True,
		):
		
	#####
	#	This is the class that does most of the legwork.
	#	First we save all the initialisation settings as class attributes.
		
	
	if debug: print "timeseriesAnalysis:\t init."	
	self.modelFiles 	= modelFiles 		
	self.dataFile		= dataFile
	self.dataType		= dataType
	self.modelcoords 	= modelcoords		
	self.modeldetails 	= modeldetails
	self.datacoords 	= datacoords
	self.datadetails 	= datadetails						
	self.datasource		= datasource
	self.model 		= model
	self.jobID		= jobID
	self.layers	 	= layers
	self.regions	 	= regions			
	self.metrics	 	= metrics						
	self.grid		= grid
	self.gridFile		= gridFile
	self.workingDir		= workingDir
  	self.imageDir 		= imageDir
	self.debug		= debug
	self.clean		= clean
		
  	self.shelvefn 		= ukp.folder(self.workingDir)+'_'.join([self.jobID,self.dataType,])+'.shelve'
	self.shelvefn_insitu	= ukp.folder(self.workingDir)+'_'.join([self.jobID,self.dataType,])+'_insitu.shelve'

	#####
	# Load Data file	
 	self.loadData()
	
	#####
	# Load Model File
  	self.loadModel()  	

	#####
	# Make the plots:
  	self.makePlots()
  	
  	
  	
  	
  def loadModel(self):
	if self.debug: print "timeseriesAnalysis:\tloadModel."		
	####
	# load and calculate the model info
	try:
		if self.clean: 
			print "User requested clean run. Wiping old data."
			assert 0		
		sh = shOpen(self.shelvefn)
		readFiles 	= sh['readFiles']
		modeldataD 	= sh['modeldata']
		sh.close()
		print "Opened shelve:", self.shelvefn, '\tread', len(readFiles)
	except:
		readFiles = []
		modeldataD = {}
		for r in self.regions:
		 for l in self.layers:
		  for m in self.metrics:
		   	modeldataD[(r,l,m)] = {}
		   	
		print "Could not open shelve:", self.shelvefn, '\tread', len(readFiles)	

	###############
	# Check whethere there has been a change in what was requested:
	for r in self.regions:
	  for l in self.layers:
	    for m in self.metrics:
	    	if self.debug:print "Checking: ",[r,l,m,],'\t',
	    	try:
	    		if self.debug: print 'has ', len(modeldataD[(r,l,m)].keys()), 'keys'
	    	except: 
	    		readFiles = []
	    		modeldataD[(r,l,m)] = {}
	    		if self.debug: print 'has no keys'
	    	try:	
	    	    	if len(modeldataD[(r,l,m)].keys()) == 0: 
	    	    		readFiles = []
	    	except: pass

	#####
	# Summarise checks
	if self.debug:	
		print "loadModel:post checks:"
		#print "modeldataD:",modeldataD
		print "shelveFn:",self.shelvefn
		print "readFiles:",readFiles

	
	###############
	# Load files, and calculate fields.
	openedFiles = 0					
	for fn in self.modelFiles:
		if fn in readFiles:continue
		print "loadModel:\tloading new file:",fn,
		nc = Dataset(fn,'r')
		ts = tst.getTimes(nc,self.modelcoords)
		meantime = np.mean(ts)
		print "\ttime:",meantime
		
		DL = tst.DataLoader(fn,nc,self.modelcoords,self.modeldetails, regions = self.regions, layers = self.layers,)
		
	
		for l in self.layers:		
	 	    counts =0		
		    for r in self.regions:
		    	
		    	#####
		    	# Check wherether you can skip loading this metric,region,layer
			skip = True
			for m in self.metrics:
				if skip == False:continue
				try: 
					a = modeldataD[(r,l,m)][meantime]
					print "Already created ",int(meantime),':\t',(r,l,m),'\t=',a
				except: 
					skip = False
					print "Need to create ",int(meantime),':\t',(r,l,m)
			if skip: continue
			
		    	#####
		    	# can't skip it, need to load it.
			layerdata = DL.load[(r,l)]
			if len(layerdata)==0:
				layerdata  = np.ma.array([-999,],mask=[True,])
		  	for m in self.metrics:
		  		try:
		  			a = modeldataD[(r,l,m)][meantime]
		  			continue
		  		except:pass
				if m == 'mean':   	modeldataD[(r,l,m)][meantime] = np.ma.mean(layerdata)
				if m == 'median':   	modeldataD[(r,l,m)][meantime] = np.ma.median(layerdata)
				if m == 'sum':   	modeldataD[(r,l,m)][meantime] = np.ma.sum(layerdata)
		  		print "Loaded metric:", int(meantime),'\t',[(r,l,m)], '\t',modeldataD[(r,l,m)][meantime]
				counts+=1
						
		readFiles.append(fn)		
		openedFiles+=1			


		nc.close()
		if openedFiles:
			print "Saving shelve:", self.shelvefn, '\tread', len(readFiles)				
			sh = shOpen(self.shelvefn)
			sh['readFiles']		= readFiles
			sh['modeldata'] 	= modeldataD
			sh.close()
			openedFiles=0	
		#self.modeldataD = modeldataD			
  		#if len(readFiles)>1: self.makePlots()		
	if openedFiles:
		print "Saving shelve:", self.shelvefn, '\tread', len(readFiles)				
		sh = shOpen(self.shelvefn)
		sh['readFiles']		= readFiles
		sh['modeldata'] 	= modeldataD
		sh.close()
	
	self.modeldataD = modeldataD
	if self.debug: print "timeseriesAnalysis:\tloadModel.\t Model loaded:",	self.modeldataD.keys()[:3], '...', len(self.modeldataD.keys())	



  def loadData(self):
  	
	if self.debug: print "timeseriesAnalysis:\t loadData.",self.dataFile		
	
  	if not self.dataFile: 
 		if self.debug: print "timeseriesAnalysis:\t No data File provided:",self.dataFile		 		
		self.dataD = {}
		return

  	if not os.path.exists(self.dataFile): 
 		if self.debug: print "timeseriesAnalysis:\tWARNING:\t No such data File:",self.dataFile		 
		self.dataD = {}
		return
				
	###############
	# load and calculate the real data info
	try:
		if self.clean: 
			print "timeseriesAnalysis:\t loadData\tUser requested clean run. Wiping old data."
			assert 0		
		sh = shOpen(self.shelvefn_insitu)
		dataD 	= sh['dataD']
		sh.close()
		print "timeseriesAnalysis:\t loadData\tOpened shelve:", self.shelvefn_insitu
		self.dataD = dataD
	except:
		dataD = {}
		print "timeseriesAnalysis:\t loadData\tCould not open shelve:", self.shelvefn_insitu

	###############
	# Test to find out if we need to load the netcdf, or if we can just return the dict as a self.object.
	needtoLoad = False
	for r in self.regions:
	    for l in self.layers:
	    	try:	
	    		print r,l, len(self.dataD[(r,l)]),self.dataD[(r,l)].shape
	    	except: 
			needtoLoad=True
			
	if needtoLoad: pass	
	else:
		self.dataD = dataD	
		return
	###############
	# Loading data for each region.
	print "timeseriesAnalysis:\t loadData,\tloading ",self.dataFile
	#nc = Dataset(self.dataFile,'r')
	#data = tst.loadData(nc, self.datadetails)
				
	
	###############
	# Loading data for each region.
	dl = tst.DataLoader(self.dataFile,'',self.datacoords,self.datadetails, regions = self.regions, layers = self.layers,)
	
	for r in self.regions:
	    for l in self.layers:
	    	dataD[(r,l)] = dl.load[(r,l,)]	
	    	dataD[(r,l,'lat')] = dl.load[(r,l,'lat')]		    	
	    	dataD[(r,l,'lon')] = dl.load[(r,l,'lon')]
		if len(dataD[(r,l)])==0:
			dataD[(r,l)]  = np.ma.array([-999,],mask=[True,])	
			dataD[(r,l,'lat')]  = np.ma.array([-999,],mask=[True,])	    	
			dataD[(r,l,'lon')]  = np.ma.array([-999,],mask=[True,])	    	
									    	
    		print "timeseriesAnalysis:\t loadData,\tloading ",(r,l),  dataD[(r,l)].min(),  dataD[(r,l)].max(),  dataD[(r,l)].mean()
    		
	###############
	# Savng shelve		
	print "timeseriesAnalysis:\t loadData.\tSaving shelve:", self.shelvefn_insitu			
	sh = shOpen(self.shelvefn_insitu)
	sh['dataD'] 	= dataD
	sh.close()
	 	
	self.dataD = dataD






  def mapplots(self, layer,filename):
  	"""	Makes a surface plot of model vs data. 
  	"""	
  	mnc = Dataset(self.modelFiles[-1],'r')
  	
  	modeldata =  tst.getHorizontalSlice(mnc,self.modelcoords,self.modeldetails,layer,data = '').squeeze()
  	modellat = ukp.extractData(mnc,[],key=self.modelcoords['lat'])
  	modellon = ukp.extractData(mnc,[],key=self.modelcoords['lon']) 
  	
	if self.dataFile:
	  	dnc = Dataset(self.dataFile,'r')
  		datadata =  tst.getHorizontalSlice(dnc,self.datacoords,self.datadetails,layer,data = '').squeeze()
  		datalat = ukp.extractData(dnc,[],key=self.datacoords['lat'])
  		datalon = ukp.extractData(dnc,[],key=self.datacoords['lon']) 
		
	else:
		datadata = np.ma.array([-1000,],mask=[True,])
		datalat  = np.ma.array([-1000,],mask=[True,])
		datalon  = np.ma.array([-1000,],mask=[True,])

	if modeldata.ndim ==4:   modeldata=modeldata.mean(0)
	if modeldata.ndim ==3:   modeldata=modeldata.mean(0)  	
	if datadata.ndim ==4:   datadata=datadata.mean(0)
	if datadata.ndim ==3:   datadata=datadata.mean(0)  
	
	
	titles = [' '.join([self.model,'('+self.jobID+')',str(layer),self.modeldetails['name']]),
		  ' '.join([self.datasource,str(layer),self.datadetails['name']])]
  	tsp.mapPlotPair(modellon, modellat, modeldata,
  			datalon,datalat,datadata,
  			filename,
  			titles	= titles,
  			lon0=0.,drawCbar=True,cbarlabel='',dpi=100,)
  	

  def mapplotsRegionsLayers(self, region,layer,filename):
  	"""	Makes a surface plot of model vs data. 
  	"""	
  	(r,l) = (region,layer)

  	modeldata =  self.dataD[(r,l,)]
  	modellat = self.dataD[(r,l,'lat')]
  	modellon = self.dataD[(r,l,'lon')]

  	print "mapplotsRegionsLayers:\t",r,l, "model contains",len(modeldata),'model data'
  	print "mapplotsRegionsLayers:\t",r,l, "model lat:",modellat.min(),modellat.mean(),modellat.max()
  	print "mapplotsRegionsLayers:\t",r,l, "model lon:",modellon.min(),modellon.mean(),modellon.max() 
  	  	
	if self.dataFile:
	  	dnc = Dataset(self.dataFile,'r')
  		datadata =  tst.getHorizontalSlice(dnc,self.datacoords,self.datadetails,layer,data = '').squeeze()
  		datalat = ukp.extractData(dnc,[],key=self.datacoords['lat'])
  		datalon = ukp.extractData(dnc,[],key=self.datacoords['lon']) 
		
	else:
		datadata = np.ma.array([-1000,],mask=[True,])
		datalat  = np.ma.array([-1000,],mask=[True,])
		datalon  = np.ma.array([-1000,],mask=[True,])
  	print "mapplotsRegionsLayers:\t",r,l, "contains",len(datadata),'in situ data'
  	print "mapplotsRegionsLayers:\t",r,l, "data lat:",datalat.min(),datalat.mean(),datalat.max()
  	print "mapplotsRegionsLayers:\t",r,l, "data lon:",datalon.min(),datalon.mean(),datalon.max() 	
	
	titles = [' '.join([self.model,'('+self.jobID+')',str(layer),self.modeldetails['name']]),
		  ' '.join([self.datasource,str(layer),self.datadetails['name']])]
  	tsp.mapPlotPair(modellon, modellat, modeldata,
  			datalon,datalat,datadata,
  			filename,
  			titles	= titles,
  			lon0=0.,drawCbar=True,cbarlabel='',dpi=100,)

	
	
  def makePlots(self):
	if self.debug: print "timeseriesAnalysis:\t makePlots."		  


	#####
	# Hovmoeller plots
	for r in self.regions:
	    for m in self.metrics: 
	   
	   	#####
	   	# Load data layers:

		data = {}
		modeldata = {}
		
	  	for l in self.layers:
	  		if type(l) == type('str'):continue	# no strings, only layers.
			#####
			# Test for presence/absence of in situ data.
			try:	dataslice = self.dataD[(r,l)]	  
			except:	dataslice = []
			try:	dataslice = dataslice.compressed()
			except:	pass
			
			if m == 'mean': 	data[l] = dataslice.mean()
			if m == 'median': 	data[l] = np.median(dataslice)
			
			modeldata[l] = self.modeldataD[(r,l,m)]
			print "loading hov data: model",l, r,m
		#####
		# check that multiple layers were requested.
		if len(data)<1: continue
		

		#####
		# create a dictionary of model depths and layers.
	  	mnc = Dataset(self.modelFiles[-1],'r')		
		modelZcoords = {i:z for i,z in enumerate(mnc.variables[self.modelcoords['z']][:])}
	  	mnc.close()  	

		if self.dataFile:
		  	dnc = Dataset(self.dataFile,'r')	
		  	dataZcoords = {i:z for i,z in enumerate(dnc.variables[self.datacoords['z']][:])}
		  	dnc.close()  	
		else: 	dataZcoords = {}
		
	    	hovfilename = ukp.folder('images/timeseries/'+self.jobID+'/'+self.dataType)+'_'.join(['hov',self.jobID,self.dataType,r,m,])+'.png'
	    	
		title = ' '.join([r,m,self.dataType])
		tsp.hovmoellerPlot(modeldata,data,hovfilename, modelZcoords = modelZcoords, dataZcoords= dataZcoords, title = title,)		
	
	
	
	
	
	
	
	
	
	
	
	
	#####
	# Trafficlight and map plots:
	for r in self.regions:
	  for l in self.layers:

	    #####
	    # Don't make pictures for each integer or float layer, only the ones that are strings. 
	    if type(l) in [type(0),type(0.)]:continue
	    
	    #####
	    # Test for presence/absence of in situ data.
    	    try:	dataslice = self.dataD[(r,l)]	  
    	    except:	dataslice = []
    	    try:	dataslice = dataslice.compressed()
    	    except:	pass
    	    


	    for m in self.metrics:  
	    	modeldataDict = self.modeldataD[(r,l,m)]

		times = sorted(modeldataDict.keys())
		modeldata = [modeldataDict[t] for t in times]
		filename = ukp.folder('images/timeseries/'+self.jobID+'/'+self.dataType)+'_'.join(['trafficlight',self.jobID,self.dataType,r,str(l),m,])+'.png'
		title = ' '.join([r,str(l),m,self.dataType])
		tsp.trafficlightsPlot(times,modeldata,dataslice,metric = m, title = title,filename=filename)		
		
		#filename = ukp.folder('images/timeseries/'+self.jobID)+'_'.join([self.jobID,self.dataType,r,l,m,])+'_trafficlights.png'		
		#tsp.trafficlightsPlots(times,modeldata,dataslice,title = title,filename=filename)
			
    	    mapfilename = ukp.folder('images/timeseries/'+self.jobID+'/'+self.dataType)+'_'.join(['map',self.jobID,self.dataType,str(l),])+'.png'
    	    if ukp.shouldIMakeFile([self.shelvefn,self.shelvefn_insitu],mapfilename,debug=False):
	    	    self.mapplots( l, mapfilename)
	

    	    mapfilename = ukp.folder('images/timeseries/'+self.jobID+'/'+self.dataType)+'_'.join(['map',self.jobID,self.dataType,str(l),r,])+'.png'
    	    if ukp.shouldIMakeFile([self.shelvefn,self.shelvefn_insitu],mapfilename,debug=False):
	    	    self.mapplotsRegionsLayers( r,l, mapfilename)
				
			
			
			
			
			
			
			
			
			
			
			
			
			
