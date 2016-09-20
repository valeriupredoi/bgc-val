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

import numpy as np
from shelve import open as shOpen
from netCDF4 import Dataset,num2date
import os
import shutil

#Specific local code:
import UKESMpython as ukp
from pftnames import getLongName
import timeseriesTools as tst 
import timeseriesPlots as tsp 
from makeEORCAmasks import makeMaskNC
#getTimes, loadData



        


class profileAnalysis:
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
		
	
	if debug: print "profileAnalysis:\t init."	
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

  	self.gridmaskshelve 	= ukp.folder(self.workingDir)+'_'.join([self.grid,])+'_masks.shelve'		
  	self.shelvefn 		= ukp.folder(self.workingDir)+'_'.join(['profile',self.jobID,self.dataType,])+'.shelve'
	self.shelvefn_insitu	= ukp.folder(self.workingDir)+'_'.join(['profile',self.jobID,self.dataType,])+'_insitu.shelve'

	self._masksLoaded_ 	= False
	
	#####
	# Load Data file	
 	self.loadData()
	#assert 0
	
	#####
	# Load Model File
  	self.loadModel()  	

	#####
	# Make the plots:
  	self.makePlots()
  	
  	
  	
  	
  def loadModel(self):
	if self.debug: print "profileAnalysis:\tloadModel."
	####
	# load and calculate the model info
	
	try:
		if self.clean: 
			print "profileAnalysis:\tloadModel:\tUser requested clean run. Wiping old data."
			assert 0		
		sh = shOpen(self.shelvefn)
		readFiles 	= sh['readFiles']
		modeldataD 	= sh['modeldata']
		sh.close()
		print "OprofileAnalysis:\tloadModel:\tpened shelve:", self.shelvefn, '\tread', len(readFiles)
	except:
		readFiles = []
		modeldataD = {}
		for r in self.regions:
		 for l in self.layers:
		  for m in self.metrics:
		   	modeldataD[(r,l,m)] = {}
		   	
		print "profileAnalysis:\tloadModel:\tCould not open shelve:", self.shelvefn, '\tread', len(readFiles)

	###############
	# Check whethere there has been a change in what was requested:
	for r in self.regions:
	  for l in self.layers:
	    for m in self.metrics:
	    	if self.debug:print "profileAnalysis:\tloadModel:\tChecking: ",[r,l,m,],'\t',
	    	try:
	    		if self.debug: print 'has ', len(modeldataD[(r,l,m)].keys()), 'keys'
	    	except: 
	    		readFiles = []
	    		modeldataD[(r,l,m)] = {}
	    		if self.debug: print 'has no keys'
	    	try:	
	    	    	if len(modeldataD[(r,l,m)].keys()) == 0: 
	    	    		print "profileAnalysis:\tloadModel:\tmodeldataD[",(r,l,m),"] has no keys"
	    	    		readFiles = []
	    	    		assert 0
	    	    		
	    	except: pass

	#####
	# Summarise checks
	if self.debug:	
		print "profileAnalysis:\tloadModel:\tloadModel:post checks:"
		#print "modeldataD:",modeldataD
		print "profileAnalysis:\tloadModel:\tshelveFn:",self.shelvefn
		print "profileAnalysis:\tloadModel:\treadFiles:",
		try:	print readFiles[-1]
		except: print '...'

	###############
	# Load files, and calculate fields.
	openedFiles = 0					
	for fn in self.modelFiles:
		if fn in readFiles:continue
		
		if not self._masksLoaded_: 
			self.loadMasks()		
		
		print "loadModel:\tloading new file:",self.dataType,fn,
		nc = Dataset(fn,'r')
		ts = tst.getTimes(nc,self.modelcoords)
		meantime = np.mean(ts)
		print "\ttime:",meantime
		
		#DL = tst.DataLoader(fn,nc,self.modelcoords,self.modeldetails, regions = self.regions, layers = self.layers,)
		nc = Dataset(fn,'r')
		dataAll = ukp.extractData(nc,self.modeldetails).squeeze()
		
		for r in self.regions:
		  for m in self.metrics:
			if m =='mean':
				data = ukp.mameanaxis(np.ma.masked_where((self.modelMasks[r] != 1) + dataAll.mask,dataAll), axis=(1,2))
				
				if self.debug:print "profileAnalysis:\tloadModel.",r,m,self.dataType,'\tyear:',int(meantime), 'mean:',data.mean()
				#if self.debug:print "profileAnalysis:\tloadModel.",self.dataType, data.shape, data.min(),data.max(), dataAll.shape ,self.modelMasks[r].shape, dataAll.min(),dataAll.max()
				
				alllayers = []
				for l,d in enumerate(data):
					#print "Saving model data profile",r,m,l,d
					modeldataD[(r,l,m)][meantime] = d
					alllayers.append(l)
					
				#####
				# Add a masked value in layers where there is no data.
				
				for l in self.layers:
					if l in alllayers:continue
					modeldataD[(r,l,m)][meantime] = np.ma.masked
					
			else:
				print 'ERROR:',m, "not implemented in profile"
				assert 0
								
		readFiles.append(fn)		
		openedFiles+=1			

		nc.close()
		if openedFiles:
			print "Saving shelve:",self.dataType, self.shelvefn, '\tread', len(readFiles)				
			sh = shOpen(self.shelvefn)
			sh['readFiles']		= readFiles
			sh['modeldata'] 	= modeldataD
			sh.close()
			openedFiles=0	
	if openedFiles:
		print "Saving shelve:",self.dataType, self.shelvefn, '\tread', len(readFiles)				
		sh = shOpen(self.shelvefn)
		sh['readFiles']		= readFiles
		sh['modeldata'] 	= modeldataD
		sh.close()
	
	self.modeldataD = modeldataD
	if self.debug: print "profileAnalysis:\tloadModel.\t Model loaded:",	self.modeldataD.keys()[:3], '...', len(self.modeldataD.keys())	

  def loadMasks(self):
  	#####
	# Here we load the masks file.
	self.maskfn = 'data/'+self.grid+'_masks.nc'
	
	if not os.path.exists(self.maskfn):
		print "Making mask file",self.maskfn

		makeMaskNC(self.maskfn, self.regions, self.grid)
	
	self.modelMasks= {}
	
	ncmasks = Dataset(self.maskfn,'r')
	
	for r in self.regions:
		if r in ncmasks.variables.keys():
			print "Loading mask",r
			self.modelMasks[r] = ncmasks.variables[r][:]
			
		else:
			newmask = 'data/'+self.grid+'_masks_'+r+'.nc'
			makeMaskNC(newmask, [r,], self.grid)
			nc = Dataset(newmask,'r')
			self.modelMasks[r] = nc.variables[r][:]
			nc.close()			
			
	print "Loaded masks",self.modelMasks.keys()

	ncmasks.close()
	self._masksLoaded_ = True



	
  def loadData(self):
  	
	if self.debug: print "profileAnalysis:\t loadData.",self.dataFile		
	
  	if not self.dataFile: 
 		if self.debug: print "profileAnalysis:\t No data File provided:",self.dataFile		 		
		self.dataD = {}
		return

  	if not os.path.exists(self.dataFile): 
 		if self.debug: print "profileAnalysis:\tWARNING:\t No such data File:",self.dataFile		 
		self.dataD = {}
		return
				
	###############
	# load and calculate the real data info
	try:
		if self.clean: 
			print "profileAnalysis:\t loadData\tUser requested clean run. Wiping old data."
			assert 0		
		sh = shOpen(self.shelvefn_insitu)
		dataD 	= sh['dataD']
		sh.close()
		print "profileAnalysis:\t loadData\tOpened shelve:", self.shelvefn_insitu
		self.dataD = dataD
	except:
		dataD = {}
		print "profileAnalysis:\t loadData\tCould not open shelve:", self.shelvefn_insitu


	###############
	# Test to find out if we need to load the netcdf, or if we can just return the dict as a self.object.
	needtoLoad = False
	for r in self.regions:
	    #if needtoLoad:continue
	    for l in self.layers:
		#if needtoLoad:continue
	    	try:	
	    		dat = self.dataD[(r,l)]
	    		#test = (len(),self.dataD[(r,l)].shape)
	    		print "profileAnalysis:\t loadData\t",(r,l),dat
	    	except: 
			needtoLoad=True
			print "profileAnalysis:\t loadData\tUnable to load",(r,l)

	if needtoLoad: pass	
	else:
		self.dataD = dataD	
		return
		
	###############
	# Loading data for each region.
	print "profileAnalysis:\t loadData,\tloading ",self.dataFile
	nc = Dataset(self.dataFile,'r')
	data = tst.loadData(nc, self.datadetails)

	
	
	###############
	# Loading data for each region.
	dl = tst.DataLoader(self.dataFile,'',self.datacoords,self.datadetails, regions = self.regions, layers = self.layers,)
	
	for r in self.regions:
	    for l in self.layers:
	    	dataD[(r,l)] = dl.load[(r,l,)]	
	    	dataD[(r,l,'lat')] = dl.load[(r,l,'lat')]		    	
	    	dataD[(r,l,'lon')] = dl.load[(r,l,'lon')]
		if len(dataD[(r,l)])==0  or np.ma.is_masked(dataD[(r,l)]):
			dataD[(r,l)]  = np.ma.masked
			dataD[(r,l,'lat')]  = np.ma.masked
			dataD[(r,l,'lon')]  = np.ma.masked
									    	
    		print "profileAnalysis:\t loadData,\tloading ",(r,l),  dataD[(r,l)].min(),  dataD[(r,l)].max(),  dataD[(r,l)].mean()
    		
	###############
	# Savng shelve		
	print "profileAnalysis:\t loadData.\tSaving shelve:", self.shelvefn_insitu			
	try:
		sh = shOpen(self.shelvefn_insitu)
		sh['dataD'] 	= dataD
		sh.close()
	except:
		print "profileAnalysis:\t WARNING.\tSaving shelve failed, trying again.:", self.shelvefn_insitu			
		print "Data is", dataD
		shutil.move(self.shelvefn_insitu, self.shelvefn_insitu+'.broken')
		sh = shOpen(self.shelvefn_insitu)
		sh['dataD'] 	= dataD
		sh.close()		
	 	
	self.dataD = dataD


 	

  def mapplotsRegionsLayers(self,):
  
  	"""	Makes a map plot of model vs data for each string-named layer (not numbered layers). 
  	"""
  	newlayers = [l for l in self.layers if type(l) not in [type(0),type(0.) ]]
	mDL = tst.DataLoader(self.modelFiles[-1],'',self.modelcoords,self.modeldetails, regions = self.regions, layers = newlayers,)
	for r in self.regions:
	    for l in self.layers:	
		if type(l) in [type(0),type(0.)]:continue
 		mapfilename = ukp.folder(self.imageDir+'/'+self.dataType)+'_'.join(['map',self.jobID,self.dataType,str(l),r,])+'.png'
   		modeldata	= mDL.load[(r,l)]
   		modellat	= mDL.load[(r,l,'lat')]
   		modellon	= mDL.load[(r,l,'lon')]
		  	
		if not len(modeldata): continue
		
	  	print "mapplotsRegionsLayers:\t",r,l, "model contains",len(modeldata),'model data'
	  	print "mapplotsRegionsLayers:\t",r,l, "model lat:",modellat.min(),modellat.mean(),modellat.max()
	  	print "mapplotsRegionsLayers:\t",r,l, "model lon:",modellon.min(),modellon.mean(),modellon.max() 
  	  	
		if self.dataFile:
		    	datadata	= self.dataD[(r,l)] 
		    	datalat		= self.dataD[(r,l,'lat')]
		    	datalon		= self.dataD[(r,l,'lon')]
		    	
		else:
			datadata = np.ma.array([-1000,],mask=[True,])
			datalat  = np.ma.array([-1000,],mask=[True,])
			datalon  = np.ma.array([-1000,],mask=[True,])

	  	print "mapplotsRegionsLayers:\t",r,l, "contains",len(datadata),'in situ data'
	  	print "mapplotsRegionsLayers:\t",r,l, "data lat:",len(datalat),datalat.min(),datalat.mean(),datalat.max()
	  	print "mapplotsRegionsLayers:\t",r,l, "data lon:",len(datalon),datalon.min(),datalon.mean(),datalon.max() 	
	
		titles = [' '.join([getLongName(t) for t in [self.model,'('+self.jobID+')',str(l),self.modeldetails['name']]]),
			  ' '.join([getLongName(t) for t in [self.datasource,str(l),self.datadetails['name']]])]
			  
	  	tsp.mapPlotPair(modellon, modellat, modeldata,
	  			datalon,datalat,datadata,
	  			mapfilename,
	  			titles	= titles,
	  			lon0=0.,drawCbar=True,cbarlabel='',dpi=100,)

	
	
  def makePlots(self):
	if self.debug: print "profileAnalysis:\t makePlots."	  


	#####
	# create a dictionary of model and data depths and layers.
  	mnc = Dataset(self.modelFiles[-1],'r')		
	modelZcoords = {i:z for i,z in enumerate(mnc.variables[self.modelcoords['z']][:])}
  	mnc.close()  

	if self.dataFile:
	  	dnc = Dataset(self.dataFile,'r')	
	  	dataZcoords = {i:z for i,z in enumerate(dnc.variables[self.datacoords['z']][:])}
	  	print 
	  	dnc.close()  	
	else: 	dataZcoords = {}

	#####
	# Hovmoeller plots
	for r in self.regions:
	    for m in self.metrics: 
	    	if m not in ['mean','median','min','max',]:continue

	   	#####
	   	# Load data layers:
		data = {}
		if self.dataFile:		
	  	    for l in self.layers:
	  		#print "Hovmoeller plots:",r,m,l
	  		
	  		if type(l) == type('str'):continue	# no strings, only numbered layers.
	  		if l > max(dataZcoords.keys()): continue
			#####
			# Test for presence/absence of in situ data.
			try:	
				dataslice = self.dataD[(r,l)]	  
				dataslice = dataslice.compressed()				
			except:	dataslice = np.ma.array([-1000,],mask=[True,])
		
			if m == 'mean': 
				try:	data[l] = np.ma.mean(dataslice)
				except:	data[l] = np.ma.array([-1000,],mask=[True,])				
			elif m == 'median':
				try: 	data[l] = np.ma.median(dataslice)
				except:	data[l] = np.ma.array([-1000,],mask=[True,])
			elif m == 'min': 
				try: 	data[l] = np.ma.min(dataslice)
				except:	data[l] = np.ma.array([-1000,],mask=[True,])
			elif m == 'max':
				try:	data[l] = np.ma.max(dataslice)
				except:	data[l] = np.ma.array([-1000,],mask=[True,])				
			
			#print "makePlots:\tHovmoeller plots:",r,m,l,'\tdata'#,data[l]								

	   	#####
	   	# Load model layers:
		modeldata = {}	   	
	  	for l in self.layers:
	  		if type(l) == type('str'):continue	# no strings, only numbered layers.
	  		if l > max(modelZcoords.keys()): continue
			modeldata[l] = self.modeldataD[(r,l,m)]
			
		#####
		# check that multiple layers were requested.
		#if len(data.keys())<1: continue
		if len(modeldata.keys())<1: continue
	


		title = ' '.join([getLongName(t) for t in [r,m,self.dataType]])	
	    	profilefn = ukp.folder(self.imageDir+'/'+self.dataType)+'_'.join(['profile',self.jobID,self.dataType,r,m,])+'.png'
	    	axislabel = getLongName(self.modeldetails['name'])+', '+getLongName(self.modeldetails['units'])
		if  ukp.shouldIMakeFile([self.shelvefn, self.shelvefn_insitu],profilefn,debug=False):					    	
			tsp.profilePlot(modeldata,data,profilefn, modelZcoords = modelZcoords, dataZcoords= dataZcoords, xaxislabel = axislabel,title = title,)			
			
			
	    	hovfilename = ukp.folder(self.imageDir+'/'+self.dataType)+'_'.join(['profilehov',self.jobID,self.dataType,r,m,])+'.png'
		if  ukp.shouldIMakeFile([self.shelvefn, self.shelvefn_insitu],hovfilename,debug=False):				
			tsp.hovmoellerPlot(modeldata,data,hovfilename, modelZcoords = modelZcoords, dataZcoords= dataZcoords, title = title,zaxislabel =axislabel, diff=False)		
	
	    	hovfilename_diff = ukp.folder(self.imageDir+'/'+self.dataType)+'_'.join(['profileDiff',self.jobID,self.dataType,r,m,])+'.png'
		if  ukp.shouldIMakeFile([self.shelvefn, self.shelvefn_insitu],hovfilename_diff,debug=False):					    	
			tsp.hovmoellerPlot(modeldata,data,hovfilename_diff, modelZcoords = modelZcoords, dataZcoords= dataZcoords, title = title,zaxislabel =axislabel,diff=True)		
	


			
			
			
			
			
			
