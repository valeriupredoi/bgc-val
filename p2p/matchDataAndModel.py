#!/usr/bin/ipython 
#
# Copyright 2014, Plymouth Marine Laboratory
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
#

from sys import argv
from os.path import exists,split, getmtime, basename
from glob import glob
from shelve import open as shOpen
from shutil import copy2
from math import radians, cos, sin, asin, sqrt
from netCDF4 import num2date,Dataset
from datetime import datetime
import numpy as np

######
# local imports
import UKESMpython as ukp 
from pftnames import getmt, CMIP5models


#####	
# These are availalble in the module:
#	https://gitlab.ecosystem-modelling.pml.ac.uk/ledm/netcdf_manip
from pruneNC import pruneNC
from convertToOneDNC import convertToOneDNC
from mergeNC import mergeNC
from changeNC import changeNC, AutoVivification


#TO DO
#	This still requires the netcdf_manip library, the ORCA1bathy file

class matchDataAndModel:
  """	matchDataAndModel: 
  	This code takes the real data from in situ  measurments in netcdf format and the Model data and created two 1D matched netcdfs. 
	The 1D matched netcdfs are then used to make plots and perform statistical analysis (not in this code).
	The first step is to produce lightweight "pruned" versions of the files, which have the unused fields stripped out.
	Some of the datasets are too large to run this code on desktop machine, so in those cases we request a specific depthLevel, ie "Surface".
	Debug: prints more statements.
  """


  def __init__(self,DataFile,ModelFile,dataType, workingDir = '',DataVars='',ModelVars='',  model = '',jobID='', year='clim',depthLevel='', grid='ORCA1',gridFile='',debug = True,):

	if debug:
		print "matchDataAndModel:\tINFO:\tStarting matchDataAndModel"
		print "matchDataAndModel:\tINFO:\tData file:  \t",DataFile
		print "matchDataAndModel:\tINFO:\tModel file: \t",ModelFile
		print "matchDataAndModel:\tINFO:\tData Type:  \t",dataType
			
	self.DataFile=DataFile
	self.ModelFile = ModelFile 
		
	self.DataVars=DataVars	
	self.ModelVars=ModelVars	

	self.model = model 
	self.jobID = jobID 
	self.year = year
	self.depthLevel = depthLevel
	self.debug = debug	
		
	self.dataType = dataType
	self._meshLoaded_ = False
	
	if debug: print  "matchDataAndModel:\tINFO:\t",self.dataType, '\tModelfile:', self.ModelFile
		
	self.compType= 'MaredatMatched-'+self.model+'-'+self.jobID+'-'+self.year
		
	if workingDir =='':
		self.workingDir = ukp.folder('/data/euryale7/scratch/ledm/ukesm_postProcessed/ukesm/outNetCDF/'+'/'.join([self.compType,self.dataType+self.depthLevel]) )
	else: 	self.workingDir = workingDir	
	self.grid = grid
	if gridFile =='':
		self.gridFile = ukp.getGridFile(grid)
	else:	self.gridFile = gridFile
	print "matchDataAndModel:\tINFO:\tGrid:  \t",grid			
	print "matchDataAndModel:\tINFO:\tGrid File:  \t",gridFile

	#if grid.upper() in ['ORCA1',]:
	#	self.grid = 'ORCA1'
	#	self.gridFile    = "data/mesh_mask_ORCA1_75.nc"
	#if grid.upper() in ['ORCA025',]:	
	#	self.grid = 'ORCA025'
		#####
		# Please add files to link to 
	#	for orcafn in [ "/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/mesh_mask_ORCA025_75.nc",	# PML
	#			"/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_ORCA025_75.nc",]:	# JASMIN
	#		if exists(orcafn):	self.gridFile  = orcafn
	#	
	#	try: 
	#		if exists(self.gridFile):pass
	#	except: 
	#		print "matchDataAndModel:\tERROR:\tIt's not possible to load the ORCA025 grid on this machine. Please add the ORCA025 file to the orcafn list to p2p/matchDataAndModel.py"
	#		assert False
	#if grid in ['Flat1deg',]:	
	#	self.grid = 'Flat1deg'
	#	self.gridFile = 'data/Flat1deg.nc'
		
		#####
		# Please add files to link to 
		#for orcafn in [ "/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/mesh_mask_ORCA025_75.nc",	# PML
		#		"/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_ORCA025_75.nc",]:	# JASMIN
		#	if exists(orcafn):	self.gridFile  = orcafn
		
		#try: 
		#	if exists(self.gridFile):pass
		#except: 
		#	print "matchDataAndModel:\tERROR:\tIt's not possible to load the ORCA025 grid on this machine. Please add the ORCA025 file to the orcafn list to p2p/matchDataAndModel.py"
		#	assert False
						
		
	self.matchedShelve 	= ukp.folder(self.workingDir)+self.model+'-'+self.jobID+'_'+self.year+'_'+'_'+self.dataType+'_'+self.depthLevel+'_matched.shelve'
	self.matchesShelve 	= ukp.folder(['shelves','ModelMatches',])+model+jobID+year+self.dataType+self.depthLevel+self.grid+'.shelve'

	self.workingDirTmp = 	ukp.folder(self.workingDir+'tmp')
	self.DataFilePruned=	self.workingDirTmp+'Data_' +self.dataType+'_'+self.depthLevel+'_'+self.model+'-'+self.jobID+'-'+self.year+'_pruned.nc'
	self.ModelFilePruned=	self.workingDirTmp+'Model_'+self.dataType+'_'+self.depthLevel+'_'+self.model+'-'+self.jobID+'-'+self.year+'_pruned.nc'	
	
	self.DataFile1D  	= self.workingDirTmp + basename(self.DataFilePruned).replace('pruned.nc','1D.nc') 
	self.maskedData1D	= self.workingDir    + basename(self.DataFile1D)
	self.Model1D     	= self.workingDir    + basename(self.ModelFilePruned).replace('pruned.nc','1D.nc')

	self.MatchedModelFile = self.Model1D
	self.MatchedDataFile  = self.maskedData1D
	self.run()


	
  def run(self,):
	"""There are two methods written for manipulating data.
	   One is designed to work with WOA formats, the other with MAREDAT formats.
	   Other data formats are run manually.
	"""
	if not ukp.shouldIMakeFile(self.DataFile,self.MatchedDataFile,debug=False) and not ukp.shouldIMakeFile(self.ModelFile,self.MatchedModelFile,debug=False):
		print "matchDataAndModel:\trun:\talready created:\t",self.maskedData1D, '\n\t\t\tand\t',self.Model1D
		return
	

	self._pruneModelAndData_()	
	self._convertDataTo1D_()	
	self._matchModelToData_()
	self._convertModelToOneD_()
	self._applyMaskToData_()




  def _pruneModelAndData_(self,):
   	""" This routine reduces the full 3d netcdfs by pruning the unwanted fields.
  	""" 
  	
	if ukp.shouldIMakeFile(self.ModelFile,self.ModelFilePruned,debug=False):
		print "matchDataAndModel:\tpruneModelAndData:\tMaking ModelFilePruned:", self.ModelFilePruned
		p = pruneNC(self.ModelFile,self.ModelFilePruned,self.ModelVars, debug = self.debug) 	
	else:	
		print "matchDataAndModel:\tpruneModelAndData:\tModelFilePruned already exists:",self.ModelFilePruned
  	 
	if ukp.shouldIMakeFile(self.DataFile,self.DataFilePruned,debug=False):
		print "matchDataAndModel:\tpruneModelAndData:\tMaking DataFilePruned:", self.DataFilePruned
		p = pruneNC(self.DataFile,self.DataFilePruned,self.DataVars, debug = self.debug) 	
	else:	
		print "matchDataAndModel:\tpruneModelAndData:\tDataFilePruned already exists:",self.DataFilePruned
	 



  def _convertDataTo1D_(self,):
   	""" This routine reduces the In Situ data into a 1D array of data with its lat,lon,depth and time components.
  	"""
		
	if not ukp.shouldIMakeFile(self.DataFilePruned,self.DataFile1D,debug=False):
		print "matchDataAndModel:\tconvertDataTo1D:\talready exists: (DataFile1D):\t",self.DataFile1D
		return

	print "matchDataAndModel:\tconvertDataTo1D:\topening DataFilePruned:\t",self.DataFilePruned		
	#nc = ncdfView(self.DataFilePruned,Quiet=True)
	nc = Dataset(self.DataFilePruned,'r')
		
	#WOADatas = [a+self.depthLevel for a in ['salinity','temperature','temp','sal','nitrate','phosphate','silicate',]]	   	 
	
	if self.depthLevel in ['Surface','100m','200m','500m','1000m','2000m','Transect',]:	
	    if nc.variables[self.DataVars[0]].shape in [(12, 14, 180, 360), (12, 24, 180, 360)]: # WOA format
		mmask = np.ones(nc.variables[self.DataVars[0]].shape)
		#####
		# This could be rewritten to just figure out which level is closest, but a bit much effort for not a lot of gain.
		
		if self.depthLevel in ['Surface','200m','100m','500m','1000m',]: 
			if self.depthLevel in ['Surface',]:	k = 0
			if self.depthLevel == '100m': 	k = 6			
			if self.depthLevel == '200m': 	k = 9
			if self.depthLevel == '500m': 	k = 13
			if self.depthLevel == '1000m': 	k = 18					
			mmask[:,k,:,:] = 0
						
		if self.depthLevel == 'Transect':	mmask[:,:,:,332] = 0   # Pacific Transect.	
		if self.depthLevel == 'PTransect':	mmask[:,:,:,200] = 0   # Pacific Transect.			
		mmask +=nc.variables[self.DataVars[0]][:].mask
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking mask shape:',mmask.shape		
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking WOA style flat array:',self.DataFilePruned,'-->',self.DataFile1D	
	  	convertToOneDNC(self.DataFilePruned,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)		  	
	    elif nc.variables[self.DataVars[0]].shape in [(1, 102, 180, 360),]: # WOA Oxygen
		mmask = np.ones(nc.variables[self.DataVars[0]].shape)
		#####
		# This could be rewritten to just figure out which level is closest, but a bit much effort for not a lot of gain.
		
		if self.depthLevel in ['Surface','200m','100m','500m','1000m','2000m']: 
			if self.depthLevel in ['Surface',]:	k = 0
			#if self.depthLevel == '100m': 	k = 6			
			#if self.depthLevel == '200m': 	k = 9
			#if self.depthLevel == '500m': 	k = 13
			#if self.depthLevel == '1000m': k = 18					
			if self.depthLevel == '2000m': 	k = 46
			mmask[:,k,:,:] = 0
						
		if self.depthLevel == 'Transect':	mmask[:,:,:,332] = 0   # Pacific Transect.	
		if self.depthLevel == 'PTransect':	mmask[:,:,:,200] = 0   # Pacific Transect.			
		mmask +=nc.variables[self.DataVars[0]][:].mask
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking mask shape:',mmask.shape		
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking WOA style flat array:',self.DataFilePruned,'-->',self.DataFile1D	
	  	convertToOneDNC(self.DataFilePruned,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)		  	
	  	
	  		  	
	    elif nc.variables[self.DataVars[0]].shape in [(12, 33, 180, 360), (12,1,180,360), ]: # Chl format
		mmask = np.ones(nc.variables[self.DataVars[0]].shape)
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking Chl-style flat array:',self.DataFilePruned,'-->',self.DataFile1D			
		if self.depthLevel in ['Surface',]:	
			mmask[:,0,:,:] = 0
		else: 
			print "matchDataAndModel:\tERROR:\t Depth level not recognises. (12, 33, 180, 360), (12,1,180,360)" ,self.depthLevel		
			assert False
		if self.depthLevel == 'Transect':assert 0
		mmask +=nc.variables[self.DataVars[0]][:].mask
		
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking mask shape:',mmask.shape		
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking Chl style flat array:',self.DataFilePruned,'-->',self.DataFile1D	
	  	convertToOneDNC(self.DataFilePruned,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)
	  	
	  	
	    elif nc.variables[self.DataVars[0]].shape in [(12,57,180, 360), ]: # O2 format
		mmask = np.ones(nc.variables[self.DataVars[0]].shape)
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking Oxygen-style flat array:',self.DataFilePruned,'-->',self.DataFile1D			
		if self.depthLevel in ['Surface','200m','100m','500m','1000m',]: 		
			if self.depthLevel in ['Surface',]:	k = 0
			elif self.depthLevel == '100m': 	k = 20			
			elif self.depthLevel == '200m': 	k = 24
			elif self.depthLevel == '500m': 	k = 36
			elif self.depthLevel == '1000m': 	k = 46							
			else:
				print "matchDataAndModel:\tERROR:\t Depth level not recognises. (12,57,180, 360)" ,self.depthLevel
				assert False
			mmask[:,k,:,:] = 0		
		if self.depthLevel == 'Transect':
			mmask[:,:,:,155]  =0 
		mmask +=nc.variables[self.DataVars[0]][:].mask				
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking mask shape:',mmask.shape		
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking Chl style flat array:',self.DataFilePruned,'-->',self.DataFile1D	
	  	convertToOneDNC(self.DataFilePruned,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)
	  		  	
	    elif nc.variables[self.DataVars[0]].ndim ==1:
	    	# This kind of file is already 1D, but we wantto apply a depth cut to it?
	    	mmask = np.ones(nc.variables[self.DataVars[0]].shape)
	    	for i,dep in enumerate(np.abs(nc.variables['DEPTH'][:])):
	    		if dep >20:continue # assume surface layer of 20?
			mmask[i] = 0
		mmask +=nc.variables[self.DataVars[0]][:].mask	
		convertToOneDNC(self.DataFilePruned,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)			    		
	    else:
	    
		print 'matchDataAndModel:\tconvertDataTo1D:\tYou need to add more file spcific depthLevels here.', nc.variables[self.DataVars[0]].shape, self.DataFilePruned,'-->',self.DataFile1D
		assert False
			
	else:
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking',self.DataFilePruned,'-->',self.DataFile1D
	  	if len(self.DataVars):	convertToOneDNC(self.DataFilePruned, self.DataFile1D, debug=True, variables = self.DataVars)
	  	else:			convertToOneDNC(self.DataFilePruned, self.DataFile1D, debug=True)
  	nc.close()	
  	
	
	
		
  	
  def _matchModelToData_(self,):
  	print "matchModelToData:\tOpened MAREDAT netcdf:", self.DataFile1D
  	
  	ncIS = Dataset(self.DataFile1D,'r')
  	#ncIS = ncdfView(self.DataFile1D,Quiet=True)  	
	is_i	= ncIS.variables['index'][:]
	
	try:
		s = shOpen(self.matchedShelve)
		maxIndex = s['maxIndex']
		self.maremask = s['maremask']
		self.matches = s['matches']
		self.imatches = s['imatches']		
		s.close()		
		print "matchModelToData:\tOpened shelve:", self.matchedShelve
		print "matchModelToData:\tStarting from maxindex:",maxIndex," and ",len(self.matches), " already matched. Mask:",self.maremask.sum()
	except:
		self.matches = {}
		self.imatches = {}
		maxIndex = 0
		self.maremask = np.zeros(is_i.shape) # zero array same length as in situ data.

		print "matchModelToData:\tStarting from maxindex",maxIndex,"\tfinished:",len(self.matches), " already matched. Mask:",self.maremask.sum()
		print "matchModelToData:\tCreating shelve:", self.matchedShelve

	try:
		s = shOpen(self.matchesShelve)		
		lldict  = s['lldict']
		s.close()
	except:
		lldict={}
	finds = 0	
	mt = getmt()
	
	#####
	# Figure out which type of data this is.
	ytype = []
	Models = [m.upper() for m in ['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10','NEMO','IMARNET','CMIP5',]] # skip these to find in situ data types.
	Models.extend(['IMARNET_' +m.upper() for m in ['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10','NEMO',]])	
	Models.extend(['CMIP5_' +m.upper() for m in CMIP5models])		
	for key in mt.keys():
		#key = key.upper()
		if key.upper() in Models:continue
		try:
			if self.dataType in mt[key].keys() and key not in ytype:
				ytype.append(key)
		except:pass
	if len(ytype) == 1:
		ytype = ytype[0]		
	else:
		print "matchModelToData:\tUnable to determine in situ data dataset type (ie, Maredat, WOA, Takahashi etc...)", ytype , (self.dataType)
		print "matchModelToData:\tYou need to add the new data dataset type informationg to getmt() in pftnames.py"		
		print "matchModelToData:\tor remove it from the list of options for ytype"
		assert False

	#####
	# Check if there is any data left to match. 
	# This makes it easier to stop and start the longer analyses.
	if maxIndex+1 >=len(is_i):
		ncIS.close()
		print "matchModelToData:\tNo need to do further matches, Finsished with ",maxIndex+1,"\tfinished:",len(self.matches)
		return 		
		
	#if maxIndex+1 <len(is_i):
	zdict={}
	tdict={}
	print 'mt[',ytype,']:', mt[ytype]
  	is_t	= ncIS.variables[mt[ytype]['t']][:]
  	is_z 	= ncIS.variables[mt[ytype]['z']][:]
  	is_la	= ncIS.variables[mt[ytype]['lat']][:]
	is_lo 	= ncIS.variables[mt[ytype]['lon']][:]	
	tdict   = mt[ytype]['tdict'] 	
	#tdict   = {i:i for i in xrange(12)}
	ncIS.close()	     

	print "tdict:", tdict
	#####
	# This list could be removed by adding a check dimensionality of data after it was been pruned.	
	flatDataOnly =  ['pCO2','seawifs','Seawifs','mld_DT02', 'mld_DR003','mld_DReqDTm02','mld',
			  'dms_and','dms_ara','dms_hal','dms_sim',]
	for d in ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim',]:
	  for i in ['','1','2']:
		flatDataOnly.append(d+i)
	if self.dataType in flatDataOnly: 
  	   	is_z 	= np.ma.zeros(len(is_t))[:]
	   	zdict = {0:0, 0.:0}  

	if not self._meshLoaded_:self.loadMesh()
	    
	for i,ii in enumerate(is_i[maxIndex:]):
		i+=maxIndex
		wt  = is_t[i]
		wz  = is_z[i]
		wla = is_la[i]
		wlo = is_lo[i] 

		#####		
		# Match Latitude and Longitude
		try:
			la,lo = lldict[(wla,wlo)]
		except:
			la,lo = self.getOrcaIndexCC(wla,wlo,debug=False)
			lldict[(wla,wlo)] = la,lo
			finds+=1
			if la == lo == -1:
				print "STRICT ERROR: Could not find, ",wla,wlo
				assert False
				return
			if self.debug:
			    print "matchModelToData:\t",i,'New match:\tlon:',[wlo,self.loncc[la,lo]],'\tlat:',[wla,self.latcc[la,lo]],[finds,len(lldict)]

		#####
		#Match Depth
		try:
			z = zdict[wz]
		except:	 
			z = getORCAdepth(wz,self.depthcc,debug=True)
			zdict[wz]	= z
			if self.debug: print "matchModelToData:\t",i, 'Found new depth:',wz,'m-->',self.depthcc[z], ['z=',z]

		#####			
		#Match Time	
		try:
			t = tdict[wt]	
		except:
			print "matchModelToData:\tunable to find time match in pftnames, mt[",ytype,"]['tdict']", wt
			print "tdict:",tdict
			assert False 
			
			t = getMonthFromSecs(wt)
			tdict[wt] = t
			if self.debug:	print "matchModelToData:\t",i, 'Found new month:', wt, '-->',t

		#####
		# Add match into array
		try:
			tmp = self.matches[(t,z,la,lo)][0]
			self.matches[(t,z,la,lo)].append(i)
			self.imatches[i] = (t,z,la,lo)
			self.maremask[i] = tmp
			#if self.debug: print "matchModelToData:\tWARNING:",i,[wt,wz,wla,wlo], '-->',(t,z,la,lo),'already matched', self.matches[(t,z,la,lo)]
				#for a in self.matches[(t,z,la,lo)]:
					#print '\t',i, a, self.imatches[a]
		except:
			# if this location in the model grid (t,z,la,lo) has not yet been found,
			# self.matches gets a list of all the in situ points that match that location.
			# Conversely, self.maremask's i-th value location of the first time it was found.
			self.matches[(t,z,la,lo)]=[i,]
			self.imatches[i] = (t,z,la,lo)			
			self.maremask[i] = i
			#if self.debug: print "matchModelToData:\tfirst match:",i,[wt,wz,wla,wlo], '-->',(t,z,la,lo)
		
		#####
		# test match up:
		fail=0
		if abs(wz-self.depthcc[z]) >500.:
			print 'depth DOESNT MATCH:',wz,self.depthcc(z)
			fail+=1
		if abs(wla-self.latcc[la,lo]) >2.:
			print 'Latitude DOESNT MATCH:',wla, self.latcc[la,lo]
			fail+=1		
		if abs(wlo-self.loncc[la,lo]) >2.:
			print 'Longitude DOESNT MATCH:',wlo,self.loncc[la,lo]
			fail+=1						
		if fail>0:assert False
		
		
		#####			
		#increment by 1 to save/ end, as it has finished, but i is used as a starting point.
		i+=1
		if i%10000==0:
			if self.debug: print "matchModelToData:\t", i,ii,self.dataType,self.depthLevel,':\t',[wt,wz,wla,wlo] ,'--->',[t,z,la,lo]
			
		if i>1 and i%50000000==0:
			if self.debug: print "matchModelToData:\tSaving Shelve: ", self.matchedShelve
			s = shOpen(self.matchedShelve)
			s['matches']  = self.matches
			s['maxIndex'] = i
			s['maremask'] = self.maremask	
			s['imatches']  = self.imatches
			s.close()
			s = shOpen(self.matchesShelve)		
			s['lldict'] = lldict 
			s.close()				
	
	#assert False
		
	print "matchDataAndModel:\tSaving Shelve", self.matchedShelve	
	s = shOpen(self.matchedShelve)
	s['matches']  = self.matches
	s['imatches']  = self.imatches	
	s['maxIndex'] = i
	s['maremask'] = self.maremask
	s.close()
	    
	s = shOpen(self.matchesShelve)		
	s['lldict'] = lldict 
	s.close()
	print "matchModelToData:\tFinsished with ",maxIndex+1,"\tfinished:",len(self.matches) #, "Mask:",self.maremask.sum()
	
	
  def _convertModelToOneD_(self,):
	if not ukp.shouldIMakeFile(self.ModelFilePruned,self.Model1D,debug=True):
		print "convertModelToOneD:\tconvertModelToOneD:\talready exists:",self.Model1D
		return	
	
	print "convertModelToOneD:\tconvertModelToOneD:\tMaking 1D Model file:", self.ModelFilePruned,'-->', self.Model1D
  	convertToOneDNC(self.ModelFilePruned, self.Model1D,newMask='',debug=self.debug,dictToKeep=self.matches)



	



  def _applyMaskToData_(self,):
  	""" This routine applies the mask of the model to the data. 
  	    It is needed because there are some points where two WOA grid cells fall into the same ORCA1 grid, these excess points should be meaned/medianed.
  	    Similarly, some data points fall into a masked grid cell in the model and need to be masked in the data.
  	"""
	
	if not ukp.shouldIMakeFile(self.ModelFilePruned,self.maskedData1D,debug=True):
		print "applyMaskToData:\tapplyMaskToData:\t", "already exists:",self.maskedData1D
		return	
		
	maremask = self.maremask #np.array([int(a) for a in self.maremask] )
		# zero shouldn't happen
		# some number:
			# need to take the median value of for everytime that the same number appears.
			
	maremaskDict={}
	for i,m in enumerate(maremask):
		try:	maremaskDict[m].append(i)
		except:	maremaskDict[m] = [i,]
	def getMedianVal(arr):
		print "applyMaskToData:\tgetMedianVal: 1d:",arr.shape,maremask.shape,'-->',
		out = []
		
		for m in sorted(maremaskDict.keys()):
		    out.append(np.median([arr[i] for i in maremaskDict[m] ]))	
		arr= np.array(out).squeeze()
		print arr.shape
		return arr

	def getMeanVal(arr):
		print "applyMaskToData:\tgetMeanVal: 1d:",arr.shape,maremask.shape,'-->',
		out = []
		
		for m in sorted(maremaskDict.keys()):
		    out.append(np.mean([arr[i] for i in maremaskDict[m] ]))	
		arr= np.array(out).squeeze()
		print arr.shape
		return arr

	def getFirstVal(arr):
		print "\napplyMaskToData:\tgetFirstVal: 1d:",arr.shape,maremask.shape,len(maremaskDict.keys()),'-->',
		out = []
		for m in sorted(maremaskDict.keys()):
		    out.append(arr[maremaskDict[m]][0])	
		    
		arr= np.array(out).squeeze()
		print arr.shape
		return arr
						
	def applyMask(arr):
		print "applyMask: 1d:",arr.shape,maremask.shape,'-->',
		out = []
		for i,m in enumerate(maremask):
			if not m: continue
			out.append(arr[i])	
		arr= np.array(out).squeeze()
		print arr.shape
		return arr
		

  	ncIS = Dataset(self.DataFile1D,'r')		
  	#ncIS = ncdfView(self.DataFile1D,Quiet=True)		
	av = AutoVivification()
	for v in ncIS.variables.keys():
		if ncIS.variables[v].ndim != 1: 
			if self.debug: print "matchDataAndModel:\tapplyMaskToData:\tERROR:\tthis is suppoed to be the one D file"
			assert False
		print  "matchDataAndModel:\tapplyMaskToData:AutoViv:", v ,len(ncIS.variables[v][:]), len(self.maremask)
		if len(ncIS.variables[v][:])== len(self.maremask):
			if self.debug: print  "matchDataAndModel:\tapplyMaskToData:AutoViv:", v ,'is getting a mask.'
			av[v]['convert'] = getMedianVal
			#av[v]['convert'] = getMeanVal			
			#av[v]['convert'] = getFirstVal			
						
	ncIS.close()
	
	if self.debug: 
		print "matchDataAndModel:\tapplyMaskToData:\tNEW MASK for 1D Maredat:",maremask.shape, maremask.sum(), maremask.size, maremask.size-maremask.sum()
		print "matchDataAndModel:\tapplyMaskToData:\tMaking 1D Maredat file:", self.DataFile1D,'-->', self.maskedData1D
		
	c = changeNC(self.DataFile1D, self.maskedData1D,av,debug = self.debug)
	

  def loadMesh(self,):
      	# This won't work unless its the 1 degree grid.
      	#f self.ORCA == "ORCA1":
 	print "matchModelToData:\tOpened Model netcdf mesh.", self.grid, '(',self.gridFile,')'
  	ncER = Dataset(self.gridFile,'r')
  	
  	#ncER = ncdfView("data/mesh_mask_ORCA1_75.nc",Quiet=True)
  	if 'nav_lat' in ncER.variables.keys():	self.latcc    = ncER.variables['nav_lat'][:].squeeze()
	elif 'lat' in ncER.variables.keys():	self.latcc    = ncER.variables['lat'][:].squeeze()

  	if 'nav_lon' in ncER.variables.keys():	self.loncc    = ncER.variables['nav_lon'][:].squeeze()
	elif 'lon' in ncER.variables.keys():	self.loncc    = ncER.variables['lon'][:].squeeze()

  	if 'deptht' in ncER.variables.keys():	self.depthcc  = ncER.variables['deptht'][:].squeeze()
	elif 'gdept_0' in ncER.variables.keys():self.depthcc  = ncER.variables['gdept_0'][:].squeeze()
	elif 'lev' in ncER.variables.keys():	self.depthcc  = ncER.variables['lev'][:].squeeze()	
	
	if self.loncc.ndim ==1 and  self.loncc.shape!= self.latcc.shape:
		self.loncc,self.latcc = np.meshgrid(self.loncc,self.latcc)
			
	ncER.close()
 	print "matchModelToData:\tloaded mesh.", self.grid, 'lat:',self.latcc.shape, 'lon:',self.loncc.shape,'depth:',self.depthcc.shape
	self._meshLoaded_ = 1
	
  def getOrcaIndexCC(self,lat,lon,debug=True,slowMethod=False,llrange=5.):
	""" takes a lat and long coordinate, an returns the position of the closest coordinate in the NemoERSEM grid.
	    uses the bathymetry file.
	"""
	km = 10.E20
	la_ind, lo_ind = -1,-1
	rangeCutoff=2.
	lat = ukp.makeLatSafe(lat)
	lon = ukp.makeLonSafe(lon)	
	
	if not self._meshLoaded_:self.loadMesh()
	c = (self.latcc - lat)**2 + (self.loncc - lon)**2

	(la_ind,lo_ind) =  np.unravel_index(c.argmin(),c.shape)

	#km2 = abs(haversine((lon, lat), (locc,lacc)))
	if debug: print 'location ', [la_ind,lo_ind],'(',self.latcc[la_ind,lo_ind],self.loncc[la_ind,lo_ind],') is closest to:',[lat,lon]
	if abs(self.latcc[la_ind,lo_ind] -lat ) > rangeCutoff and abs(self.loncc[la_ind,lo_ind] -lon ) > rangeCutoff:
		d = myhaversine(self.loncc[la_ind,lo_ind],self.latcc[la_ind,lo_ind],lon,lat)
		print "getOrcaIndexCC:\tERROR:\tDISTANCE IS TOO FAR:",[self.latcc[la_ind,lo_ind],'-->',  lat,], [self.loncc[la_ind,lo_ind], '-->', lon],'distance:',d
		if d>300.:
			#return -1,-1
			print "distance too great:",d
			assert False
	return la_ind,lo_ind
		

#########################################
# Trimming masked values out of 1d files:


	
#########################################
# Coords and Depth:
def myhaversine(lon1, lat1, lon2, lat2):
	"""
	    Calculate the great circle distance between two points 
	    on the earth (specified in decimal degrees)
	"""
	# convert decimal degrees to radians 
	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
	# haversine formula 
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = sin(dlat/2.)**2 + cos(lat1) * cos(lat2) * sin(dlon/2.)**2
	c = 2. * asin(sqrt(a)) 
	km = 6367. * c
	return km 

def quadraticDistance(lon1, lat1, lon2, lat2):
	"""
	    Calculate the flat quadratic diatance in degrees. Its a rough approximation. But the maps don't include the poles, so it's okay. 
	"""
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	return sqrt(dlon*dlon + dlat*dlat)

def getORCAdepth(z,depth,debug=True):
	d = 1000.
	best = -1
	for i,zz in enumerate(depth.squeeze()):
		d2 = abs(abs(z)-abs(zz))
		if d2<d:
		   d=d2
		   best = i
		   print 'getORCAdepth:',i,z,zz,depth.shape, 'best:',best
	if debug: print 'depth: in situ:', z,'index:', best, 'distance:',d,', closest model:',depth.shape, depth[best]
	return best
		
def makeLonSafe(lon):
	while True:
		if -180<lon<=180:return lon
		if lon<=-180:lon+=360.
		if lon> 180:lon-=360.	
			
def makeLonSafeArr(lon):
	if lon.ndim == 2:
	 for l,lon1 in enumerate(lon):
	  for ll,lon2 in enumerate(lon1):
	   lon[l,ll] = makeLonSafe(lon2)
	 return lon
	if lon.ndim == 1:
	 for l,lon1 in enumerate(lon):
	   lon[l] = makeLonSafe(lon1)
	 return lon
	 	 
	assert False		
	
	  

def getMonthFromSecs(secs):
	if int(secs/30.4) in xrange(0,12):return int(secs/30.4)
	day = int(secs/(24*60*60.))
	#very approximate
	return int(day/30.4) #returns month between 0 and 11 ...
	
	

	




def main():
	assert False
	



	

	pCO2 =True
	if pCO2:
		datafile = "/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/./outNetCDF/MaredatMatched-xhonp-clim/pCO2/takahashi2009_month_flux_pCO2_2006c_noHead.nc"
		#datafile = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/takahashi2009_month_flux_pCO2_2006c_noHead.nc"
			#"outNetCDF/MaredatMatched/pCO2/takahashi2009_month_flux_pCO2_2006c_noHead.nc"
		#datafile = "outNetCDF/MaredatMatched/pCO2/takahashi2009_annual_flux_pCO2_2006c_noHead.nc"
		if jobID.upper() == 'MEDUSA':		
			pco2vars = ['xxxx']
		else:		
			pco2vars = ["chl","fAirSeaC","pCO2w","netPP",]
		b = matchDataAndModel(datafile, ModeldiagFile,pco2vars,key=key)	

			
			

	seawifs = 0#False
	if seawifs:
		datafile = "/data/euryale7/scratch/ledm/seawifs_monthly/SeaWiFs_climatology_1997-2007_tiny.nc"
		#datafile = "outNetCDF/MaredatMatched/pCO2/takahashi2009_annual_flux_pCO2_2006c_noHead.nc"

		Modelvars = ['P1c','Chl1','P2c','Chl2','P3c','Chl3','P4c','Chl4',] #'N1p','N3n','N4n','N5s','N7f',
				
		b = matchDataAndModel(datafile, ModelBGCFile,Modelvars,key=key)		# for PFT Chl	

		Modelvars = ["chl",]
		b = matchDataAndModel(datafile, ModeldiagFile,Modelvars,key=key) 	# for total Chl



	PP = 0#True
	if PP:
		datafile = MareDatFold+"PP100108.nc"
		#modelFile= "outNetCDF/Climatologies/"+jobID+"_clim_PP.nc"
		pco2vars = ["PP",]
		b = matchDataAndModel(datafile, ModeldiagFile,pco2vars,key=key)		
		
	
	


		
	intPP =True
	if intPP:
		datafile = "/data/euryale7/scratch/ledm/LestersReportData/PPint_1deg.nc"
		#modelFile= "outNetCDF/Climatologies/"+jobID+"_clim_IntPP.nc"
		pco2vars = ["netPP","IntPP"]
		b = matchDataAndModel(datafile, ModeldiagFile,pco2vars,key=key)	
	




		

	



	

				

			



if __name__=="__main__":
	main() 
	print 'The end.'
