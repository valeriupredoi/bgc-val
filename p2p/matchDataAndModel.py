#!/usr/bin/ipython
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
from pftnames import getmt


#####	
# These are availalble in the module:
#	https://gitlab.ecosystem-modelling.pml.ac.uk/ledm/netcdf_manip
from pruneNC import pruneNC
from convertToOneDNC import convertToOneDNC
from mergeNC import mergeNC
from changeNC import changeNC, AutoVivification


#TO DO
#	Remove ncdfView requirement
#	This still requires the netcdf_manip library, the ncdfView code, the ORCA1bathy file

class matchDataAndModel:
  """	matchDataAndModel: 
  	This code takes the real data from in situ  measurments in netcdf format and the Model data and created two 1D matched netcdfs. 
	The 1D matched netcdfs are then used to make plots and perform statistical analysis (not in this code).
	The first step is to produce lightweight "pruned" versions of the files, which have the unused fields stripped out.
	Some of the datasets are too large to run this code on desktop machine, so in those cases we request a specific region, ie "Surface".
	Debug: prints more statements.
  """


  def __init__(self,DataFile,ModelFile,dataType, workingDir = '',DataVars='',ModelVars='',  model = '',jobID='', year='clim',region='',debug = True,):

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
	self.region = region
	self.debug = debug	
		
	self.dataType = dataType
	self._meshLoaded_ = False
	
	if debug: print  "matchDataAndModel:\tINFO:\t",self.dataType, '\tModelfile:', self.ModelFile
		
	self.compType= 'MaredatMatched-'+self.model+'-'+self.jobID+'-'+self.year
		
	if workingDir =='':
		self.workingDir = ukp.folder('/data/euryale7/scratch/ledm/ukesm_postProcessed/ukesm/outNetCDF/'+'/'.join([self.compType,self.dataType+self.region]) )
	else: 	self.workingDir = workingDir	

	self.matchedShelve 	= ukp.folder(self.workingDir)+self.model+'-'+self.jobID+'_'+self.year+'_'+'_'+self.dataType+'_'+self.region+'.shelve'
	self.matchesShelve 	= ukp.folder(['shelves','MaredatModelMatch',])+'WOAtoORCA1.shelve'

	self.workingDirTmp = 	ukp.folder(self.workingDir+'tmp')
	self.DataFilePruned=	self.workingDirTmp+'Data_' +self.dataType+'_'+self.region+'_'+self.model+'-'+self.jobID+'-'+self.year+'_pruned.nc'
	self.ModelFilePruned=	self.workingDirTmp+'Model_'+self.dataType+'_'+self.region+'_'+self.model+'-'+self.jobID+'-'+self.year+'_pruned.nc'	
	
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
		
	#WOADatas = [a+self.region for a in ['salinity','temperature','temp','sal','nitrate','phosphate','silicate',]]	   	 
	
	if self.region in ['Surface','200m','100m','500m','1000m','Transect',]:	
	    if nc.variables[self.DataVars[0]].shape in [(12, 14, 180, 360), (12, 24, 180, 360)]: # WOA format
		mmask = np.ones(nc.variables[self.DataVars[0]].shape)
		#####
		# This could be rewritten to just figure out which level is closest, but a bit much effort for not a lot of gain.
		
		if self.region in ['Surface','200m','100m','500m','1000m',]: 
			if self.region in ['Surface',]:	k = 0
			if self.region == '100m': 	k = 6			
			if self.region == '200m': 	k = 9
			if self.region == '500m': 	k = 13
			if self.region == '1000m': 	k = 18					
			mmask[:,k,:,:] = 0
						
		if self.region == 'Transect':	mmask[:,:,:,200] = 0   # Pacific Transect.	
		mmask +=nc.variables[self.DataVars[0]][:].mask
		print mmask.shape

		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking WOA style flat array:',self.DataFilePruned,'-->',self.DataFile1D	
	  	convertToOneDNC(self.DataFilePruned,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)		  	
	    else:
		print 'matchDataAndModel:\tconvertDataTo1D:\tYou need to add more file spcific regions here.', nc.variables[self.DataVars[0]].shape
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
		s.close()		
		print "matchModelToData:\tOpened shelve:", self.matchedShelve
		print "matchModelToData:\tStarting from maxindex:",maxIndex," and ",len(self.matches), " already matched. Mask:",self.maremask.sum()
	except:
		self.matches = {}
		maxIndex = 0
		self.maremask = np.zeros(is_i.shape) # zero: 

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
	for key in mt.keys():
		key = key.upper()
		if key in ['ERSEM','NEMO','MEDUSA']:continue
		try:
			if self.dataType in mt[key].keys() and key not in ytype:
				ytype.append(key)
		except:pass
	if len(ytype) == 1:
		ytype = ytype[0]		
	else:
		print "matchModelToData:\tUnable to determine data dataset type (ie, Maredat, WOA, etc...)", ytype

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
	tdict = {i+1:i for i in xrange(12)}
	ncIS.close()	     

	#####
	# This list could be removed by adding a check dimensionality of data after it was been pruned.	    
	if self.dataType in ['pCO2','seawifs','Seawifs','mld_DT02', 'mld_DR003','mld_DReqDTm02','mld',]: 
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
			z = getORCAdepth(wz,self.deptht,debug=True)
			zdict[wz]	= z
			if self.debug: print "matchModelToData:\t",i, 'Found new depth:',wz,'-->',z								

		#####			
		#Match Time	
		try:
			t = tdict[wt]
		except:
			t = getMonthFromSecs(wt)
			tdict[wt] = t
			if self.debug:	print "matchModelToData:\t",i, 'Found new month:', wt, '-->',t

		#####		
		# Add match into array
		try:
			tmp = self.matches[(t,z,la,lo)][0]
			self.maremask[i] = tmp
			self.matches[(t,z,la,lo)].append(i)
			print "matchModelToData:\tWARNING:",i,[wt,wz,wla,wlo], '-->',(t,z,la,lo),'already matched', tmp
		except:
			self.matches[(t,z,la,lo)]=[i,]
			self.maremask[i] = i

		#####			
		#increment by 1 to save/ end, as it has finished, but i is used as a starting point.
		i+=1
		if i%100==0:	
			if self.debug: print "matchModelToData:\t", i,ii,self.dataType,self.region,':\t',[wt,wz,wla,wlo] ,'--->',[t,z,la,lo]
			
		if i>1 and i%500000==0:
			if self.debug: print "matchModelToData:\tSaving Shelve: ", self.matchedShelve
			s = shOpen(self.matchedShelve)
			s['matches']  = self.matches
			s['maxIndex'] = i

			s['maremask'] = self.maremask			
			s.close()
			s = shOpen(self.matchesShelve)		
			s['lldict'] = lldict 
			s.close()				
			
	print "matchDataAndModel:\tSaving Shelve", self.matchedShelve	
	s = shOpen(self.matchedShelve)
	s['matches']  = self.matches
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
		
	maremask = np.array([float(a) for a in self.maremask] )
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
		
		for m in maremaskDict.keys():
		    values = []
		    for i in maremaskDict[m]:
		    	values.append(arr[i])
		    	
		    out.append(np.median(values))	
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
	ncIS.close()
	
	if self.debug: 
		print "matchDataAndModel:\tapplyMaskToData:\tNEW MASK for 1D Maredat:",maremask.shape, maremask.sum(), maremask.size, maremask.size-maremask.sum()
		print "matchDataAndModel:\tapplyMaskToData:\tMaking 1D Maredat file:", self.DataFile1D,'-->', self.maskedData1D
		
	c = changeNC(self.DataFile1D, self.maskedData1D,av,debug = self.debug)
	

  def loadMesh(self,):
      	# This won't work unless its the 1 degree grid.
  	print "matchModelToData:\tOpened Model netcdf: ~/data/mesh_mask_ORCA1_75.nc"
  	ncER = Dataset("data/mesh_mask_ORCA1_75.nc",'r')
  	#ncER = ncdfView("data/mesh_mask_ORCA1_75.nc",Quiet=True)  	
	self.latcc    = ncER.variables['nav_lat'][:]
	self.loncc    = makeLonSafeArr(ncER.variables['nav_lon'][:])
  	self.deptht   = ncER.variables['gdept_0'][:]
	ncER.close()
	self._meshLoaded_ = 1
	
  def getOrcaIndexCC(self,lat,lon,debug=True,slowMethod=False,llrange=5.):
	""" takes a lat and long coordinate, an returns the position of the closest coordinate in the NemoERSEM (ORCA1) grid.
	    uses the bathymetry file.
	"""
	km = 10.E20
	la_ind, lo_ind = -1,-1
	lat = ukp.makeLatSafe(lat)
	lon = ukp.makeLonSafe(lon)	
	
	if not self._meshLoaded_:self.loadMesh()
	c = (self.latcc - lat)**2 + (self.loncc - lon)**2

	(la_ind,lo_ind) =  np.unravel_index(c.argmin(),c.shape)

	#km2 = abs(haversine((lon, lat), (locc,lacc)))
	if debug: print 'location ', [la_ind,lo_ind],'(',self.latcc[la_ind,lo_ind],self.loncc[la_ind,lo_ind],') is closest to:',[lat,lon]	
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
	#km = 6367. * c
	return c 

def quadraticDistance(lon1, lat1, lon2, lat2):
	"""
	    Calculate the flat quadratic diatance in degrees. Its a rough approximation. But the maps don't include the poles, so it's okay. 
	"""
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	return sqrt(dlon*dlon + dlat*dlat)

def getORCAdepth(z,depth,debug=True):
	d = 10000.
	best = -1
	for i,zz in enumerate(depth):
		d2 = abs(abs(z)-abs(zz))
		if d2<d:
		   d=d2
		   best = i
	if debug: print 'depth: in situ:', z,', closest model:',depth[best], 'index:', best, 'distance:',d
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
