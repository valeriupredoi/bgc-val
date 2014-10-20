#!/usr/bin/ipython
from sys import argv
from os.path import exists,split, getmtime, basename
from glob import glob
from shelve import open as shOpen
from shutil import copy2

from matplotlib.colors import LogNorm
from matplotlib import pyplot, ticker
#from numpy.ma import  min as mmin, max as mmax, masked_where, array as ma,zeros as mzeros, masked,divide,nonzero,array as marray
#from numpy import isnan as npNaN, isinf as npInf, array,arange, delete, tile,NaN,power,zeros,median
import numpy as np
from numpy import array,ones,unravel_index
from math import radians, cos, sin, asin, sqrt
from haversine import haversine
from netCDF4 import num2date
from datetime import datetime

#local imports
from ncdfView import ncdfView
from deMoraTools import getLogTicks,folder
#import iMarNetPython as impy
#from iMarNetPython import noLeapConvert, shouldIMakeFile,depthKeys,alwaysInclude,getFileList,getMareDatatype,getpco2,getHenryConstant,getORCAdepth#getJobID
#from iMarNetPython import makeLatSafe, makeLonSafe,makeLonSafeArr
from UKESMpython import shouldIMakeFile

# These are availalble in the module:
#	https://gitlab.ecosystem-modelling.pml.ac.uk/ledm/netcdf_manip
from pruneNC import pruneNC
from convertToOneDNC import convertToOneDNC
from mergeNC import mergeNC
from changeNC import changeNC, AutoVivification
#from addERSEMtoXL import getORCAdepth# getOrcaIndexCC
#from xtraPlots import hovmoeller



class matchDataAndModel:
  """	matchDataAndModel: 
  	This code takes the real data from in situ  measurments in netcdf format and the Model data and created two 1D matched netcdfs. 
	The 1D matched netcdfs are then used to make plots and perform statistical analysis.
	Some of the datasets are too large to run this code on desktop machine, so in those cases we request a specific region, ie "Surface".
	Debug = True prints more statements.
  """


  def __init__(self,DataFile,ModelFile,dataType, workingDir = '',DataVars='',  jobID='', year='clim',region='',debug = True,):

	if debug:
		print "matchDataAndModel:\tINFO:\tStarting matchDataAndModel"
		print "matchDataAndModel:\tINFO:\tData file: \t",DataFile
		print "matchDataAndModel:\tINFO:\tModel file: \t",ModelFile
		print "matchDataAndModel:\tINFO:\tdataType: \t",dataType	
	self.DataFile=DataFile
	self.ModelFile = ModelFile 	
	self.DataVars=DataVars	


	#self.ModelVars.extend( ['time_counter','nav_lat','nav_lon','deptht'])


	self.jobID = jobID 
	self.year = year
	self.region = region
	self.debug = debug	
		
	self.dataType = dataType
		#getMareDatatype(self.DataFile) 
		# this means that the maredat file input determines which data to process, which we don't want.
		# ie, what if there were two different datasets in one file.
	
	if debug: print  "matchDataAndModel:\tINFO:\t",self.dataType, '\tModelfile:', self.ModelFile
		
	ekey = self.jobID+'-'+self.year
	self.compType= 'MaredatMatched-'+ekey
		
	if workingDir =='':
		self.workingDir = folder('/data/euryale7/scratch/ledm/ukesm_postProcessed/ukesm/outNetCDF/'+'/'.join([self.compType,self.dataType+self.region]) )
	else: 	self.workingDir = workingDir	
	self.workingDirTmp = 	folder(self.workingDir+'/tmp/')
	self.matchedShelve 	= folder(['/tmp','shelves','MaredatModelMatch',self.dataType+self.region])+self.jobID+'_'+self.dataType+'.shelve'
	self.matchesShelve 	= folder(['shelves','MaredatModelMatch',])+'WOAtoORCA1.shelve'

	self.DataFile1D  	= self.workingDirTmp +basename(self.DataFile).replace('.nc','') +self.region+'.nc'	
	self.maskedData1D	= self.workingDir + basename(self.DataFile1D).replace('.nc','')+'_2.nc'
	self.Model1D     	= self.workingDir +basename(self.ModelFile).replace('.nc','')+'_'+self.dataType+self.region+'_1D.nc' 


	
	self.run()

	
  def run(self,):
	"""There are two methods written for manipulating data.
	   One is designed to work with WOA formats, the other with MAREDAT formats.
	   Other data formats are run manually.
	"""


	
	self._convertDataTo1D_()	
	self._matchModelToData_()
	self._convertModelToOneD_()
	self._applyMaskToData_()





  def _convertDataTo1D_(self,):
   	""" This routine reduces the In Situ data into a 1D array of data with its lat,lon,depth and time components.
  	"""
		
	if not shouldIMakeFile(self.DataFile,self.DataFile1D,debug=False):
		print "matchDataAndModel:\tconvertDataTo1D:\talready exists:",self.DataFile1D
		return

	if self.dataType in ['pCO2','iron',]: 
		 # pCO2, iron Files area already 1D
		self.DataFile1D = self.DataFile 
		print "matchDataAndModel:\tconvertDataTo1D:\tpCO2/Iron File has already been converted to 1D", 'Making',self.DataFile
		return
		
  	
	if self.dataType in ['temp','sal','nit','nitrate','phosphate','silicate']:	#World Ocean Atlas format
		nc = ncdfView(self.DataFile,Quiet=True)
		mmask = ones(nc(self.DataVars[0]).shape)
		
		if self.region in ['Surface','200m','100m','500m','1000m',]:
			if self.region in ['Surface',]:
				k = 0
		        if self.dataType in ['silicate','nitrate','phosphate','sal','temp',]:			
				if self.region == '100m': k = 6			
				if self.region == '200m': k = 9
				if self.region == '500m': k = 13
				if self.region == '1000m': k = 18					
			mmask[:,k,:,:] = 0
						
		if self.region == 'Transect':	mmask[:,:,:,200] = 0   # Pacific Transect.	
		if self.region in ['All','']:	mmask[:,] = 0   # Entire Dataset
		mmask +=nc(self.DataVars[0]).mask
		print mmask.shape
		nc.close()
		
	  	convertToOneDNC(self.DataFile,self.DataFile1D ,newMask=mmask, variables = self.DataVars, debug=True)	
	  	
	else:
		print 'matchDataAndModel:\tconvertDataTo1D:\tMaking',self.DataFile,'-->',self.DataFile1D
	  	if len(self.DataVars):	convertToOneDNC(self.DataFile, self.DataFile1D, debug=True, variables = self.DataVars)
	  	else:			convertToOneDNC(self.DataFile, self.DataFile1D, debug=True)
  		
	
	
		
  	
  def _matchModelToData_(self,):
  	print "matchModelToData:\tOpened MAREDAT netcdf:", self.DataFile1D
  	
  	ncIS = ncdfView(self.DataFile1D,Quiet=True)
	is_i	= ncIS('index')[:]
	
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
	if maxIndex+1 <len(is_i):
	    zdict={}
	    tdict={}
	    	
	    if self.dataType in ['temp','sal','nitrate','phosphate','silicate',]:
  	   	 is_t	= ncIS('time')[:]	
  	   	 is_z 	= ncIS('depth')[:]
  	   	 is_la	= ncIS('lat')[:]
	   	 is_lo 	= ncIS('lon')[:]	    
	    elif self.dataType in ['pCO2',]:
  	   	 is_t	= ncIS('MONTH')[:]	
  	   	 is_z 	= ncIS('index_z')[:]
  	   	 is_la	= ncIS('LAT')[:]
	   	 is_lo 	= ncIS('LON')[:]
	   	 tdict = {i+1:i for i in xrange(12)}	
	   	 zdict = {0:0, 0.:0}  
	    elif self.dataType in ['nit',]:
  	   	 is_t	= ncIS('index_t')[:]
  	   	 is_z 	= ncIS('index_z')[:]
  	   	 is_la	= ncIS('latitude')[:]
	   	 is_lo 	= ncIS('longitude')[:]
	   	 tdict = {i:i for i in xrange(12)}
	    elif self.dataType in ['seawifs','Seawifs',]: #'chl','Tchl','Micro_percent_Tchl','Nano_percent_Tchl','Pico_percent_Tchl',
  	   	 is_t	= ncIS('month')[:]
  	   	 is_z 	= ncIS('deptht')[:] 
  	   	 is_la	= ncIS('latitude')[:]
	   	 is_lo 	= ncIS('longitude')[:]
	   	 tdict = {i+1:i for i in xrange(12)}
	   	 zdict = {0:0, 0.:0}  
	    elif self.dataType in ['intPP',]:
  	   	 is_t	= ncIS('index_t')[:]
  	   	 is_z 	= ncIS('index_z')[:]
  	   	 is_la	= ncIS('LATITUDE')[:]
	   	 is_lo 	= ncIS('LONGITUDE')[:]
	   	 tdict = {i:i for i in xrange(12)}
	    elif self.dataType in ['iron',]:
  	   	 #fulltime	= num2date(ncIS('time')[:],ncIS.variables['time'].units)
  	   	 #is_t  = array([t.month for t in fulltime])
  	   	 is_t	= ncIS('MONTH')[:]
  	   	 is_z 	= ncIS('DEPTH')[:]
  	   	 is_la	= ncIS('Latitude')[:]
	   	 is_lo 	= ncIS('Longitude')[:]
	   	 tdict = {i+1:i for i in xrange(12)}
	    elif self.dataType in ['mld_DT02', 'mld_DR003','mld_DReqDTm02']:
  	   	 #fulltime	= num2date(ncIS('time')[:],ncIS.variables['time'].units)
  	   	 #is_t  = array([t.month for t in fulltime])
  	   	 is_t	= ncIS('time')[:]
  	   	 is_z 	= np.ma.zeros(len(is_t))[:]
  	   	 is_la	= ncIS('lat')[:]
	   	 is_lo 	= ncIS('lon')[:]
	   	 tdict = {i+1:i for i in xrange(12)}
	    else:	    
  	   	 #is_t	= ncIS('index_t')[:] # range 0->11	  	   	 
  	   	 is_t	= ncIS('TIME')[:] # range 1->12	
  	   	 is_z 	= ncIS('DEPTH')[:]
  	   	 is_la	= ncIS('LATITUDE')[:]
	   	 is_lo 	= ncIS('LONGITUDE')[:]
	   	 #tdict = {i:i for i in xrange(12)}
	   	 tdict = {i+1:i for i in xrange(12)}		   	 
	    ncIS.close()

  	    print "matchModelToData:\tOpened Model netcdf: /users/modellers/ledm/workspace/iMarNet/outNetCDF/mesh_mask_ORCA1_75.nc"
  	    ncER = ncdfView("/users/modellers/ledm/workspace/iMarNet/outNetCDF/mesh_mask_ORCA1_75.nc",Quiet=True)
	    self.latcc    = ncER('nav_lat')[:]
	    self.loncc    = makeLonSafeArr(ncER('nav_lon')[:])
	    
  	    deptht   = ncER('gdept_0')[:]
	    ncER.close()


		
	    for i,ii in enumerate(is_i[maxIndex:]):
		i+=maxIndex
		wt  = is_t[i]
		wz  = is_z[i]
		wla = is_la[i]
		wlo = is_lo[i] 
		
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

		#Match Depth	
		try:
			z = zdict[wz]

		except:	 
			z = getORCAdepth(wz,deptht,debug=True)
			zdict[wz]	= z
			if self.debug: print "matchModelToData:\t",i, 'Found new depth:',wz,'-->',z								
			
		#Match Time	
		try:
			t = tdict[wt]
		except:
			t = getMonthFromSecs(wt)
			tdict[wt] = t
			if self.debug: print "matchModelToData:\t",i, 'Found new month:', wt, '-->',t
		

		# Add match into array
		try:
			tmp = self.matches[(t,z,la,lo)][0]
			self.maremask[i] = tmp
			self.matches[(t,z,la,lo)].append(i)
			print "matchModelToData:\tWARNING:",i,[wt,wz,wla,wlo], '-->',(t,z,la,lo),'already matched', tmp
		except:
			self.matches[(t,z,la,lo)]=[i,]
			self.maremask[i] = i
			
		#increment by 1 to save/ end, as it has finished, but i is used as a starting point.
		i+=1
		if i%100==0:	
			if self.debug: print "matchModelToData:\t", i,ii,self.dataType,':\t',[wt,wz,wla,wlo] ,'--->',[t,z,la,lo]
			
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
				    
	try:ncIS.close()
	except:pass
	print "matchModelToData:\tFinsished with ",maxIndex+1,"\tfinished:",len(self.matches) #, "Mask:",self.maremask.sum()
	
	
  def _convertModelToOneD_(self,):
	if not shouldIMakeFile(self.ModelFile,self.Model1D,debug=True):
		print "convertModelToOneD:\tconvertModelToOneD:\talready exists:",self.Model1D
		return	
	
	print "convertModelToOneD:\tconvertModelToOneD:\tMaking 1D Model file:", self.ModelFile,'-->', self.Model1D
  	convertToOneDNC(self.ModelFile, self.Model1D,newMask='',debug=self.debug,dictToKeep=self.matches)



	



  def _applyMaskToData_(self,):
  	""" This routine applies the mask of the model to the data. 
  	    It is needed because there are some points where two WOA grid cells fall into the same ORCA1 grid, these excess points should be meaned/medianed.
  	    Similarly, some data points fall into a masked grid cell in the model and need to be masked in the data.
  	"""
	
	if not shouldIMakeFile(self.ModelFile,self.maskedData1D,debug=True):
		print "applyMaskToData:\tapplyMaskToData:\t", "already exists:",self.maskedData1D
		return	
		
	maremask = array([float(a) for a in self.maremask] )
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
		arr= array(out).squeeze()
		print arr.shape
		return arr
		
	def applyMask(arr):
		print "applyMask: 1d:",arr.shape,maremask.shape,'-->',
		out = []
		for i,m in enumerate(maremask):
			if not m: continue
			out.append(arr[i])	
		arr= array(out).squeeze()
		print arr.shape
		return arr


  	ncIS = ncdfView(self.DataFile1D,Quiet=True)		
	av = AutoVivification()
	for v in ncIS.variables.keys():
		if ncIS(v).ndim != 1: 
			if self.debug: print "matchDataAndModel:\tapplyMaskToData:\tERROR:\tthis is suppoed to be the one D file"
			assert False

		if len(ncIS(v)[:])== len(self.maremask):
			if self.debug: print  "matchDataAndModel:\tapplyMaskToData:AutoViv:", v ,'is getting a mask.'
			av[v]['convert'] = getMedianVal
	ncIS.close()
	
	if self.debug: 
		print "matchDataAndModel:\tapplyMaskToData:\tNEW MASK for 1D Maredat:",maremask.shape, maremask.sum(), maremask.size, maremask.size-maremask.sum()
		print "matchDataAndModel:\tapplyMaskToData:\tMaking 1D Maredat file:", self.DataFile1D,'-->', self.maskedData1D
		
	c = changeNC(self.DataFile1D, self.maskedData1D,av,debug = self.debug)
	

  def getOrcaIndexCC(self,lat,lon,debug=True,slowMethod=False,llrange=5.):
	""" takes a lat and long coordinate, an returns the position of the closest coordinate in the NemoERSEM (ORCA1) grid.
	    uses the bathymetry file.
	"""
	km = 10.E20
	la_ind, lo_ind = -1,-1
	lat = makeLatSafe(lat)
	lon = makeLonSafe(lon)	
	
	c = (self.latcc - lat)**2 + (self.loncc - lon)**2

	(la_ind,lo_ind) =  unravel_index(c.argmin(),c.shape)

	#km2 = abs(haversine((lon, lat), (locc,lacc)))
	if debug: print 'location ', [la_ind,lo_ind],'(',self.latcc[la_ind,lo_ind],self.loncc[la_ind,lo_ind],') is closest to:',[lat,lon]	
	return la_ind,lo_ind
		
	
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
	# don't call main. This should be moved elsewhere.
	try: 	
		jobID = argv[1]
		key   = argv[2]
		print "Using command line arguments:", argv[1], argv[2]
	except:
		jobID= 'xhonp'
		key = "clim"
		print "Not using command line arguments,Defaults:",jobID,key
		

		##jobID = "xjeza"	
		#jobID = "xhono"	

	if jobID.upper() == 'MEDUSA':
	
		ppdataFold="/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
		ModelBGCFile 	= ppdataFold+"medusa_bio_"+key+".nc"
		#SalTempWindFile	= ppdataFold+"outNetCDF/"+jobID+"-"+key+"/"+jobID+"_"+key+"_SalTempWind.nc"
		ModeldiagFile 	= ModelBGCFile#ppdataFold+"outNetCDF/"+jobID+"-"+key+"/"+jobID+"_"+key+"_Diag.nc"
				
	else:
		ppdataFold="/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/"
		ModelBGCFile 	= ppdataFold+"outNetCDF/"+jobID+"-"+key+"/"+jobID+"_"+key+"_ERSEM.nc"
		SalTempWindFile	= ppdataFold+"outNetCDF/"+jobID+"-"+key+"/"+jobID+"_"+key+"_SalTempWind.nc"
		ModeldiagFile 	= ppdataFold+"outNetCDF/"+jobID+"-"+key+"/"+jobID+"_"+key+"_Diag.nc"

	MareDatFold="/data/perseus2/scratch/ledm/MAREDAT/MAREDAT/"		
	#ERSEMBGCFile 	= "outNetCDF/Climatologies/xhono_clim.nc"
	#SalTempWindFile= "outNetCDF/Climatologies/xhono_clim_SalTempWind.nc"
	#ERSEMdiagFile 	= "outNetCDF/Climatologies/xhono_clim_Diag.nc"


	mld = 0	#True
	if mld:
	    mldvars = ['somxl010',]
	    for datafile in ["/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DT02_c1m_reg2.0.nc",
	    			"/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DR003_c1m_reg2.0.nc",
	    			"/data/euryale7/scratch/ledm/IFREMER-MLD/mld_DReqDTm02_c1m_reg2.0.nc", ]:
	    
		b = matchDataAndModel(datafile, SalTempWindFile,mldvars,key=key)
			

	

	



	
	
	#### MAREDAT files:
		
	diatoms = True
	if diatoms:
		diatomsfn 	= MareDatFold+"MarEDat20120716Diatoms.nc"
		if jobID.upper() == 'MEDUSA':				
			diatomsvars = ["PHD",]
		else:		
			diatomsvars 	= ['N1p','N3n','N4n','N5s','N7f','P1c','P1f','P1n','P1p','P1s','Chl1']
		b = matchDataAndModel(diatomsfn, ModelBGCFile,diatomsvars,key=key)	


	chl =True
	if chl:
		chlfn 	= MareDatFold+"MarEDat20121001Pigments.nc"		
		#datafile = "/data/euryale7/scratch/ledm/LestersReportData/Seawifs_chl.nc"
		if jobID.upper() == 'MEDUSA':				
			chlvars = ["CHL",]
		else:
			chlvars = ["chl",]
		
		b = matchDataAndModel(chlfn, ModeldiagFile,chlvars,key=key)	
		
	if jobID.upper() != 'MEDUSA':		
		bac = True
		if bac:
			bacfn 		= MareDatFold+"MarEDat20120214Bacteria.nc"
			bacvariables 	= ['B1c','B1n','B1p', 'N1p','N3n','N4n','N5s','N7f']		
			b = matchDataAndModel(bacfn, ModelBGCFile,bacvariables,key=key)

		picophyto= True
		if picophyto:
			picophytofn 	= MareDatFold+"MarEDat20111206Picophytoplankton.nc"
			picovars 	= ['N1p','N3n','N4n','N5s','N7f','P3c','P3f','P3n','P3p','P3s','Chl3']			
			b = matchDataAndModel(picophytofn, ModelBGCFile,picovars,key=key)
			
	microzoo= True
	if microzoo :
	
		microzoofn 	= MareDatFold+"MarEDat20120424Microzooplankton.nc"
		if jobID.upper() == 'MEDUSA':				
			microzoovars = ["ZMI",]
		else:
			microzoovars 	= ['N1p','N3n','N4n','N5s','N7f','Z4c','Z5c','Z5n', 'Z5p', 'Z6c', 'Z6n', 'Z6p',]			
		b = matchDataAndModel(microzoofn, ModelBGCFile,microzoovars,key=key)	
	
	macrozoo = False
	if macrozoo :
		macrozoofn 	= MareDatFold+"MarEDat20120216Macrozooplankton.nc"
		macrozoovars 	= ['N1p','N3n','N4n','N5s','N7f','Z4c','Z5c','Z5n', 'Z5p', 'Z6c', 'Z6n', 'Z6p',]
		b = matchDataAndModel(macrozoofn, modelFile,macrozoovars,key=key)

	mesozoo = True
	if mesozoo:
		mesozoofn 	= MareDatFold+"MarEDat20120705Mesozooplankton.nc"
		if jobID.upper() == 'MEDUSA':				
			mesozoovars = ["ZME",]
		else:		
			mesozoovars 	= ['N1p','N3n','N4n','N5s','N7f','Z4c','Z5c','Z5n', 'Z5p', 'Z6c', 'Z6n', 'Z6p',]
		b = matchDataAndModel(mesozoofn, ModelBGCFile,mesozoovars,key=key)	
			
	#return


	iron = True
	if iron:
		datafile = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
		if jobID.upper() == 'MEDUSA':		
			ironvars = ['FER',]
		else:
			ironvars = ['N1p','N3n','N4n','N5s','N7f',]
		b = matchDataAndModel(datafile, ModelBGCFile,ironvars,key=key)	



				

#	assert False
	
	

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

		#also need temperature for pCO2 calculation

		#modelFile= SalTempWindFile
		#pco2vars = temperaturevars
		#b = matchDataAndModel(datafile, modelFile,pco2vars,key=key)	
		
		#dfile  =  "outNetCDF/MaredatMatched/pCO2/xhono_clim_Diag_pCO2_1D.nc"
		#tfile  =  "outNetCDF/MaredatMatched/pCO2/xhono_clim_SalTempWind_pCO2_1D.nc"
		#outfile=  "outNetCDF/MaredatMatched/pCO2/xhono_clim_pCO2_1D.nc"
		#r = recalculate_pCO2(dfile,tfile,outfile)
			
			
			

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
		

	#pCO2 = True
	#if pCO2:
#		datafile = "outNetCDF/MaredatMatched-"+jobID+"-"+key+"/pCO2/takahashi2009_month_flux_pCO2_2006c_noHead.nc"
		#try:
#		copy2("outNetCDF/MaredatMatched-xhont-clim/pCO2/takahashi2009_month_flux_pCO2_2006c_noHead.nc", datafile)
	#	except: pass
		#datafile = "outNetCDF/MaredatMatched/pCO2/takahashi2009_annual_flux_pCO2_2006c_noHead.nc"

#		pco2vars = ["chl","fAirSeaC","pCO2w","netPP",]
#		b = matchDataAndModel(datafile, ModeldiagFile,pco2vars,key=key)



		

	for woaKey in ['silicate','nitrate','phosphate','salinity','temperature',]:
		if woaKey in ['salinity','temperature',]:
			continue
			modelFile = SalTempWindFile
		else:
			modelFile = ModelBGCFile
		datafile = '/data/euryale7/scratch/ledm/WOA/'+woaKey+'_monthly_1deg.nc'
		
		l = woaKey[0]
		if woaKey == 'silicate': l='i' # silicate key is i, not s.
		m_temp_vars = [l+'_mn','depth','lat','lon','time'] 
		temperaturevars = ['votemper',] 
		for s in ['Surface','500m','100m','200m',]:#'Transect',]:#'All',
		    
		    b = matchDataAndModel(datafile, modelFile,temperaturevars,DataVars = m_temp_vars, region = s,key=key)

		for s in ['1000m',]:#'Transect','Surface',]:#'All',
		    if  woaKey in ['silicate','nitrate','phosphate',]:#'salinity','temperature',]: 
		    	continue
		    
		    b = matchDataAndModel(datafile, modelFile,temperaturevars,DataVars = m_temp_vars, region = s,key=key)
	#assert False

	


	
	
	
	
	chl =0#True
	if chl:
		print "DON'T USE THIS ONE, THERE IS A "
		datafile = "/data/euryale7/scratch/ledm/LestersReportData/Seawifs_chl.nc"

		pco2vars = ["chl",]
		b = matchDataAndModel(datafile, ModeldiagFile,pco2vars,key=key)	

	
	GLODAPpCO2 =0#True
	if GLODAPpCO2:
		datafile = "/data/euryale7/scratch/ledm/LestersReportData/GLODAP_TCO2.nc"

		modelFile= ModeldiagFile
		pco2vars = ["chl","fAirSeaC","pCO2w","netPP",]
		b = matchDataAndModel(datafile, modelFile,pco2vars,key=key)	



	
	intPP =True
	if intPP:
		datafile = "/data/euryale7/scratch/ledm/LestersReportData/PPint_1deg.nc"
		#modelFile= "outNetCDF/Climatologies/"+jobID+"_clim_IntPP.nc"
		pco2vars = ["netPP","IntPP"]
		b = matchDataAndModel(datafile, ModeldiagFile,pco2vars,key=key)	
	




		

	



	

				

			



if __name__=="__main__":
	main() 
	print 'The end.'
