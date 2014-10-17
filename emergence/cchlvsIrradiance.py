#!/usr/bin/ipython 
#Standard Python modules:
from netCDF4 import Dataset
from sys import argv
from matplotlib import pyplot
import numpy as np
#Specific local code:
from pftnames import pftnames as pn
from UKESMpython import folder,getFileList,sliceA
from analysis import analysis


#	I'm running this as a attempt before trying to do something that counts.



class cchlvsIrradiance(analysis):
  def __init__(self,fileIn):
	analysis.__init__(self, fileIn)
	print self.fileIn, self.jobID
	self.loadcchlData()
	
	for r in ['Surface', 'SurfaceNoArtics', 'Global']:
		self.makePlot(region = r)


  def loadcchlData(self,):
  	"""	While this method can be memory intensive, it may be faster to do it this way than to load them from file each time.
  	"""
  	print "cchlvsIrradiance:\tINFO:\tloadData"
  	if not self.ncOpened: self.loadData()
	self.ncP = self.nc
	if self.model=='MEDUSA':
		self.autoDate()
		#fnD = '/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/xhont-'+self.date+'/'
		fnD = '/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/xhonp-'+self.date+'/xhonp_'+self.date+'_Diag.nc'
		# loading clim file for the same period 
		
	if self.model=='ERSEM':

		fnD = self.fileIn.replace('P.nc','D.nc')		
		
	ncD = self.loadFile(fnD)
	
	self.eir = ncD.variables['EIR'][:]
	
	if self.model=='ERSEM':
  		print "cchlvsIrradiance:\tINFO:\tloading ERSEM data"
		self.chl = ncD.variables['chl'][:]
		self.c   = self.nc.variables['P1c'][:]+self.nc.variables['P2c'][:]+self.nc.variables['P3c'][:]+self.nc.variables['P4c'][:]
	if self.model=='MEDUSA':
  		print "cchlvsIrradiance:\tINFO:\tloading MEDUSA data"		
		self.chl = self.nc.variables['CHL'][:]
		self.c   = (self.nc.variables['PHD'][:]+self.nc.variables['PHN'][:])*79.573

	
  def makePlot(self,region='Surface'):
  	
  	self.autoFilename('cchlLight','cchlLight_'+region)
	#self.filename = folder('images/analysis')+'_'.join(['analysis',self.date,self.jobID])+'.png'
  	print "cchlvsIrradiance:\tINFO\tMaking plot:", self.filename 	
  	fig = pyplot.figure()
  	ax = fig.add_subplot(111)
  	
  	#pyplot.pcolormesh

  	#etc
  	
  	cchl = sliceA(self.c,region)/sliceA(self.chl,region)
  	light = sliceA(self.eir,region)

	cchl  = np.ma.masked_where((cchl>10E10) +(cchl<10E-05),cchl)
	light = np.ma.masked_where(light>10E10,light)	
	
	cchl = 	np.ma.masked_where(cchl.mask + light.mask,cchl)
	light = np.ma.masked_where(cchl.mask + light.mask,light)
			
  	h = pyplot.hexbin(cchl.compressed(), light.compressed(),  bins='log',mincnt=1,) #xscale='log', yscale='log',
  	
  	pyplot.title(self.model +' '+self.date+' '+region+' Carbon:Chlorophyll vs Light ')
	pyplot.xlabel('Carbon:Chlorophyll ratio')
	pyplot.ylabel('Light')
		
  	print "analysis:\tINFO\tSaving: " + self.filename
	pyplot.savefig(self.filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	






def main():
	try:
		filesIn = getFileList(argv[1:])
		if len(filesI)==0:assert False
		print "Using command line arguments:", filesIn	
		
	except:
		print "cchlvsIrradiance:\tERROR:\tNo file specified"
		#return
		filesIn = ['/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc',]
		
	
	for fn in filesIn:
		a = cchlvsIrradiance(fn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
	
	
