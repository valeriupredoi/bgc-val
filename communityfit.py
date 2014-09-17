#!/usr/bin/ipython 
#Standard Python modules:
from netCDF4 import Dataset
from sys import argv
from matplotlib import pyplot, colors
import numpy as np
from scipy import optimize
#Specific local code:
from pftnames import pftnames as pn
from UKESMpython import folder,getFileList,sliceA
from analysis import analysis
#import C2Chl
import communitystructure as cs

"""	THIS REQUIRES LEE'S Community structure MODULE:
	https://gitlab.ecosystem-modelling.pml.ac.uk/ledm/communitystructure

"""

def fitAndplot(chl, carbon, dx):
	if chl.min() < 0. or carbon.min()<0.:
		print "is this log scale or something? carbon:", carbon.min(), '\tchl:',chl.min()
		assert False
		
	communityfit_func = lambda chl,m,p: 10**m*chl**(p-1)
	popt, pcov = optimize.curve_fit(communityfit_func, chl,carbon)
	print 'FitAndPlot:',popt,pcov
	return communityfit_func(np.ma.array(dx), popt[0],popt[1])
	
	
class communityfit(analysis):
  def __init__(self,fileIn):
	analysis.__init__(self, fileIn)
	print self.fileIn, self.jobID
	
	#if self.model == 'ERSEM':
	pfts = ['Total', 'Diatoms', 'NonDiatoms' ]
	
	#for self.pft in pfts:
	self.loadcommunityfitData()
	
	for r in ['SurfaceNoArtics',]:#'Surface', 'Top40m','Top200m', 'Global']:
		self.makePlot(region = r)


  def loadcommunityfitData(self,):
  	"""	While this method can be memory intensive, it may be faster to do it this way than to load them from file each time.
  	"""
  	print "communityfit:\tINFO:\tloadData"
  	if not self.ncOpened: self.loadData()

	#c and chl have different names in each model
	if self.model=='ERSEM':
  		print "communityfit:\tINFO:\tloading ERSEM data"
		#fnD = self.fileIn.replace('P.nc','D.nc')		
		#ncD = self.loadFile(fnD)  		
		#self.chl = ncD.variables['chl'][:]
		#if self.pft == 'Total':
		#	self.c   = self.nc.variables['P1c'][:]+self.nc.variables['P2c'][:]+self.nc.variables['P3c'][:]+self.nc.variables['P4c'][:]
		
		#if self.pft == 'Diatoms':
		self.dia   = self.nc.variables['P1c'][:]
		#if self.pft == 'Flagellates':
		self.flag   = self.nc.variables['P2c'][:]
		#if self.pft == 'Picophytoplankton':
		self.pico   = self.nc.variables['P3c'][:]		
		#if self.pft == 'Large Phytoplankton':
		self.large   = self.nc.variables['P4c'][:]
			

		
	if self.model=='MEDUSA':
  		print "communityfit:\tINFO:\tloading MEDUSA data"		
		#self.chl = self.nc.variables['CHL'][:]
		#if self.pft == 'Total':		
		#	self.c   = self.nc.variables['PHD'][:]+self.nc.variables['PHN'][:]
		#if self.pft == 'Diatoms':		
		self.dia   = self.nc.variables['PHD'][0,0]*79.573
		self.nond   = self.nc.variables['PHN'][0,0]*79.573
		

	
  def makePlot(self,region='Surface',doLog=True):
  	
  	self.autoFilename('communityfit','communityfit_'+region)
  	
  	if self.existsFilename():
  		print "communityfit:\tINFO\tPlot already exists:", self.filename  
  		return
  	print "communityfit:\tINFO\tMaking plot:", self.filename 	
  	  		


	
	if self.model == 'MEDUSA':
	  	print "communityfit:\tINFO\tMaking slices:", region
	  	dia   = sliceA(self.dia, region)
	  	nond  = sliceA(self.nond,region)
		m = dia.mask + nond.mask
	  	dia  = np.ma.masked_where(m, dia ).compressed()
	  	nond = np.ma.masked_where(m, nond).compressed()
  		totalc = dia+nond		
	  	carbonpc = [100.*dia/totalc,100.*nond/totalc]
		#totalc = 
	  	titles = ['Diatoms', 'Non Diatoms']
		subplots = [211,212]
		

  	
	x_range = [-3.,1.]
	
  	fig = pyplot.figure()
  	for i,sp in enumerate(subplots):
	  	ax = fig.add_subplot(sp)
	  	
		if i ==0: 	
		  	pyplot.title(self.model +' '+self.date+' '+region+' Community Structure')	  	
	  	if sp == subplots[-1]:	
	  		pyplot.xlabel(r'log$_{10}$ Phytoplankton Carbon')
	  		
		pyplot.ylabel(titles[i]+' %')
		pyplot.xlim(x_range)		
		pyplot.ylim([-4., 102.])
		
		pyplot.hist2d(np.ma.log10(totalc),carbonpc[i],bins=(1000,52),cmin=1, norm=colors.LogNorm(),cmap=pyplot.get_cmap('Greys'))	
		
		csf = cs.comstrucFit(totalc,carbonpc[i], pft='pico', chlRange=[10.**x_range[0],10.**x_range[1]], fitTypes =['NoFit',])
		cx,cy = csf.doRunningMean()		
		#cbar=pyplot.colorbar()
		
				
	  	


  	print "analysis:\tINFO\tSaving: " + self.filename
	pyplot.savefig(self.filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	






def main():
	try:
		filesIn = getFileList(argv[1:])
		if len(filesIn)==0:assert False
		print "Using command line arguments:", filesIn	
		
	except:
		print "communityfit:\tERROR:\tNo file specified"
		#return
		filesIn = ['/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc',]
		
	
	for fn in filesIn:
		a = communityfit(fn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
	
	
