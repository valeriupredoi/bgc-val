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
# ledm@pml.ac.uk
#

#Standard Python modules:
from netCDF4 import Dataset
from sys import argv
from matplotlib import pyplot
import numpy as np
from scipy import optimize
#Specific local code:
#from pftnames import pftnames as pn
from UKESMpython import folder,getFileList,sliceA
from analysis import analysis
from bgcvaltools import C2Chl


"""	THIS REQUIRES MOMME'S C2CHL MODULE:
	https://gitlab.ecosystem-modelling.pml.ac.uk/momm/pml-python-tools/blob/master/C2Chl.py
"""

def fitAndplot(chl, carbon, dx):
	if chl.min() < 0. or carbon.min()<0.:
		print "is this log scale or something? carbon:", carbon.min(), '\tchl:',chl.min()
		assert False
		
	cchl_func = lambda chl,m,p: 10**m*chl**(p-1)
	popt, pcov = optimize.curve_fit(cchl_func, chl,carbon)
	print 'FitAndPlot:',popt,pcov
	return cchl_func(np.ma.array(dx), popt[0],popt[1])
	
	
class cchl(analysis):
  def __init__(self,fileIn):
	analysis.__init__(self, fileIn)
	self.dataloaded = False
	print self.fileIn, self.jobID
	
	#if self.model == 'ERSEM':
	pfts = ['Total', 'Diatoms', 'NonDiatoms' ]
	
	for self.pft in pfts:
	
		for r in ['Surface', 'SurfaceNoArtics','Top40m','Top200m', 'Global']:
			self.makePlot(region = r)


  def loadcchlData(self,):
  	"""	While this method can be memory intensive, it may be faster to do it this way than to load them from file each time.
  	"""
  	print "cchl:\tINFO:\tloadData"
  	if not self.ncOpened: self.loadData()
	#self.ncP = self.nc
	#if self.model=='MEDUSA':
	#	self.autoDate()
	#	fnD = '/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/xhonp-'+self.date+'/xhonp_'+self.date+'_Diag.nc'
	#	ncD = self.loadFile(fnD)	
	#if self.model=='ERSEM':
	#	fnD = self.fileIn.replace('P.nc','D.nc')		
	#	ncD = self.loadFile(fnD)
	
	# Same in both:
	#self.eir = ncD.variables['EIR'][:]
	
	#c and chl have different names in each model
	if self.model=='ERSEM':
  		print "cchl:\tINFO:\tloading ERSEM data"
		fnD = self.fileIn.replace('P.nc','D.nc')		
		ncD = self.loadFile(fnD)  		
		self.chl = ncD.variables['chl'][:]
		if self.pft == 'Total':
			self.c   = self.nc.variables['P1c'][:]+self.nc.variables['P2c'][:]+self.nc.variables['P3c'][:]+self.nc.variables['P4c'][:]
		if self.pft == 'Diatoms':
			self.c   = self.nc.variables['P1c'][:]
		if self.pft == 'NonDiatoms':
			self.c   = self.nc.variables['P2c'][:]+self.nc.variables['P3c'][:]+self.nc.variables['P4c'][:]
		#self.poc   = 
	if self.model=='MEDUSA':
  		print "cchl:\tINFO:\tloading MEDUSA data"		
		self.chl = self.nc.variables['CHL'][:]
		if self.pft == 'Total':		
			self.c   = self.nc.variables['PHD'][:]+self.nc.variables['PHN'][:]
		if self.pft == 'Diatoms':		
			self.c   = self.nc.variables['PHD'][:]
		if self.pft == 'NonDiatoms':		
			self.c   = self.nc.variables['PHN'][:]
		
		self.c = self.c*79.573
	self.dataloaded = True
	
  def makePlot(self,region='Surface',doLog=True,fullrange=False):
  	if fullrange:
	  	self.autoFilename('cchl','cchl_'+region+'_'+self.pft)
	else:
	  	self.autoFilename('cchl','cchl_'+region+'_'+self.pft+'_tight')	
	  		
  	if self.existsFilename():
  		print "cchl:\tINFO\tPlot already exists:", self.filename  
  		return
  		
	#self.filename = folder('images/analysis')+'_'.join(['analysis',self.date,self.jobID])+'.png'
  	print "cchl:\tINFO\tMaking plot:", self.filename 	
  	fig = pyplot.figure()
  	ax = fig.add_subplot(111)
  	
  	#pyplot.pcolormesh

  	#etc
  	if not self.dataloaded:self.loadcchlData()
  	c   = sliceA(self.c,region)
  	chl = sliceA(self.chl,region)
  	

	c    = np.ma.masked_where((c>10E10) +(c<10E-05),c)
  	if fullrange:	
		chl  = np.ma.masked_where((chl>10E10) +(chl<10E-05),chl)

  	if not fullrange:
		chl  = np.ma.masked_where((chl>20.) +(chl<0.01),chl)  	


	m = c.mask + chl.mask	
	c = 	np.ma.masked_where(m,c).compressed()
	chl = 	np.ma.masked_where(m,chl).compressed()


	if doLog:
		# find the model fit:
		
		x_range = [chl.min(), chl.max()]

		x_r_log10 = np.log10(x_range)
		dxnolog =  np.linspace(x_range[0], x_range[1],50.)	
		fitdy = fitAndplot(chl, c, dxnolog)
			
		c = np.log10(c)
		chl= np.log10(chl)
		pyplot.xlabel(r'log$_{10}$ '+self.pft+' Phytoplankton Carbon')
		pyplot.ylabel(r'log$_{10}$ Total Chlorophyll')

		dx  = np.linspace(x_r_log10.min(), x_r_log10.max(),50.)	
		bdx = np.linspace(x_r_log10.min(), 0.,50.)

		
		pyplot.plot(dx,  C2Chl.SH(dx), label = 'Sathyendranath 2009 (HPLC)')
		pyplot.plot(dx,  C2Chl.ST(dx), label = 'Sathyendranath 2009 (Turner)')
		pyplot.plot(bdx,  C2Chl.B(bdx),  label = 'Buck (Turner)')
		

		pyplot.plot(np.ma.log10(dxnolog), np.ma.log10(fitdy), 'r--',lw=2, label = 'Fit to '+self.model)		
		
		#if drawMinAndMax:
		pyplot.plot(np.ma.log10(dxnolog),  np.ma.log10( 15.*dxnolog), 'k--',)# label = 'Minimum')
		pyplot.plot(np.ma.log10(dxnolog),  np.ma.log10(176.*dxnolog), 'k--',)# label = 'Maximum')
							
		pyplot.tight_layout()			
		leg = pyplot.legend(loc='best',prop={'size':10})
		leg.draw_frame(False) 
		leg.get_frame().set_alpha(0.)		
		
						
	else:
		pyplot.xlabel(self.pft+' Phytoplankton Carbon')
		pyplot.ylabel('Total Chlorophyll')
				
	
  	h = pyplot.hexbin(chl, c,  bins='log',mincnt=1,cmap=pyplot.get_cmap('Greys')) #xscale='log', yscale='log',
  	pyplot.colorbar()
  	
  	pyplot.title(self.model +' '+self.date+' '+region+' '+self.pft+' Carbon: Total Chlorophyll')

  	print "analysis:\tINFO\tSaving: " + self.filename
	pyplot.savefig(self.filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	






def main():
	try:
		filesIn = getFileList(argv[1:])
		if len(filesIn)==0:assert False
		print "Using command line arguments:", filesIn	
		
	except:
		print "cchl:\tERROR:\tNo file specified"
		#return
		filesIn = ['/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc',]
		
	
	for fn in filesIn:
		a = cchl(fn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
	
	
