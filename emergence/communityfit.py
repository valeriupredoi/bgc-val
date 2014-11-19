#!/usr/bin/ipython 
#Standard Python modules:
from netCDF4 import Dataset
from sys import argv
from matplotlib import pyplot, colors
import numpy as np
#from scipy import optimize
#Specific local code:
#from pftnames import pftnames as pn
from UKESMpython import folder,getFileList,sliceA
from analysis import analysis

from bgcvaltools import communitystructure as cs
#try:
#	import communitystructure as cs
#except:#
#	print "THIS REQUIRES The Community structure MODULE."
#	print "Please download from:"
#	print "https://gitlab.ecosystem-modelling.pml.ac.uk/ledm/communitystructure"
#	exit


	
	
class communityfit(analysis):
  def __init__(self,fileIn):
	analysis.__init__(self, fileIn)
	print self.fileIn, self.jobID
	self.dataLoaded= False
	#if self.model == 'ERSEM':
	pfts = ['Total', 'Diatoms', 'NonDiatoms' ]

	
	
	for r in ['SurfaceNoArtics','Surface', 'Top40m','Top200m', 'Global']:
	   for n in [2,3,4]:
		self.makePlot(region = r,numberPlots=n)


  def loadcommunityfitData(self,):
  	"""	While this method can be memory intensive, it may be faster to do it this way than to load them from file each time.
  	"""
  	print "communityfit:\tINFO:\tloadData"
  	if not self.ncOpened: self.loadData()

	#c and chl have different names in each model
	if self.model=='ERSEM':
  		print "communityfit:\tINFO:\tloading ERSEM data"
		self.dia    = self.nc.variables['Chl1'][:]
		self.flag   = self.nc.variables['Chl2'][:]
		self.pico   = self.nc.variables['Chl3'][:]		
		self.large  = self.nc.variables['Chl4'][:]
			

		
	if self.model=='MEDUSA':
  		print "communityfit:\tINFO:\tloading MEDUSA data"				
		self.dia   = self.nc.variables['CHD'][:]#79.573
		self.nond  = self.nc.variables['CHN'][:]#79.573
	self.dataLoaded = True

	
  def makePlot(self,region='Surface',numberPlots=2,doLog=True):
  	
  	self.autoFilename('communityfit','communityfit_'+region)
  	
  	if self.existsFilename():
  		print "communityfit:\tINFO\tPlot already exists:", self.filename  
  		return
  	print "communityfit:\tINFO\tMaking plot:", self.filename 	
  	  		
	if not self.dataLoaded: self.loadcommunityfitData()
	
	x_range = [-2.,1.]
	drawRunningMean = True
	drawfit = True
	fx = np.arange(x_range[0],x_range[1],0.1)
	fits = ["Brewin 2014", 'Hirata 2011',"Devred 2011","Brewin 2012",]#"Brewin 2011a","Brewin 2010",	
	

	if self.model == 'MEDUSA':
		if numberPlots!=2:return
	  	print "communityfit:\tINFO\tMaking slices:", region, self.model
	  	
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
		#fits:
		if drawfit:
			csfs = {}
			csfs['Non Diatoms'] = cs.comstrucFit(	totalc,		# Total Carbon
							carbonpc[1]/100., 	# PFT carbon
							pft='piconano',		# pft type
							chlRange=[10.**x_range[0],10.**x_range[1]],# range
							fitTypes =[])
			fits_xvalues= [	fx, fx]	
			tmpfy = csfs['Non Diatoms'].plotForX(10.**fx) *100.
			fits_yvalues= [	100. -tmpfy , tmpfy ]
		legx,legy = 0.475, -0.18		
			

		
	if self.model == 'ERSEM':
  	    print "communityfit:\tINFO\tMaking slices:", region, self.model
	    if numberPlots==2:	  	
	  	micro     = sliceA(self.dia , region) + sliceA(self.large , region)
	  	piconano  = sliceA(self.pico,region)  + sliceA(self.flag,region)
		m = micro.mask + piconano.mask
	  	micro  = np.ma.masked_where(m, micro ).compressed()
	  	piconano = np.ma.masked_where(m, piconano).compressed()
  		totalc = micro+piconano		
	  	carbonpc = [100.*micro/totalc,100.*piconano/totalc]
		#totalc = 
	  	titles = ['Microphytoplankton', 'Pico and Nano phytoplankton',]
		subplots = [211,212]
		#fits:
		if drawfit:
			csfs = {}
			csfs['Pico and Nano phytoplankton'] = cs.comstrucFit(	totalc,			# Total Carbon
							carbonpc[1]/100., 	# PFT carbon
							pft='piconano',		# pft type
							chlRange=[10.**x_range[0],10.**x_range[1]],# range
							fitTypes =[])
			fits_xvalues= [	fx, fx]	
			tmpfy = csfs['Pico and Nano phytoplankton'].plotForX(10.**fx) *100.
			fits_yvalues= [	100. -tmpfy , tmpfy ]
		legx,legy = 0.475, -0.18		
            else: return
  	

	
  	fig = pyplot.figure()
  	for i,sp in enumerate(subplots):
	  	ax = fig.add_subplot(sp)
	  	
		if i ==0: 	
		  	pyplot.title(self.model +' '+self.date+' '+region+' Community Structure')	  	
	  	if sp == subplots[-1]:	
	  		pyplot.xlabel(r'log$_{10}$ Phytoplankton Chl')
	  		
		pyplot.ylabel(titles[i]+' chl %')

		# Draw black to white distribution
		pyplot.hist2d(np.ma.log10(totalc),carbonpc[i],bins=(1000,52),cmin=1, norm=colors.LogNorm(),cmap=pyplot.get_cmap('Greys'))	

		# Draw plot range:
		pyplot.xlim(x_range)		
		pyplot.ylim([-4., 102.])
		
		# Drawing running means
		if drawRunningMean:
			csf = cs.comstrucFit(totalc,carbonpc[i], pft='pico', chlRange=[10.**x_range[0],10.**x_range[1]], fitTypes =['NoFit',])
			cx,cy = csf.doRunningMean()		
			cx = np.ma.log10(cx)
			pyplot.plot(cx,cy, 'k-', label='Running mean')		
		
		# Drawing published fits.
		for fit in fits:	
			if titles[i].lower().replace(' ','') in ['nondiatoms', 'piconano', 'picoandnanophytoplankton', 'picoandnanophytoplankton']:
				fy = cs.chlPercent(10.**fx, 'piconano',fit=fit)
				pyplot.plot(fx,fy,  label=fit)				
			if titles[i].lower().replace(' ','') in ['diatoms', 'micro','microphytoplankton',]:			
				fy = cs.chlPercent(10.**fx, 'micro',fit=fit)			
				pyplot.plot(fx,fy,  label=fit)
				
		# Drawing fit to model data:
		if drawfit:
			pyplot.plot(fits_xvalues[i],fits_yvalues[i],  'r--', label='Fit to '+self.model)

		# Drawing legend on last pane:						
		if sp == subplots[-1]:
				#legnd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=(len(fits)+1+1)/2,prop={'size':11},borderaxespad=0.)
				legnd = ax.legend(loc='upper center', bbox_to_anchor=(legx,legy), ncol=len(fits)+2,prop={'size':8.5},borderaxespad=0.)
				legnd.draw_frame(False) 
				legnd.get_frame().set_alpha(0.)
					  	


  	print "communityfit:\tINFO\tSaving: " + self.filename
	pyplot.savefig(self.filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	






def main():
	try:
		filesIn = getFileList(argv[1:])
		if len(filesIn)==0:assert False
		print "Using command line arguments:", filesIn	
		
	except:
		print "communityfit:\tWARNING:\tNo file specified"
		#return
		filesIn = ['/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc','/data/euryale7/scratch/ledm/iMarNet/xhonv/MEANS/xhonvo_18990901m01P.nc',]
		
	
	for fn in filesIn:
		a = communityfit(fn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
	
	
