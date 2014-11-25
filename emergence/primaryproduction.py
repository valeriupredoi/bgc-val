#!/usr/bin/ipython 

#Standard Python modules:
import os
from sys import argv
from netCDF4 import Dataset,num2date
from matplotlib import pyplot
import numpy as np
from glob import glob
from shelve import open as shOpen

#Specific local code:
from UKESMpython import folder,getFileList,sliceA,makeOneDPlot
from analysis import analysis


class primaryproduction(analysis):
  def __init__(self,filesIn):
  	fns = []
  	for f in filesIn:
  		fns.extend(glob(f))
	analysis.__init__(self, sorted(fns)[0])
	
	self.files = fns
	
	self.loadPProdData(region = 'Global')
	self.makeOneDPlot(region = 'Global')

  def loadPProdData(self,region='Global'):
	shelveFile = folder('shelves/primaryproduction/'+self.model+'/'+self.jobID)+self.jobID+'_'+region+'.shelve'	


	try:
		s = shOpen(shelveFile)
		fnDictsT	= s['fnDictsT']
		fnDictsPP	= s['fnDictsPP']				
		readFiles	= fnDictsPP.keys()
		s.close()
		savedLength = len(readFiles)
		print "makePProdPlot:\tOpened shelve:", shelveFile, '\tread'#, len(times)
		
	except:
		readFiles= []
		fnDictsPP,fnDictsT = {},{}
		savedLength = 0		
		print "makePProdPlot:\tcreating new shelve:", shelveFile
		
	#####
	# Load cell volume
	# Warning, this will fail if you're using a different grid. (ORCA1)
  	ncmesh = Dataset("data/mesh_mask_ORCA1_75.nc",'r')
	pvol =  sliceA(ncmesh.variables['pvol'][:],region)
	ncmesh.close()

	SaveEvery,count =12,0
	for fn in self.files[:]:
		if fn not in readFiles: 
		    print 'run:', fn
		    nc = Dataset(fn,'r')
		
		    dtime = num2date(nc.variables['time_counter'][:], nc.variables['time_counter'].units,calendar=self.cal)
		    for t,dt in enumerate(dtime):
			if self.model=='ERSEM':
				npp = sliceA(nc.variables['netPP'][t,],region)
				pp = np.ma.masked_where(npp > 1.E30,npp).sum()
				pp = pp/1000000000. 						# tons -> Gigatons
				pp = pp*365. 							# /day -> /year	
			if self.model=='MEDUSA':
				npp = sliceA(nc.variables['TOTPP3'][t,:],region)		# [gC/m3/d]
				pp  = np.ma.masked_where(npp > 1.E35,npp*pvol).sum() 		# [gC/d]
				pp  = pp*365. 							# [gC/year]
				pp  = pp/1.E15							# [GT/year]
				
			try:	fnDictsPP[fn].append(pp)
			except: fnDictsPP[fn] = [pp,]
			try:	fnDictsT[fn].append(dt)
			except: fnDictsT[fn] = [dt,]

	 		if fn not in readFiles: readFiles.append(fn)
			count+=1
			if count%SaveEvery ==0:
				print "primProd:\tSaving shelve:", shelveFile, '\tread'#, len(netPPs)									
				s = shOpen(shelveFile)			
				s['fnDictsPP'] = 	fnDictsPP			
				s['fnDictsT'] = 	fnDictsT							
				s.close()
				count=0	
		    nc.close()				

	if savedLength < len(self.files):
		print "primProd:\tFinal saving shelve:", shelveFile, '\tread', len(fnDictsPP.keys())									
	  	s = shOpen(shelveFile)	
		s['fnDictsPP'] 	= fnDictsPP			
		s['fnDictsT'] 	= fnDictsT
		s.close()
		
	netPPs = []
	times  = []		
	for fn in sorted(fnDictsPP.keys()):
		print fn,fnDictsPP[fn], fnDictsT[fn]
		netPPs.extend(fnDictsPP[fn])
		times.extend(fnDictsT[fn])



	annual 		= {}
	years 		= []
	annualNetPP 	= [] 
	for t,n in zip(times,netPPs):
		try:	annual[t.year].append(n)
		except: annual[t.year] = [n,]
	
	for k,i in sorted(annual.items()):
		years.append(k+0.5)
		if len(i)==12:
			annualNetPP.append(np.ma.array(i).mean())
		elif len(i)==1:
			annualNetPP.append(i)
			print k, 'One file per year'	
		else:
			annualNetPP.append(np.ma.masked)
			print k, 'Not enough months for annual mean:', len(i)	

	self.times 		= np.array(times).squeeze()
	self.netPPs 		= np.array(netPPs).squeeze()
	self.years 		= np.array(years)
	self.annualNetPP	= np.array(annualNetPP)

	print self.years, self.annualNetPP
	
	if len(self.times) != len(self.netPPs):
		print "makeOneDPlot:\tTHere is a size Mismatch between time and data", len(self.times),'!=', len(self.netPPs)
		assert False	
	
	if len(self.years) != len(self.annualNetPP):
		print "makeOneDPlot:\tTHere is a size Mismatch between years and annualNetPP", len(self.years),'!=', len(self.annualNetPP)
		assert False	

	




  def makeOneDPlot(self,region='Global',dpi=100, drawRange=True):
	analysis.autoFilename(self,'primaryproduction', 'primaryproduction_'+region,noDate=True)
	
	print "makeOneDPlot: ", self.filename
	
	fig = pyplot.figure()
	ax = fig.add_subplot(111)
	fig.set_size_inches(16, 6)
	
	print self.times,self.netPPs
	newTimes  = [t.year + (t.month-1)/12. for t in self.times]
	pyplot.plot(newTimes, self.netPPs,c='k',ls='-')
	pyplot.plot(self.years, self.annualNetPP,'b-')
	
	#pyplot.plot_date(x=self.times, y=self.netPPs,c='k',ls='-')
	if drawRange:
		ymin,ymax = pyplot.ylim()
		ymin = min([ymin,35.])
		ymax = max([ymax,75.])
		pyplot.ylim((ymin,ymax))
			
		xlims  = pyplot.xlim()
		#pyplot.axhline(y=45.,c='r',ls='-',lw=3,alpha=0.4)
		#pyplot.axhline(y=65.,c='r',ls='-',lw=3,alpha=0.4)
		lll = np.array([ymin for i in xlims]) 
		l40 = np.array([40. for i in xlims])
		l50 = np.array([50. for i in xlims])
		l60 = np.array([60. for i in xlims])		
		l70 = np.array([70. for i in xlims])
		lul = np.array([ymax for i in xlims]) 		

		ax.fill_between(xlims,lll, l40 ,color='r', alpha = 0.2)		
		ax.fill_between(xlims,l40 ,l50 ,color='DarkOrange', alpha = 0.2)				
		ax.fill_between(xlims,l50 ,l60 ,color='g', alpha = 0.2)
		ax.fill_between(xlims,l60 ,l70 ,color='DarkOrange', alpha = 0.2)
		ax.fill_between(xlims,l70 ,lul ,color='r', alpha = 0.2)		
		
		#pyplot.axhline(y=50.,c='g',ls='-',lw=2,alpha=0.5)
		#pyplot.axhline(y=60.,c='g',ls='-',lw=2,alpha=0.5)		

#	title = region + ' Primary Production - Gigatons / year'						
	title	= ' '.join([region,self.model,'Primary Production, GT/year'])
	
	pyplot.title(title) 
		
	print "makeOneDPlot:\tSaving: " +  self.filename
	pyplot.savefig( self.filename,dpi=dpi)
	pyplot.close()	


def main():

	try:
		filesIn = getFileList(argv[1:])
		if len(filesIn)==0:assert False
		print "Using command line arguments:", filesIn	
		
	except:
		print "cchl:\tERROR:\tNo file specified"
		#return
		filesIn = ['/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_199?.nc',]
	
	a = primaryproduction(filesIn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
	
	
	
