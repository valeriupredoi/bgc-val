#!/usr/bin/ipython 

#Standard Python modules:
import os
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

	#jobID = getJobID(files[0],False)
	#cal = getCalendar(files[0])
	#print jobID
	shelveFile = folder('shelves/primaryproduction/'+self.model+'/'+self.jobID)+self.jobID+'_'+region+'.shelve'
	
	annual = {}

	try:
		s = shOpen(shelveFile)
		fnDictsT	= s['fnDictsT']
		fnDictsPP	= s['fnDictsPP']				
		#netPPs 	= s['netPPs']
		#times 		= s['times']
		readFiles	= fnDictsPP.keys()
		s.close()
		savedLength = len(readFiles)
		print "makePProdPlot:\tOpened shelve:", shelveFile, '\tread'#, len(times)
		#fnDictsPP,fnDictsT = {},{}
		#for fn in readFiles:	
		#	fnDictsPP[fn] = netPPs[fn]
		#	fnDictsT[fn]  = times[fn]			
			

	except:
		#netPPs = {}
		#times = {}
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
			
		    dtime = num2date(nc.variables['time_counter'][:], nc.variables['time_counter'].units,calendar='365_day')
		    for t,dt in enumerate(dtime):
			#yr= dtime.year
			if self.model=='ERSEM':
				npp = sliceA(nc.variables['netPP'][t,],region)
				pp = np.ma.masked_where(npp > 1.E30,npp).sum()
				pp = pp/1000000000. 						# tons -> Gigatons
				pp = pp*365. 							# /day -> /year	
			if self.model=='MEDUSA':
				npp = sliceA(nc.variables['TOTPP3'][t,:],region)		# [gC/m3/d]
				pp =  np.ma.masked_where(npp > 1.E35,npp*pvol).sum() 		# [gC/d]
				pp = pp*365. 							# [gC/year]
				pp = pp/10.E15							# [GT/year]
			#time = dt	#dtime.year + (dtime.month-1)/12.
			try:	fnDictsPP[fn].append(pp)
			except: fnDictsPP[fn] = [pp,]
			try:	fnDictsT[fn].append(dt)
			except: fnDictsT[fn] = [dt,]			
			
			#netPPs.append(pp)
			#times.append(dt)

	 		if fn not in readFiles: readFiles.append(fn)
			count+=1
			if count%SaveEvery ==0:
				print "primProd:\tSaving shelve:", shelveFile, '\tread'#, len(netPPs)									
				s = shOpen(shelveFile)			
				#s['netPPs'] = netPPs
				#s['times']  = times
				#s['readFiles']= readFiles
				s['fnDictsPP'] = 	fnDictsPP			
				s['fnDictsT'] = 	fnDictsT							
				s.close()
				count=0	
		    nc.close()				

		#pps = fnDictsPP[fn] 
		#times = fnDictsT[fn] 
		
		#yr = np.mean([t.year for t in times])
		#try:annual[yr] = pps
		#except:annual[yr] = [pp,]
		#print os.path.basename(fn),time,pp,yr
		
		
	if savedLength < len(self.files):
		print "primProd:\tFinal saving shelve:", shelveFile, '\tread', len(fnDictsPP.keys())									
	  	s = shOpen(shelveFile)	
		#s['netPPs'] = netPPs
		#s['times']  = times
		s['fnDictsPP'] = 	fnDictsPP			
		s['fnDictsT'] = 	fnDictsT		
		#s['readFiles']= readFiles
		s.close()
		
	netPPs = []
	times = []		
	for fn in sorted(fnDictsPP.keys()):
		print fn,fnDictsPP[fn], fnDictsT[fn]
		netPPs.extend(fnDictsPP[fn])
		times.extend(fnDictsT[fn])




	years 		= []
	annualNetPP 	= [] #
	for k,i in sorted(annual.items()):
		years.append(k+0.5)
		if len(i)==12:
			annualNetPP.append(array(i).mean())
		elif len(i)==1:
			annualNetPP.append(i)
			print k, 'One file per year'	
		else:
			annualNetPP.append(masked)
			print k, 'Not enough months for annual mean:', len(i)	

	self.times 		= np.array(times).squeeze()
	self.netPPs 		= np.array(netPPs).squeeze()
	self.years 		= np.array(years)
	self.annualNetPP	= np.array(annualNetPP)
	
	
	if len(self.times) != len(self.netPPs):
		print "makeOneDPlot:\tTHere is a size Mismatch between time and data", len(self.times),'!=', len(self.netPPs)
		assert False	
	
	if len(self.years) != len(self.annualNetPP):
		print "makeOneDPlot:\tTHere is a size Mismatch between years and annualNetPP", len(self.years),'!=', len(self.annualNetPP)
		assert False	

	




  def makeOneDPlot(self,region='Global',dpi=100):
	analysis.autoFilename(self,'primaryproduction', 'primaryproduction_'+region,noDate=True)
	#filename = folder(['images',self.model,self.jobID,'PrimaryProduction'])+'primaryproduction_'+region+'.png'
	
	print "makeOneDPlot: ", self.filename
	
	fig = pyplot.figure()
	ax = fig.add_subplot(111)
	fig.set_size_inches(16, 6)
	
	print self.times,self.netPPs
	newTimes  = [t.year + (t.month-1)/12. for t in self.times]
	pyplot.plot(newTimes, self.netPPs,c='k',ls='-')	
	#pyplot.plot_date(x=self.times, y=self.netPPs,c='k',ls='-')
	#pyplot.plot(self.years, self.annualNetPP,'b-')	
	
	#pyplot.axhline(y=45.,c='r',ls='--')
	#pyplot.axhline(y=65.,c='r',ls='--')	

#	title = region + ' Primary Production - Gigatons / year'						
	title	= ' '.join([region,self.model,'Primary Production, GT/year'])
	
	pyplot.title(title) 
		
	print "makeOneDPlot:\tSaving: " +  self.filename
	pyplot.savefig( self.filename,dpi=dpi)#, bbox_inches='tight')
	pyplot.close()	

			
