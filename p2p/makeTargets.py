#!/usr/bin/ipython
from matplotlib import pyplot
from matplotlib import rc
from matplotlib.patches import Arrow
#from matplotlib.markers import filled_markers 
from glob import glob

import numpy as np
#from numpy import sin, arccos, sqrt, max as nmax, abs
#from numpy import isnan as npNaN, isinf as npInf
#from calendar import month_name

#from itertoolsmodule import product
from itertools import cycle, product
from operator import itemgetter
from os.path import basename,exists
from sys import argv
from shelve import open as shOpen

import UKESMpython as ukp
from pftnames import AutoVivification,getLongName

from bgcvaltools.StatsDiagram import TaylorDiagram, TargetDiagram,TaylorDiagramMulti





class makeTargets:
  def __init__(self,matchedShelves, filename, diagramTypes=['Target','Taylor'],legendKeys = ['name', 'newSlice','xkey','ykey',],debug=True): #name='', #imageDir='',

  
  	self.matchedShelves =matchedShelves
    	#self.name = name

    	self.filename = filename
	self.diagramTypes = diagramTypes 	#['Taylor','Target']
	self.debug = debug
	#self.imageDir = imageDir

	self.legendKeys = legendKeys 
	self.determineLegend()



	
  	#self.shelvedir = workingDir
  	#if self.shelvedir == '':self.shelvedir = ukp.folder(['shelves',self.xtype,self.ytype, 'Slices',self.name])
  	#else:			self.shelvedir = ukp.folder(self.shelvedir)		

	if len(self.matchedShelves)>0 and ukp.shouldIMakeFile(self.matchedShelves,self.filename,debug=False):
		self.loadShelves()
		self.makeDiagram()

  def determineLegend(self,):
  	#####
  	# Determine legend:
  	# Looks for ways to differenciate the shelves by looking through their metadata. 
  	# It generally works, but it might be better just to declare what you want in the legend when you run it.
  	self.xtypes = {}
  	self.ytypes = {}
  	self.names = {}
  	self.regions = {}  	
  	self.years = {}
  	self.xkeys = {}
  	self.ykeys = {}  
  	self.newSlices={}
  	  	  	
  	for sh in self.matchedShelves:  	
  		print "determineLegend:\tINFO:\tLOADING:",sh
  		if not exists(sh): 
  			print "determineLegend:\tWARNING:\tDoes not exist:",sh 
  			continue    
  		s = shOpen(sh,flag='r')
	  	self.xtypes[s['xtype']]	= True
	  	self.ytypes[s['ytype']]	= True
	  	self.names[s['name']]	= True
	  	self.regions[s['region']]= True	  	
	  	self.years[s['year']]	= True	
	  	self.xkeys[s['xkey']]	= True	
	  	self.ykeys[s['ykey']]	= True		  		  	
		self.newSlices[s['newSlice']] = True  		
 		s.close() 

 	if len(self.legendKeys):return
 	
	self.legendKeys = []
	if len(self.names.keys()) >1:		self.legendKeys.append('name')
	if len(self.regions.keys()) >1:		self.legendKeys.append('region')
	if len(self.newSlices.keys()) >1:	self.legendKeys.append('newSlice')
	if len(self.xtypes.keys()) >1:		self.legendKeys.append('xtype')		
	if len(self.ytypes.keys()) >1:		self.legendKeys.append('ytype')		
	if len(self.xkeys.keys()) >1:		self.legendKeys.append('xkey')		
	if len(self.ykeys.keys()) >1:		self.legendKeys.append('ykey')	
	if len(self.years.keys()) >1:		self.legendKeys.append('year')				  	


  def loadShelves(self,):
  	self.data = AutoVivification()
  	
  	for sh in self.matchedShelves:
  		print "loadShelves:\tINFO:\tLOADING:",sh
  		if not exists(sh): 
  			print "loadShelves:\tWARNING:\tDoes not exist:",sh 
  			continue
  		s = shOpen(sh,flag='r')
 		E0 = s['Taylor.E0' ]
		R  = s['Taylor.R']
		G  = s['Taylor.gamma']
		p  = s['Taylor.p']	 		
		N  = s['N']
		leg = ' - '.join([getLongName(s[i]) for i in self.legendKeys])


								
		s.close()
		breaks=0
		for func,a in product([np.isnan,np.isinf],[E0,R,G,N,p],):
			if func(a):
				if self.debug: print 'LoadShelves:\tWARNING:\t',a, 'is nan/inf'
				breaks+=1
		if breaks>0: continue
		self.data[leg]['E0'] = E0
		self.data[leg]['R'] = R
		self.data[leg]['G'] = G
		self.data[leg]['p'] = p
		self.data[leg]['N'] = N
	
				
  def makeTitle(self,):
  	"""	MakeTitle determines how you have sliced the data.
  		ie, which ever field that there is only one of gets added to the title.
  		so if this is only from one year, or only one Model, then those are added to the title.
  	"""
  	
  	self.xtype = ', '.join(self.xtypes.keys())
	self.ytype = ', '.join(self.ytypes.keys())	
	title =self.xtype + ' Model vs '+self.ytype+' Data'
	if len(self.names.keys()) ==1:
		title += ', '+', '.join([getLongName(k) for  k in self.names.keys()])
	if len(self.regions.keys()) ==1:
		title += ', '+', '.join([getLongName(k) for  k in self.regions.keys()])
		
	if len(self.ykeys.keys()) ==1:
		title += ', '+', '.join([getLongName(k) for  k in self.ykeys.keys()])
		
	if len(self.years.keys()) ==1:
		title += ', '+ ', '.join([str(k) for  k in self.years.keys()]) 

	if len(self.newSlices.keys()) ==1:
		title = ', '.join([getLongName(k) for  k in self.newSlices.keys()]) +' '+ title
		
	print 'makeTitle:\t',title
	return title
			
  def makeDiagram(self):
  	title = self.makeTitle()

						
	filled_markers  =  ('o', 'v', '^', '<', '>', '8', 's', 'p',  'h',   'd')#'*','H','D',	
	markercycler = cycle(filled_markers )

	#if self.imageDir=='':	self.imageDir = ukp.folder(['images',self.xtype.replace(', ','-'),'Targets'])
	#else: 			self.imageDir = ukp.folder(self.imageDir)
		

	for t in self.diagramTypes:
		if not len(self.data.keys()):
			continue
		    	print 'makeDiagram\t:No Plots to make'

			
		fig = pyplot.figure()		
		ax = pyplot.subplot(111, aspect='equal')
		c = pyplot.get_cmap('jet')

		proxyArt,labs = [],[]
		
		if t == 'Target':
			for leg in sorted(self.data.keys()):
			    ma = next(markercycler)
			    
			    try:
			    	Target.add(self.data[leg]['G'], self.data[leg]['E0'],self.data[leg]['R'], marker = ma, s=150, cmap=c, label=leg,)
			    	#TD.add(g[i], E0[i], R[i],  marker = ma, s=150, cmap=c, label=i,)
			    	#TD.labels(i)
			    except:
			    	print 'makeDiagram:\tFirst target diagram:\t', title
			    	Target=TargetDiagram(self.data[leg]['G'], self.data[leg]['E0'],self.data[leg]['R'],marker = ma,s=150,cmap=c, label=leg,)
			    labs.append(leg)

			    print 'Target:\t',leg,'\tGamma:',self.data[leg]['G'], '\tE0:',self.data[leg]['E0'],'\tR:',self.data[leg]['R']			    
			    proxyArt.append(pyplot.Line2D([0],[0], linestyle="none", c=c(self.data[leg]['R']), marker = ma,markersize=9,))
				

			if len(labs)<8:
				legend = pyplot.legend(proxyArt, labs,loc=4, ncol=1, borderaxespad=0., numpoints = 1, scatterpoints=1, prop={'size':8},) 
			else:	legend = pyplot.legend(proxyArt, labs,loc=8, ncol=2, borderaxespad=0., numpoints = 1, scatterpoints=1, prop={'size':8},) 
	
		
		
		if t == 'Taylor':		
			#fig.set_size_inches(10,4)		# only if antiCorrelation = True
			gams = []
			e0s  = []
			Rs   = []
			marks = []
			for leg in sorted(self.data.keys()):
				gams.append(self.data[leg]['G'])
				e0s.append(self.data[leg]['E0'])
				Rs.append(self.data[leg]['R'])
				ma = next(markercycler)				
				marks.append(ma)
			    	labs.append(leg)

			    				
			
			Taylor=TaylorDiagramMulti(gams, e0s,Rs,antiCorrelation=False,markers = marks,s=150,cmap=c, label=leg,)

			cmax=max(1,int(np.abs(np.max(e0s))+1))
			
			for i,leg in enumerate(sorted(self.data.keys())):
			    	print 'Taylor:\t',leg,'\tGamma:',self.data[leg]['G'], '\tE0:',self.data[leg]['E0'],'\tR:',self.data[leg]['R']
			    	proxyArt.append(pyplot.Line2D([0],[0], linestyle="none", c=c((1./(2*cmax))*self.data[leg]['E0'] +0.5), marker = marks[i],markersize=9,))	
			    				


			    #try:
			    #	Taylor.add(self.data[leg]['G'], self.data[leg]['E0'],self.data[leg]['R'], marker = ma, s=150, cmap=c, label=leg,)

			    #except:
			    #	print 'makeDiagram\t:First Taylor diagram:\t', title
			    #	Taylor=TaylorDiagramMulti(self.data[leg]['G'], self.data[leg]['E0'],self.data[leg]['R'],R=5,antiCorrelation=False,marker = ma,s=150,cmap=c, label=leg,)
	    
				
			legend = pyplot.legend(proxyArt, labs,loc=1, ncol=1, borderaxespad=0., numpoints = 1, scatterpoints=1, prop={'size':8},) 		
		
		legend.draw_frame(False) 
		legend.get_frame().set_alpha(0.) 
		pyplot.title(title)
		filename = self.filename.replace('.png','_'+t+'.png')
		if self.debug:print 'makeDiagram:\tsaving file:', filename 
		pyplot.savefig(filename ,dpi=200,)
		pyplot.close()
		try:
			del(TD)
			del(fig)
		except: pass
		
		







		
if __name__=="__main__":
	print 'The end.'








