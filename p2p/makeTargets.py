#!/usr/bin/ipython
from matplotlib import pyplot
from matplotlib import rc
from matplotlib.patches import Arrow
#from matplotlib.markers import filled_markers 
from StatsDiagram import TaylorDiagram, TargetDiagram
from glob import glob

import numpy as np
#from numpy import sin, arccos, sqrt, max as nmax, abs
#from numpy import isnan as npNaN, isinf as npInf
#from calendar import month_name

from itertoolsmodule import product
from itertools import cycle
from operator import itemgetter
from os.path import basename
from sys import argv
from shelve import open as shOpen

import UKESMpython as ukp
from pftnames import AutoVivification



#from deMoraTools import folder, foldToEps,AutoVivification,mnStr
#from shelveManager import shelveManager
#from BGCnames import fullname
#from iMarNetPython import getLongName,shouldIMakeFile



class makeTargets:
  def __init__(self,matchedShelves, filename,name='', imageDir='', diagramTypes=['Target',],debug=True): 

  
  	self.matchedShelves =matchedShelves
    	self.name = name
    	self.filename = filename
	self.diagramTypes = diagramTypes 	#['Taylor','Target']
	self.debug = debug
	self.imageDir = imageDir



	
	
	
  	#self.shelvedir = workingDir
  	#if self.shelvedir == '':self.shelvedir = ukp.folder(['shelves',self.xtype,self.ytype, 'Slices',self.name])
  	#else:			self.shelvedir = ukp.folder(self.shelvedir)		

	self.loadShelves()
	
	self.makeDiagram()


  def loadShelves(self,):
  	self.data = AutoVivification()
  	
  	for sh in self.matchedShelves:
  		s = shOpen(sh,flag='r')
 		E0 = s['Taylor.E0' ]
		R  = s['Taylor.R']
		G  = s['Taylor.gamma']
		p  = s['Taylor.p']	 		
		N  = s['N']
		leg = ' - '.join([s[i] for i in ['xtype','ytype', 'name', 'newSlice','xkey','ykey',]])
		#leg =  mnStr(i)+'/'+getLongName(leg[i])
		try:	self.xtype[s['xtype']]=True
		except: 
			self.xtype = {s['xtype']:True}
			
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
	
	self.xtype = ''.join(self.xtype.keys())
				
		
  def makeDiagram(self):
	title = 'Model vs Data'
	filled_markers  =  ('o', 'v', '^', '<', '>', '8', 's', 'p',  'h',   'd')#'*','H','D',	
	markercycler = cycle(filled_markers )

	if self.imageDir=='':	self.imageDir = ukp.folder(['images',self.xtype,'P2P_plots',self.name])
	else: 			self.imageDir = ukp.folder(self.imageDir)
		

	for t in self.diagramTypes:
		fig = pyplot.figure()		
		ax = pyplot.subplot(111, aspect='equal')
		c = pyplot.get_cmap('jet')
		
		proxyArt,labs = [],[]
		for leg in sorted(self.data.keys()):
		    ma = next(markercycler)
		    
		    try:
		    	TD.add(self.data[leg]['G'], self.data[leg]['E0'],self.data[leg]['R'], marker = ma, s=150, cmap=c, label=leg,)
		    	#TD.add(g[i], E0[i], R[i],  marker = ma, s=150, cmap=c, label=i,)
		    	#TD.labels(i)
		    except:
		    	print 'First target diagram:\t', leg
		    	TD=TargetDiagram(self.data[leg]['G'], self.data[leg]['E0'],self.data[leg]['R'],marker = ma,s=150,cmap=c, label=leg,)
		    labs.append(leg)
		    proxyArt.append(pyplot.Line2D([0],[0], linestyle="none", c=c(self.data[leg]['R']), marker = ma,markersize=9,))

		    
		    #print leg.ljust(maxForTable+1), '\tgamma:',round(g[i],4),'\tE0:',round(E0[i],4),'\tR:',round(R[i],4),'\tN:',N[i].rjust(maxN+1),'\tp:',round(p[i],6)


		#w = 0.33
		#rmax= 2.25
		#if t == 'Target' and max(abs([pyplot.xlim().max(),pyplot.clim().max()]))<rmax:			
		#	pyplot.xlim([-rmax,rmax])
		#	pyplot.ylim([-rmax,rmax])
		        #pyplot.plot((0,0),(-rmax,rmax),'k-')
			#pyplot.plot((rmax,-rmax),(0,0),'k-')

		#box = ax.get_position() 			
		#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
		#pyplot.legend()

		doLabels = True # False
		if doLabels:
			#labs = [l.split('/')[1] for l in labs]		
			if len(labs)<8:
				legnd = pyplot.legend(proxyArt, labs,loc=4, borderaxespad=0., numpoints = 1, ncol=1, scatterpoints=1, prop={'size':8},) #title=legtitle) #
			else:
				legnd = pyplot.legend(proxyArt, labs,loc=8, borderaxespad=0., numpoints = 1, ncol=2, scatterpoints=1, prop={'size':8},) #title=legtitle) #		
		
		
		legnd.draw_frame(False) 
		legnd.get_frame().set_alpha(0.) #v3
		pyplot.title(title)
		
		if self.debug:print 'multiTarget:\tsaving file:', self.filename 
		pyplot.savefig(self.filename ,dpi=300,bbox_inches='tight')
		pyplot.close()
		del(TD)
		del(fig)
		
		







		
if __name__=="__main__":
	print 'The end.'








