#!/usr/bin/ipython 
#
# Copyright 2014, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license. 

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
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
from string import ascii_lowercase
from re import findall
import os
#Specific local code:
#from bgcvaltools.pftnames import pftnames as pn
from UKESMpython import folder,getFileList


"""	This is a template for performing an analsysis based on a single output file.
"""

class analysis:
  def __init__(self,fileIn,):#makePlot=True):
  	self.fileIn = fileIn
	self.ncOpened = False
  	self.filenamed = False	  		
	if not self.ncOpened: self.loadData()
	  	
	self.autoDate()
	self.autoJobID()
	self.autoModel()
	self.autoCalendar()		
	#self.autoFilename()

	
	#if makePlot:self.makePlot()
	
  def __del__(self):
  	# Destructor method.
  	print "analysis:\tINFO\tClosing netcdf"
	self.nc.close()
	self.ncOpened=False
	try:
		pyplot.close()
	except: pass

	
  def loadData(self,):
  	"""	Load the netcdf given.
  	"""
  	print "analysis:\tINFO\tLoading netcdf"  
  	self.nc = Dataset(self.fileIn,'r')
  	self.ncOpened = True
  	# Loading netcdf data with pftnames:
	#self.chl = self.nc.variables[pn['E']['diatoms']['chl']][0,0,]

  def loadFile(self,filename):
  	"""	Load the netcdf given.
  	"""
  	print "analysis:\tINFO\tLoading netcdf file" , filename
  	return Dataset(filename,'r')




  def autoDate(self,):
  	self.date	= 'date'	  
	a = findall(r'\d\d\d\d\d\d\d\d',self.fileIn)
	for i in a:
	    if 17000000<int(i) <22000000:
	      self.date = i
	      return
	      	
	a = findall(r'\d\d\d\d',self.fileIn)
	for i in a:
	    if 1700<int(i) <2200:
	      self.date = i
	      return

	if fn.find('clim')>-1:
		self.date = 'Climatology'
		return	      		
	#if not self.ncOpened: self.loadData()	





  def autoJobID(self,):
	self.jobID = 'jobID'
  	for jobid in ['xhon','xjez','xkad']:
  	    for a in ascii_lowercase:
  		if self.fileIn.find(jobid+a)>-1:  	
  			self.jobID = jobid+a
  			return	  
	self.jobID = 'iMarNet'
	  	
  	
  def autoModel(self,):
  	for jobid in ['xhon','xjez','xkad']:
  		if self.fileIn.find(jobid)>-1:  	
  			self.model	= 'ERSEM'		
  			return 
  	self.model	= 'MEDUSA' 	  			

  def autoCalendar(self,):
  	for jobid in ['xhon',]:#'xjez','xkad']:
  		if self.fileIn.find(jobid)>-1:  	
  			self.cal	= '365_day'		
  			return 
  	for jobid in ['xjez',]:
  		if self.fileIn.find(jobid)>-1:  	
  			self.cal	= '360_day'		
  			return 
  	self.cal	= '365_day'
  			  			

  	  	
  def autoFilename(self,plotType,plotType2,noDate=False):
  	if noDate:
 	 	self.filename =  folder(['images',self.model,self.jobID,plotType])
	  	self.filename += '_'.join([plotType,self.model,self.jobID,plotType2])+'.png'
  	else:
 	 	self.filename =  folder(['images',self.model,self.jobID,plotType,self.date])
	  	self.filename += '_'.join([plotType,self.model,self.jobID,self.date,plotType2])+'.png'
  	print "analysis:\tINFO:\t Auto Filename:",self.filename
  	self.filenamed = True	  	

  def existsFilename(self):
	if self.filenamed == False: 		return False
  	if os.path.exists(self.filename): 	return True
  	

  	  	
  def makePlot(self,):
  	"""	This is a dummy routine to produce a blank plot,
  		as an example with some useful bits and pieces.
  	"""
  	print "analysis:\tWARNING:\t analysis.makePlot only makes a dummy plot!" 
  	
  	 
  	self.autoFilename('analysis','')
	#self.filename = folder('images/analysis')+'_'.join(['analysis',self.date,self.jobID])+'.png'
  	print "analysis:\tINFO\tMaking plot:", self.filename 	
  	fig = pyplot.figure()
  	ax = fig.add_subplot(111)
  	
  	#pyplot.pcolormesh
  	#pyplot.hexbin 
  	#etc
  	
  	pyplot.title('Title')
	pyplot.xlabel('x axis label')
	pyplot.ylabel('y axis label')
		
  	print "analysis:\tINFO\tSaving: " + self.filename
	pyplot.savefig(self.filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	








def main():
	try:
		filesIn = getFileList(argv[1:])
		print "Using command line arguments:", filesIn	
		
	except:
		print "analysis:\tERROR:\tNo file specified"
		return
		
	
	for fn in filesIn:
		a = analysis(fn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
