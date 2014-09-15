#!/usr/bin/ipython 
from netCDF4 import Dataset

from pftnames import pftnames
from UKESMpython import getCommandJobIDandTime


"""	This is a template for performing an analsysis based on a single output file.
"""

class analysistemplate:
  def __init__(self,fileIn,)
	getDateFromFile(filein)
	getJobIDFromFile(filein)

	self.loadData()
	self.makePlot()
	
  def loadData(self,):
  	self.nc = Dataset(self.fileIn,'r')


  def makePlot(self,):
  	fig = pyplot.figure()
  	ax = fig.add_subplot(111)
  	
  	pyplot.title('Title')

	pyplot.xlabel('x axis label')
	pyplot.ylabel('y axis label')
	
		
  	print "analysistemplate:\tmakePlot\tSaving: " + self.filename
	pyplot.savefig(filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	





































	
if __name__=="__main__":
	main() 
	print 'The end.'
