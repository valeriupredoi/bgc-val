#!/usr/bin/ipython 
#Standard Python modules:
from netCDF4 import Dataset
from sys import argv
from matplotlib import pyplot
#Specific local code:
from pftnames import pftnames
from UKESMpython import folder,getFileList


"""	This is a template for performing an analsysis based on a single output file.
"""

class analysistemplate:
  def __init__(self,fileIn,):
  	self.fileIn = fileIn
	self.date	= 'date'	#getDateFromFile(filein)
	self.jobID	= 'jobid'	#getJobIDFromFile(filein)

	self.loadData()
	self.makePlot()
	
	self.nc.close()
	
  def loadData(self,):
  	self.nc = Dataset(self.fileIn,'r')

  def makePlot(self,):
	
	self.filename = folder('images/analysistemplate')+'_'.join(['analysistemplate',self.date,self.jobID])+'.png'
	
  	fig = pyplot.figure()
  	ax = fig.add_subplot(111)
  	
  	pyplot.title('Title')

	pyplot.xlabel('x axis label')
	pyplot.ylabel('y axis label')
		
  	print "analysistemplate:\tINFO\tSaving: " + self.filename
	pyplot.savefig(self.filename,dpi=100,)# bbox_inches='tight')
	pyplot.close()	








def main():
	getFileList(argv[1:])
	try:
		filesIn = getFileList(argv[1:])
		print "Using command line arguments:", filesIn	
		
	except:
		print "analysistemplate:\tERROR:\tNo file specified"
		return
		
	
	for fn in filesIn:
		a = analysistemplate(fn)






	
if __name__=="__main__":
	main() 
	#print 'The end.'
