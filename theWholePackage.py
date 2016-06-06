
#####
#

"""
	In this code, we run the whole package analysis suite.
	
"""
from sys import argv,exit

from downloadFromMass import  downloadMass, findLastFinishedYear
from analysis_timeseries import analysis_timeseries
from analysis_p2p import analysis_p2p
from makeReport import html5Maker

def theWholePackage(jobID):
	print "The Whole Package:\tStarting job", jobID 
	downloadMass(jobID)
	
	
	print "The Whole Package:\tStarting Time series (surface only)", jobID 	
	suite = 'debug'
	analysis_timeseries(jobID =jobID,analysisSuite=suite, z_component = 'SurfaceOnly',)

	print "The Whole Package:\tStarting Time series (Full depth)", jobID 		
        analysis_timeseries(jobID =jobID,analysisSuite=suite, z_component = 'FullDepth',)#clean=1) 	 
        
	print "The Whole Package:\tLocating final year of model data", jobID 		        
        year = findLastFinishedYear(jobID)
	print "The Whole Package:\tFinal year of model data", jobID,"is", year


	print "The Whole Package:\tRunning point to point analysis of", jobID,"on", year
	analysis_p2p(models	= ['NEMO','MEDUSA',],
		jobID 	= jobID,
		years 	= [year,], #'2075','2076',
		modelGrid = 'eORCA1',
		annual 	= True,
		noPlots = False,
		analysisSuite=suite,)        
		
	print "The Whole Package:\tmaking Summary report"	
	html5Maker(jobID =jobID,
		   reportdir=folder('reports/'+jobID,
		   year = year,
		   clean=True,
		   )
	


if __name__=="__main__":	
	try:	jobID = argv[1]
	except:	
		print "Please provide a job ID"

	theWholePackage(jobID)
		
