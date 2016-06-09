#!/usr/bin/ipython 
#####
#

"""
	In this code, we run the whole package analysis suite.
	
"""
from sys import argv,exit
from multiprocessing import Pool

from downloadFromMass import  downloadMass, findLastFinishedYear
from analysis_timeseries import analysis_timeseries, singleTimeSeries, timeseriesDict, singleTimeSeriesProfile
from analysis_p2p import analysis_p2p
from makeReport import html5Maker
from UKESMpython import folder


def timeseriesParrallel(index):
	print "timeseriesParrallel",index, jobID, 'START'
	key = timeseriesDict[index]
	singleTimeSeries(jobID, key,)
	print "timeseriesParrallel",index, jobID, 'SUCESS'	
	
def timeseriesParrallelProfile(index):
	print "timeseriesParrallelProfile",index, jobID, 'START'
	key = timeseriesDict[index]
	singleTimeSeriesProfile((jobID, key,)
	print "timeseriesParrallelProfile",index, jobID, 'SUCESS'



def theWholePackage(jobID):
	print "########\nThe Whole Package:\tStarting job", jobID 
#	downloadMass(jobID)

	cores = 4
		
	
	print "########\nThe Whole Package:\tStarting Time series (surface only)", jobID 	
	suite = 'all'
	
	

			  

	print "########\nThe Whole Package:\tStarting Time series ", jobID 
	remaining = sorted(timeseriesDict.keys())

   	p = Pool(cores)
    	p.map(timeseriesParrallel,remaining)
	#assert 0
	
	print "########\nThe Whole Package:\tLocating final year of model data", jobID 		        
        year = findLastFinishedYear(jobID)
	print "########\nThe Whole Package:\tFinal year of model data", jobID,"is", year


	print "########\nThe Whole Package:\tRunning point to point analysis of", jobID,"on", year
	analysis_p2p(models	= ['NEMO','MEDUSA',],
		jobID 	= jobID,
		years 	= [year,], #'2075','2076',
		modelGrid = 'eORCA1',
		annual 	= True,
		noPlots = False,
		analysisSuite=suite,)        
	
		
	print "########\nThe Whole Package:\tmaking Summary report"	
	html5Maker(jobID =jobID,
		   reportdir=folder('reports/'+jobID),
		   year = year,
		   clean=True,
		   )


	print "########\nThe Whole Package:\tStarting Time series profiles", jobID 
	remaining = sorted(timeseriesDict.keys())
   	p1 = Pool(cores)
    	p1.map(timeseriesParrallelProfile,remaining)
	#assert 0
		
	print "########\nThe Whole Package:\tmaking Summary report"	
	html5Maker(jobID =jobID,
		   reportdir=folder('reports/'+jobID),
		   year = year,
		   clean=False,
		   )



if __name__=="__main__":	
	try:	jobID = argv[1]
	except:	
		print "Please provide a job ID"
		exit()
		
	theWholePackage(jobID)
		
