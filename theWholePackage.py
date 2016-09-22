#!/usr/bin/ipython 
#####
#

"""
	In this code, we run the whole package analysis suite.
	
"""


import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from sys import argv,exit
from multiprocessing import Pool

from downloadFromMass import  downloadMass, findLastFinishedYear
from analysis_timeseries import analysis_timeseries, singleTimeSeries, singleTimeSeriesProfile
from analysis_timeseries import level1KeysDict, timeseriesDict, physKeysDict
from analysis_p2p import analysis_p2p, p2pDict_level2, p2pDict_physics,single_p2p
from makeReport import html5Maker
from UKESMpython import folder


def timeseriesParrallel(index):
	print "timeseriesParrallel",index, jobID, 'START'
	key = timeseriesDict[index]
	singleTimeSeries(jobID, key,)
	print "timeseriesParrallel",index, jobID, 'SUCESS',key	
	
def timeseriesParrallelL1(index):
	print "timeseriesParrallel",index, jobID, 'START'
	key = level1KeysDict[index]
	singleTimeSeries(jobID, key,)
	print "timeseriesParrallel",index, jobID, 'SUCESS',key	

def timeseriesParrallelPhys(index):
	key = physKeysDict[index]
	print "timeseriesParrallel",index, jobID, 'START',key,index
	try:singleTimeSeries(jobID, key,)
	except:
		print "timeseriesParrallel failed for",index, jobID, key
		assert 0
	print "timeseriesParrallel",index, jobID, 'SUCESS',key
	
def p2pParrallel(index):
	print "p2pParrallel",index, jobID, 'START'
	key = p2pDict_level2[index]
	single_p2p(jobID, key, year)
	print "p2pParrallel",index, jobID, 'SUCESS',key
	
def p2pParrallel_phys(index):
	print "p2pParrallel_phys",index, jobID, 'START'
	key = p2pDict_physics[index]
	single_p2p(jobID, key, year)
	print "p2pParrallel_phys",index, jobID, 'SUCESS',key


def theWholePackage(jobID,year=False,suite = 'level1'):
        if year == False: year = '*'
	print "########\nThe Whole Package:\tStarting job", jobID , year
#	downloadMass(jobID)

	parrallel = True
	cores = 8
	#suite = 'all'
	

        print "########\nThe Whole Package:\tmaking Summary report"
        if year == False: year = '*'
        if suite =='physics':	physicsOnly=True
        else: 			physicsOnly=False
        html5Maker(jobID =jobID,
                   reportdir=folder('reports/'+jobID),
                   year = year,
                   clean=True,
                   physicsOnly=physicsOnly
                   )
#	return			  

	print "########\nThe Whole Package:\tStarting Time series (surface only)", jobID 	
	if parrallel:
		if suite =='all':	remaining = sorted(timeseriesDict.keys())[:]
		if suite =='level1':	remaining = sorted(level1KeysDict.keys())[:]
		if suite =='physics':	remaining = sorted(physKeysDict.keys())[:]		
			
	   	p = Pool(cores)
	    	if suite =='all':	p.map(timeseriesParrallel,  remaining)
	    	if suite =='level1':	p.map(timeseriesParrallelL1,remaining)
	    	if suite =='physics':	p.map(timeseriesParrallelPhys,remaining)	    	
	    	p.close()
	else:	
		analysis_timeseries(jobID =jobID,analysisSuite=suite, )#z_component = 'SurfaceOnly',)
		
		
	if year not in ['*', False]:
		if suite =='physics':pass
		else:	suite = 'level2'
		
		print "########\nThe Whole Package:\tRunning point to point analysis of", jobID,"on", year, 'in suite:',suite
		
		if parrallel:
		   	p1 = Pool(cores)
		   	if suite == 'physics':	
				remaining = sorted(p2pDict_physics.keys())[:]		   	
		   		p1.map(p2pParrallel_phys,remaining)
		    	else:	
				remaining = sorted(p2pDict_level2.keys())[:]
		    		p1.map(p2pParrallel,remaining)	
			p1.close()
		
			#####
			# And once over to make the summary target diagrams.
			analysis_p2p(models	= ['NEMO','MEDUSA',],
				jobID 	= jobID,
				years 	= [year,], #'2075','2076',
				modelGrid = 'eORCA1',
				annual 	= True,
				noPlots = True,
				analysisSuite=suite,) 
					    	
		else:	
			analysis_p2p(models	= ['NEMO','MEDUSA',],
				jobID 	= jobID,
				years 	= [year,], #'2075','2076',
				modelGrid = 'eORCA1',
				annual 	= True,
				noPlots = False,
				analysisSuite=suite,)        

	else:
		print "########\nThe Whole Package:\tNot Running point to point analysis of", jobID," because year is:", year
	






if __name__=="__main__":	

	try:	jobID = argv[1]
	except:	
		print "Please provide a job ID"
		exit()
	if 'ReportOnly' in argv[:]:ReportOnly=True
	else:	ReportOnly = False

	if 'physics' in argv[:]:physicsOnly=True
	else:	physicsOnly = False
	
        year = findLastFinishedYear(jobID,dividby=25)		
	if not ReportOnly:
		if physicsOnly:	theWholePackage(jobID,year=year,suite='physics')
		else:		theWholePackage(jobID,year=year)
		
        if year == False: year = '*'
	print "########\nThe Whole Package:\tmaking Summary report", jobID,year	
	html5Maker(jobID =jobID,
		   reportdir=folder('reports/'+jobID),
		   year = year,
		   clean=True,
		   physicsOnly=physicsOnly,
		   )
		   
		   
		
