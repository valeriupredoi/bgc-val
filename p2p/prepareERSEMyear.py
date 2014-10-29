#!/usr/bin/ipython
from glob import glob
from os.path import basename,exists
from sys import argv
import numpy as np
import UKESMpython as ukp 
#from UKESMpython import folder
#from UKESMpython import getFileList,shouldIMakeFile#,getCalendar

from mergeNC import mergeNC
"""	The goal of this code is to have a simple way to make climatology code.
"""


def run(jobID,key,runType):
		
	foldIn = '/data/euryale7/scratch/ledm/iMarNet/'+jobID+'/MEANS/'+jobID+'o_'

	try:		    
		float(key)
		yearKey=True	
		print "Key is a specific year:", key
		years = [key,]
	except:	
		yearKey=False
		years = [str(y) for y in np.arange(1997,2008)]
	
	baseline = ['deptht','nav_lat','nav_lon','time_counter',]
	if runType == 'ERSEMNuts':
		keys = ['N1p','N3n','N4n','N5s','N7f',]	

		L = 'P'
		
	if runType == 'ERSEMphytoBm':
		keys = ['P1c',	'P2c','P3c','P4c',]
		L = 'P'
	if runType == 'ERSEMphytoChl':
		keys = ['Chl1','Chl2','Chl3','Chl4',]
		L = 'P'		
	if runType == 'ERSEMzoo':
		keys = ['Z4c',	'Z5c','Z6c',]
		L = 'P'		

	if runType == 'ERSEMbac':
		keys = ['B1c',]
		L = 'P'		
						
	if runType == 'ERSEMMisc':
		keys = ['pCO2w', 'netPP', 'fAirSeaC','chl','EIR' ]
		L = 'D'
		
	if runType == 'NEMO':
		keys = ['vosaline','votemper', 'somxl010']
		L = 'T'
		
	#if runType == 'NEMOSalTempWind':
	#	keys = ['deptht','nav_lat','nav_lon','time_counter','vosaline','votemper', 'sowindsp','sosstsst','somxl010',]
	#	L = 'T'	
		
		

	months = sorted(['0131', '0801','0302','0901', '0402','1001', '0502','1101', '0602','1201', '0702','1231'])
	cal = '365_day'		

	#if jobID[:4]=='xjez' and key in ['2001', ]:
	#	print jobID, key
	#	months = sorted([ '20010301', '20010501', '20010701', '20010901', '20011101', '20010130', '20010330', '20010530', '20010730', '20010930', '20011130','20020101',])	
	#	cal = '360_day'
		
	mergedFiles = []
	

	for month in months:
		#filesIn.extend(getFileList([foldIn+'19[789]?0101m01'+L+'.nc',foldIn+'200[0-9]0101m01'+L+'.nc',]))
		
		for ykey in years:
			if month not in ['1231',]:
			    filesIn = ukp.getFileList([foldIn+ykey+month+'m01'+L+'.nc',])
			else:
			    file0101 = foldIn+str(int(ykey)+1)+'0101'+'m01'+L+'.nc'
			    print "trying file0101 instead:",file0101			    
			    if exists(file0101):
			    	print "Using file0101 instead:",file0101
			 	filesIn = ukp.getFileList([file0101,])
			    else:
			    	filesIn = ukp.getFileList([foldIn+ykey+month+'m01'+L+'.nc',])
			    	
				    
		#if jobID=='xhonp' and key in ['HighResp', ]:
		 #   filesIn = ukp.getFileList([foldIn+'189[3]'+month+'m01'+L+'.nc',])
		    
		print "filesIn:", filesIn
		fileOut = ukp.folder('/tmp/outNetCDF/tmp-Clims')+basename(filesIn[0])[:-3]+'_'+key+'_'+runType+'.nc'
		print fileOut
		
		mergedFiles.append(fileOut)
		if exists(fileOut):continue
		m = mergeNC( filesIn, fileOut, keys, timeAverage=True,debug=True,calendar=cal)
		del m
		
	#if key == 'clim':
	#	filenameOut = folder('outNetCDF/Climatologies')+jobID+'_'+key+'_'+runType+'.nc'
	#if key in [ '2006','2001', '1982','1948', '1894', "HighResp",]:
	
	filenameOut = ukp.folder('/data/euryale7/scratch/ledm/UKESM/ERSEM/'+jobID+'/'+key)+jobID+'_'+key+'_'+runType+'.nc'
	
	if not ukp.shouldIMakeFile(mergedFiles,filenameOut): return		
	
	m = mergeNC( mergedFiles, filenameOut, keys, timeAverage=False,debug=True,calendar=cal)

def main():
	jobID='xhonp'#xhono'xjeza' #
	key ='clim'#'1894'# '2001'#'clim'#'fullClim'#'2006'#'clim'#'2001' #'1982'#'1948' #'HighResp'#'1894'#'clim'
	runTypes= ['NEMO','ERSEMNuts','ERSEMphytoBm','ERSEMphytoChl','ERSEMzoo', 'ERSEMMisc','ERSEMbac']
	#'SalTempWind','ERSEMFull','ERSEMphyto','Detritus', ]#'SalTempWind', ]# ]#]#]
		
	try: 	
		jobID = argv[1]
		key   = argv[2]
		print "Using command line arguments:", jobID,key
	except:
		jobID = 'xhonp'
		key = 'clim'
		print "Not using command line arguments"
		jobs = ['xhonp','xhont','xhonu','xhonv','xjezd', ]
		for j in jobs:
		  for r in runTypes:
		  	run(j,'1899',r)
		return
		


	
	for r in runTypes: run(jobID,key,r)
main()	
