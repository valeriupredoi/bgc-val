#!/usr/bin/ipython
from mergeNC import mergeNC
from deMoraTools import folder
from glob import glob
from iMarNetPython import getFileList,shouldIMakeFile,getCalendar
from os.path import basename,exists
from sys import argv
"""	The goal of this code is to have a simple way to make climatology code.
"""

# the path of least resistance is mergeNC, first for each month.
def run(jobID,key,runType):
	#jobID='xhonp'#xhono
	#key = 'HighResp'#'1894'#'clim'
	#runType	= 'Diag'#'ERSEM'#'SalTempWind'
		
	foldIn = '/data/euryale7/scratch/ledm/iMarNet/'+jobID+'/MEANS/'+jobID+'o_'

	yearKey=False
	try:	#		if key in ['1891', '1894', '1895','1899','1982','2006','1948', ]:		    
		float(key)
		yearKey=True	
		print "Key is a specific year:", key
	except:	pass
	
	if runType == 'ERSEM':
		keys = ['deptht','nav_lat','nav_lon','time_counter','Z4c','Z5c', 'Z6c', 'B1c','N1p','N3n','N4n','N5s','N7f','P1c','Chl1','P2c','Chl2','P3c','Chl3','P4c','Chl4',]
		L = 'P'
	if runType == 'ERSEMFull':
		keys = ['deptht','nav_lat','nav_lon','time_counter','Z4c','Z5c', 'Z6c',
			 'B1c','N1p','N3n','N4n','N5s','N7f','P1c','Chl1','P2c','Chl2','P3c','Chl3','P4c','Chl4','R1c','R2c','R4c','R6c','R8c','O3c','O2o',]
		L = 'P'
	if runType == 'Detritus':
		keys = ['deptht','nav_lat','nav_lon','time_counter',
			 'R1c','R2c','R4c','R6c','R8c',]
		L = 'P'	
	if runType == 'ERSEMphyto':
		keys = ['deptht','nav_lat','nav_lon','time_counter','N1p','N3n','N4n','N5s','N7f','O3c',
			'P1c','Chl1','P1s','P1n','P1p','P1f',
			'P2c','Chl2','P2n','P2p','P2f', 
			'P3c','Chl3','P3n','P3p','P3f',
			'P4c','Chl4','P4n','P4p','P4f',]
		L = 'P'
		
				
	if runType == 'SalTemp':
		keys = ['deptht','nav_lat','nav_lon','time_counter','vosaline','votemper', ]
		L = 'T'
	if runType == 'SalTempWind':
		keys = ['deptht','nav_lat','nav_lon','time_counter','vosaline','votemper', 'sowindsp','sosstsst','somxl010',]
		L = 'T'		
	if runType == 'Diag':
		keys = ['deptht','nav_lat','nav_lon','time_counter','pCO2w', 'netPP', 'fAirSeaC','chl','EIR' ]
		L = 'D'
		
		
		
		
	if jobID[:4]=='xhon':
		months = sorted(['0131', '0801','0302','0901', '0402','1001', '0502','1101', '0602','1201', '0702','1231'])

	if jobID[:4]=='xjez' and key in ['2001', ]:
		print jobID, key
		months = sorted([ '20010301', '20010501', '20010701', '20010901', '20011101', '20010130', '20010330', '20010530', '20010730', '20010930', '20011130','20020101',])	
	mergedFiles = []
	
	cal = getCalendar(jobID)
	for month in months:
		if jobID[:4]=='xhon':	
			if key =='clim':
			    filesIn = getFileList([foldIn+'199[789]'+month+'m01'+L+'.nc',foldIn+'200[0-9]'+month+'m01'+L+'.nc',])
			    print filesIn

			if key =='fullClim':		
			    filesIn = getFileList([foldIn+'19[789]?'+month+'m01'+L+'.nc',foldIn+'200[0-9]'+month+'m01'+L+'.nc',])
			    if month == '1231':
				filesIn.extend(getFileList([foldIn+'19[789]?0101m01'+L+'.nc',foldIn+'200[0-9]0101m01'+L+'.nc',]))
		
		
			if yearKey:
			
				if month not in ['1231',]:
				    filesIn = getFileList([foldIn+key+month+'m01'+L+'.nc',])
				else:
				    file0101 = foldIn+str(int(key)+1)+'0101'+'m01'+L+'.nc'
				    print "trying file0101 instead:",file0101			    
				    if exists(file0101):
				    	print "Using file0101 instead:",file0101
				 	filesIn = getFileList([file0101,])
				    else:
				    	filesIn = getFileList([foldIn+key+month+'m01'+L+'.nc',])
				    	
		if jobID[:4]=='xjez':
			print foldIn+'*'+month+'*m01'+L+'.nc'
		    	filesIn = getFileList([foldIn+'*'+month+'*m01'+L+'.nc',])

		    

		
		if jobID=='xhonp' and key in ['HighResp', ]:
		    filesIn = getFileList([foldIn+'189[3]'+month+'m01'+L+'.nc',])
		    
		print "filesIn:", filesIn
		fileOut = folder('/tmp/outNetCDF/tmp-Clims')+basename(filesIn[0])[:-3]+'_'+key+'_'+runType+'.nc'
		print fileOut
		
		mergedFiles.append(fileOut)
		if exists(fileOut):continue
		m = mergeNC( filesIn, fileOut, keys, timeAverage=True,debug=True,calendar=cal)
		del m
		
	#if key == 'clim':
	#	filenameOut = folder('outNetCDF/Climatologies')+jobID+'_'+key+'_'+runType+'.nc'
	#if key in [ '2006','2001', '1982','1948', '1894', "HighResp",]:
	
	filenameOut = folder('/data/euryale7/scratch/ledm/iMarNet_postProcessed/iMarNetp2p/outNetCDF/'+jobID+'-'+key)+jobID+'_'+key+'_'+runType+'.nc'
	
	if not shouldIMakeFile(mergedFiles,filenameOut): return		
	
	m = mergeNC( mergedFiles, filenameOut, keys, timeAverage=False,debug=True,calendar=cal)

def main():
	jobID='xhonp'#xhono'xjeza' #
	key ='clim'#'1894'# '2001'#'clim'#'fullClim'#'2006'#'clim'#'2001' #'1982'#'1948' #'HighResp'#'1894'#'clim'
	runTypes= ['Diag','ERSEM',]#'SalTempWind','ERSEMFull','ERSEMphyto','Detritus', ]#'SalTempWind', ]# ]#]#]
		
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
