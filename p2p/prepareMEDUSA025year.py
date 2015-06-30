#!/usr/bin/ipython
#
# Copyright 2014, Plymouth Marine Laboratory
#
# This file is part of the ukesm-validation library.
#
# ukesm-validation is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license. 

# ukesm-validation is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with ukesm-validation.
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


from glob import glob
from os.path import basename,exists
from sys import argv
import numpy as np
import UKESMpython as ukp 
#from UKESMpython import folder
#from UKESMpython import getFileList,shouldIMakeFile#,getCalendar
from netCDF4 import Dataset
from changeNC import changeNC,AutoVivification
from mergeNC import mergeNC
from pruneNC import pruneNC
"""	The goal of this code is to have a simple way to make climatology code.
"""


def run(jobID,key,runType,foldIn):

	try:		    
		float(key)
		yearKey=True	
		print "Key is a specific year:", key
		years = [key,]
	except:	
		assert false
		#yearKey=False
		#years = [str(y) for y in np.arange(1997,2008)]
	
	baseline = ['deptht','nav_lat','nav_lon','time_counter',]

	if runType == 'CH':
		keys = ['CHD','CHN',]
		finalKeys = ['CHL',]	
		L = '_ptrc_T'

	if runType == 'DIN':
		keys = ['DIN',]
		finalKeys = ['DIN',]	
		L = '_ptrc_T'
				
	#months = sorted(['0121', '0821','0321','0921', '0421','1021', '0521','1121', '0621','1221', '0721','1221'])
	cal = '365_day'		

	#if jobID[:4]=='xjez' and key in ['2001', ]:
	#	print jobID, key
	#	months = sorted([ '20010301', '20010501', '20010701', '20010901', '20011101', '20010130', '20010330', '20010530', '20010730', '20010930', '20011130','20020101',])	
	#	cal = '360_day'
		
	mergedFiles = []
	

	filesIn = sorted(glob(foldIn+'/xjwkio_1m_1979*'+L+'*.nc'))
		
	print "filesIn:"
		
	for fn in filesIn:
	
		prunedfn = ukp.folder('/tmp/outNetCDF/tmp-Clims')+basename(fn)[:-3]+'_'+key+'_'+runType+'.nc'
		print fn, '--->', prunedfn
		

		if not exists(prunedfn):
			m = pruneNC( fn, prunedfn, keys, debug=True)#,calendar=cal)
		if runType == 'CH':
			nc = Dataset(prunedfn,'r')
			
			fileOut = prunedfn.replace('.nc','_chl.nc')
			if not exists(fileOut):
				av = AutoVivification()
				av['CHN']['name']='False'
				av['CHD']['name']	='CHL'
				av['CHD']['units']	='[mg Chl/m3]'
				av['CHD']['long_name']	='Total Chlorophyll'
				av['CHD']['newDims']	=(u'time_counter', u'deptht', u'y', u'x') 
				av['CHD']['newData'] = nc.variables['CHD'][:] + nc.variables['CHN'][:]
				nc.close()
				c = changeNC(prunedfn,fileOut,av,debug=True)
			print fileOut
			#print prunedfn
		if runType in ['DIN',]:fileOut = prunedfn
		
		mergedFiles.append(fileOut)		
		#del m
		

	
	filenameOut = ukp.folder('/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/'+jobID+'_postProc/'+key)+jobID+'_'+key+'_'+runType+'.nc'
	
	if not ukp.shouldIMakeFile(mergedFiles,filenameOut): return		
	
	m = mergeNC( mergedFiles, filenameOut, finalKeys, timeAverage=False,debug=True,calendar=cal)

def main():
	jobID='xhonp'#xhono'xjeza' #
	key ='clim'#'1894'# '2001'#'clim'#'fullClim'#'2006'#'clim'#'2001' #'1982'#'1948' #'HighResp'#'1894'#'clim'
	runTypes= ['DIN',]#'CH','
	#'ERSEMNuts','ERSEMphytoBm','ERSEMphytoChl','ERSEMzoo', 'ERSEMMisc','ERSEMbac']
	#'SalTempWind','ERSEMFull','ERSEMphyto','Detritus', ]#'SalTempWind', ]# ]#]#]
	

			
	try: 	
		jobID = argv[1]
		key   = argv[2]
		print "Using command line arguments:", jobID,key
	except:
		jobID = 'xjwki'
		key = '1979'
		print "Not using command line arguments"
		jobs = ['xjwki', ]
		for j in jobs:
		  for r in runTypes:
			#foldIn = '/data/euryale7/scratch/ledm/iMarNet/'+j+'/MEANS/'		  
			foldIn = '/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/'+j		  			
		  	run(j,key,r,foldIn)
		return
		


	
	for r in runTypes: 
		#foldIn = '/data/euryale7/scratch/ledm/iMarNet/'+jobID+'/MEANS/'	
		foldIn = '/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/'+j		  					
		run(jobID,key,r,foldIn)
main()	
