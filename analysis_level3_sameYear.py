#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
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
"""
.. module:: analysis_level3_sameYear
   :platform: Unix
   :synopsis: A script to produce a level 3 analysis comparing the same year in two jobs.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""


#####	
# Load Standard Python modules:
from sys import argv,exit
from os.path import exists
from calendar import month_name
from socket import gethostname
from netCDF4 import Dataset
from glob import glob
from scipy.interpolate import interp1d
import numpy as np
import os,sys
from getpass import getuser
from matplotlib import pyplot
from shelve import open as shOpen
import re

#####	
# Load specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis
from timeseries import profileAnalysis
from timeseries import timeseriesPlots as tsp 

#####
# User defined set of paths pointing towards the datasets.
import paths	

def listModelDataFiles(jobID, filekey, datafolder, annual):
	if annual:
		return sorted(glob(datafolder+jobID+"/"+jobID+"o_1y_*_"+filekey+".nc"))
	else:
		return sorted(glob(datafolder+jobID+"/"+jobID+"o_1m_*_"+filekey+".nc"))
		
def overlapyears(files1,files2):
	""" 
		Finds out which years overlap.
	"""	
	filepairs = {}	
	for f1 in sorted(files1):
	    f1b= os.path.basename(f1)
	    f1y = re.findall(r'\d+', f1b)[2]
	    for f2 in files2:
		f2b =  os.path.basename(f2)
	 	f2y = re.findall(r'\d+', f2b)[2]
	 	
	 	if f2y == f1y: 	filepairs[f2y] = [f1,f2]
	return filepairs

def maskAndCompress(la,lo,d1,d2):
	m  = np.ma.array(d1).mask
	m += np.ma.array(d2).mask	
	
	return np.ma.masked_where(m,la).compressed(), np.ma.masked_where(m,lo).compressed(), np.ma.masked_where(m,d1).compressed(), np.ma.masked_where(m,d2).compressed()
			

def analysis_sy(jobID1 = 'u-af983',jobID2 = 'u-ah531', ):
	annual = True
	
	analysisKeys = []
	analysisKeys.append('SST')		
	analysisKeys.append('SSS')					
        analysisKeys.append('Ice')
	
	analysisDict = {}
	imagedir	= ukp.folder(paths.imagedir +'/'+jobID1+'-'+jobID2+'/Level3/')
	#shelvedir 	= ukp.folder(paths.shelvedir+'/'+jobID+'/Level3/'+jobID1+'-'+jobID2)

	dataD = {}		
	modeldataD = {}


	files1 = listModelDataFiles(jobID1, 'grid_T', paths.ModelFolder_pref, annual)
	files2 = listModelDataFiles(jobID2, 'grid_T', paths.ModelFolder_pref, annual)	
	
	filepairs = overlapyears(files1,files2)

	plotDetails = {}
	plotDetails['SST'] = {'name':'SST', 'key':'votemper', 'ndim':4,'longname':'Sea Surface Temperature'}
	plotDetails['SSS'] = {'name':'SSS', 'key':'vosaline', 'ndim':4,'longname':'Sea Surface Salinity'}
        plotDetails['Ice'] = {'name':'Ice', 'key':'soicecov', 'ndim':4, 'longname':'Ice fraction'}

	
	for ystr, [fp1,fp2] in filepairs.items():
		print ystr,[fp1,fp2]
		nc1 = Dataset(fp1, 'r')
		nc2 = Dataset(fp2, 'r')
	
		lons_cc = nc1.variables['nav_lon'][:]
		lats_cc = nc1.variables['nav_lat'][:]		
	
	
		for n in analysisKeys:

			filename = imagedir+plotDetails[n]['name']+'_'+ystr+'.png'
			if plotDetails[n]['ndim']==4:
				data1 = nc1.variables[plotDetails[n]['key']][0,0]
				data2 = nc2.variables[plotDetails[n]['key']][0,0]
                        if plotDetails[n]['ndim']==3:
                                data1 = nc1.variables[plotDetails[n]['key']][0]
                                data2 = nc2.variables[plotDetails[n]['key']][0]
			

			lons, lats, data1,data2 = maskAndCompress(lons_cc,lats_cc,data1,data2)
			
			ukp.robinPlotQuad(lons, lats, data1,data2,
					filename,
					titles=[jobID1,jobID2], 
					title=plotDetails[n]['longname']+' ' + ystr[:4]+'-'+ystr[5:7]+'-'+ystr[7:],
					vmin='',vmax='',)#maptype='Basemap')
	
		

def main():
	try:	
		jobID1 = argv[1]
		jobID2 = argv[2]
	except:	
		jobID = "u-ad371"
		jobID = "u-ad371"
		
	if 'debug' in argv[1:]:	suite = 'debug'
	else:	suite = 'normal'	
		
	
	analysis_sy(jobID1 =jobID1,jobID2 =jobID2, )

if __name__=="__main__":
	main()	
                    


