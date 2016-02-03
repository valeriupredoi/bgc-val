#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
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

#Standard Python modules:
from sys import argv,exit
from os.path import exists
from calendar import month_name
from socket import gethostname

from glob import glob



import numpy as np


#Specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis#trafficlightsPlots, getMeanSurfaceChl,getHorizontalSlice




		

def analysis_timeseries(jobID = "u-ab671",
			clean = 0,
			):

	"""
		The role of this code is to produce time series analysis.
		
	"""	

	
	#doIntPP = True
	#if doIntPP:
	# code plan:
	# for each metric.

	#	for each model file:
	#	for r in regions:	
	#		iterate over a list of analysis	
	#		extract mean, depth int, etc for each region.
	#		do it monthly, then annual.
	#		Save it into a 
	#	load data
	#	for r in regions:
	#		iterate over a list of analysis
	#		produce a distribution of values.
	#		do it monthly, then annual.	
	#		

	#jobID = "u-ab671"

	

	
	name = 'Chlorophyll'
	modelcoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '365_day',}	# model doesn't need time dict.
	datacoords 	= {'t':'index_t', 'z':'DEPTH',  'lat': 'LATITUDE', 'lon': 'LONGITUDE', 'cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
	
	modeldetails = 	{'name': 'Chlorophylla', 'vars':['CHN','CHD'], 'convert': ukp.sums,'units':'mg C/m^3'}
	datadetails  =  {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
	
	layers 		= ['Surface','100m','200m',]
	regions 	= ['Global',]
	metrics 	= ['mean','median']
	modelFiles  	= sorted(glob("/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
	dataFile 		= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"+"MarEDat20121001Pigments.nc"	
	datasource 	= 'MAREDAT'
	model		= 'MEDUSA'



	shelvedir 	= ukp.folder('shelves/timeseries/'+jobID)
	imagedir	 = ukp.folder('images/timeseries/'+jobID)
	
	timeseriesAnalysis(
			modelFiles, 
			dataFile,
			dataType	= name,
  			modelcoords 	= modelcoords,
  			modeldetails 	= modeldetails,
  			datacoords 	= datacoords,
  			datadetails 	= datadetails,								
			datasource	= datasource,
			model 		= model,
			jobID		= jobID,
			layers	 	= layers,
			regions	 	= regions,			
			metrics	 	= metrics,
			workingDir	= shelvedir,
			imageDir	= imagedir,					
			#grid		= grid,
			#gridFile	= gridFile,
			clean 		= clean,
		)


if __name__=="__main__":	
	analysis_timeseries(jobID = "u-ab671")		
	analysis_timeseries(jobID = "u-ab749")			
	
	
	
	
	 
