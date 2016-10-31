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
.. module:: analysis_level3_omz
   :platform: Unix
   :synopsis: A script to produce a level 3 analysis for Oxygen Minimum Zones.

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

#####	
# Load specific local code:
import UKESMpython as ukp
from timeseries import timeseriesAnalysis
from timeseries import profileAnalysis
from timeseries import timeseriesPlots as tsp 

#####
# User defined set of paths pointing towards the datasets.
import paths	


def analysis_omz(jobID=''):	
	annual = True
	
	analysisKeys = []
        analysisKeys.append('O2')                       # WOA Oxygen	
	analysisKeys.append('OMZMeanDepth')		# OMZ mean depth
	analysisKeys.append('OMZThickness')             # Oxygen Minimum Zone Thickness					
	analysisKeys.append('TotalOMZVolume')		# Total OMZ volume
			
				
	analysisDict = {}
	imagedir	= ukp.folder(paths.imagedir +'/'+jobID+'/Level3/OMZ')
	shelvedir 	= ukp.folder(paths.shelvedir+'/'+jobID+'/Level3/OMZ')

	#####
	# make a link to the time series
	for a in analysisKeys:
		level1shelveFold =  paths.shelvedir+'/timeseries/'+jobID
		files = glob(level1shelveFold+'/*'+a+'*')
		for f in files:
			lnfile = shelvedir+os.basename(f)
			if os.path.exists(lnfile):continue
			print "linking ",f,lnfile
			assert 0
			os.symlink(f,lnfile)

	
	regionList	= ['Global', 'ignoreInlandSeas',
	  		'SouthernOcean','ArcticOcean',
			'Equator10', 'Remainder',
			'NorthernSubpolarAtlantic','NorthernSubpolarPacific',
			]
  	metricList = ['mean',]			
	dataD = {}		
	modeldataD = {}



	def listModelDataFiles(jobID, filekey, datafolder, annual):
		if annual:
			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1y_*_"+filekey+".nc"))
		else:
			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1m_*_"+filekey+".nc"))
					        			

	
	medusaCoords 	= {'t':'time_counter', 'z':'deptht', 'lat': 'nav_lat',  'lon': 'nav_lon',   'cal': '360_day',}	# model doesn't need time dict.
	
	
	if 'O2' in analysisKeys:
		name = 'Oxygen'
		if annual:
			av[name]['modelFiles']  = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref, annual)		
			av[name]['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
				
		av[name]['modelcoords'] 	= medusaCoords 	
		av[name]['datacoords'] 		= woaCoords
	
		av[name]['modeldetails'] 	= {'name': name, 'vars':['OXY',], 'convert': ukp.NoChange,'units':'mmol O2/m^3'}	
		av[name]['datadetails']  	= {'name': name, 'vars':['o_an',], 'convert': ukp.oxconvert,'units':'mmol O2/m^3'}

		av[name]['layers'] 		=  layerList
		av[name]['regions'] 		= regionList
		av[name]['metrics']		= metricList

		av[name]['datasource'] 		= 'WOA'
		av[name]['model']		= 'MEDUSA'

		av[name]['modelgrid']		= 'eORCA1'
		av[name]['gridFile']		= paths.orcaGridfn
		av[name]['Dimensions']		= 3

	if 'OMZMeanDepth' in analysisKeys:
		if annual:
			av['OMZMeanDepth']['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av['OMZMeanDepth']['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
		else:
			print "OMZ Thickness not implemented for monthly data"
			assert 0
			
		nc = Dataset(paths.orcaGridfn,'r')
		depths   	= np.abs(nc.variables['gdepw' ][:])
		tmask 		= nc.variables['tmask'][:]						
		nc.close()			

		omzthreshold = 20.
				
		def modelMeanOMZdepth(nc,keys):
			o2 = nc.variables[keys[0]][:].squeeze()
			meandepth = np.ma.masked_where((o2>omzthreshold)+o2.mask + (tmask==0),depths).mean(0)
			if meandepth.max() in [0.,0]: return np.array([0.,])						
			return np.ma.masked_where(meandepth==0., meandepth)

		def woaMeanOMZdepth(nc,keys):
			o2 = nc.variables[keys[0]][:].squeeze() *44.661
			pdepths = np.zeros_like(o2) 
			lons = nc.variables['lon'][:]
			lats = nc.variables['lat'][:]			
			wdepths = np.abs(nc.variables['depth'][:])
			
			for y,lat in enumerate(lats):
			    for x,lon in enumerate(lons):			  
				pdepths[:,y,x] = wdepths
			wmeanDepth = np.ma.masked_where((o2>omzthreshold)+o2.mask,pdepths).mean(0).data
			print "woaMeanOMZdepth",wmeanDepth.min(),wmeanDepth.mean(),wmeanDepth.max()
			#assert 0
			
			if wmeanDepth.max() in [0.,0]: return np.array([1000.,])
			return np.ma.masked_where(wmeanDepth==0., wmeanDepth)
				
		av['OMZMeanDepth']['modelcoords'] 	= medusaCoords 	
		av['OMZMeanDepth']['datacoords'] 	= woaCoords
	
		av['OMZMeanDepth']['modeldetails'] 	= {'name': 'OMZMeanDepth', 'vars':['OXY',],  'convert': modelMeanOMZdepth,'units':'m'}
		av['OMZMeanDepth']['datadetails']  	= {'name': 'OMZMeanDepth', 'vars':['o_an',], 'convert': woaMeanOMZdepth,'units':'m'}		
	
		av['OMZMeanDepth']['layers'] 		= ['layerless',] 
		av['OMZMeanDepth']['regions'] 		= regionList
		av['OMZMeanDepth']['metrics']		= metricList

		av['OMZMeanDepth']['datasource'] 	= 'WOA'
		av['OMZMeanDepth']['model']		= 'MEDUSA'

		av['OMZMeanDepth']['modelgrid']		= 'eORCA1'
		av['OMZMeanDepth']['gridFile']		= paths.orcaGridfn	
		av['OMZMeanDepth']['Dimensions']	= 2		
		
		
		
		
 
	if 'OMZThickness' in analysisKeys or 'OMZThickness50' in analysisKeys:
		if 'OMZThickness' in analysisKeys and 'OMZThickness50' in analysisKeys:
			print "Only run one of these at a time"
			assert 0
			
		
		if annual:
			av['OMZThickness']['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av['OMZThickness']['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
		else:
			print "OMZ Thickness not implemented for monthly data"
			assert 0
			
		nc = Dataset(paths.orcaGridfn,'r')
		thickness   	= nc.variables['e3t' ][:]
		tmask 		= nc.variables['tmask'][:]								
		nc.close()			

		if 'OMZThickness' in analysisKeys: 	omzthreshold = 20.
		if 'OMZThickness50' in analysisKeys: 	omzthreshold = 50.
				
		def modelOMZthickness(nc,keys):
			o2 = nc.variables[keys[0]][:].squeeze()
			totalthick = np.ma.masked_where((o2>omzthreshold)+o2.mask+ (tmask==0),thickness).sum(0).data
			if totalthick.max() in [0.,0]: return np.array([0.,])
			
			return np.ma.masked_where(totalthick==0., totalthick)
			#return np.ma.masked_where((arr>omzthreshold) + (arr <0.) + arr.mask,thickness).sum(0)

			
		def woaOMZthickness(nc,keys):
			o2 = nc.variables[keys[0]][:].squeeze() *44.661
			pthick = np.zeros_like(o2) 
			lons = nc.variables['lon'][:]
			lats = nc.variables['lat'][:]			
			zthick  = np.abs(nc.variables['depth_bnds'][:,0] - nc.variables['depth_bnds'][:,1])

			for y,lat in enumerate(lats):
			    for x,lon in enumerate(lons):			  
				pthick[:,y,x] = zthick
			totalthick = np.ma.masked_where((o2>omzthreshold)+o2.mask,pthick).sum(0).data
			
			if totalthick.max() in [0.,0]: return np.array([0.,])
			return np.ma.masked_where(totalthick==0., totalthick)

		av['OMZThickness']['modelcoords'] 	= medusaCoords 	
		av['OMZThickness']['datacoords'] 		= woaCoords
	
		av['OMZThickness']['modeldetails'] 	= {'name': 'OMZThickness', 'vars':['OXY',], 'convert': modelOMZthickness,'units':'m'}
		av['OMZThickness']['datadetails']  	= {'name': 'OMZThickness', 'vars':['o_an',], 'convert': woaOMZthickness,'units':'m'}
	
		av['OMZThickness']['layers'] 		= ['layerless',] 
		av['OMZThickness']['regions'] 		= regionList
		av['OMZThickness']['metrics']		= metricList

		av['OMZThickness']['datasource'] 		= 'WOA'
		av['OMZThickness']['model']		= 'MEDUSA'

		av['OMZThickness']['modelgrid']		= 'eORCA1'
		av['OMZThickness']['gridFile']		= paths.orcaGridfn	
		av['OMZThickness']['Dimensions']		= 2		
		
		
		

	if 'TotalOMZVolume' in analysisKeys or 'TotalOMZVolume50' in analysisKeys:
		if 'TotalOMZVolume' in analysisKeys and 'TotalOMZVolume50' in analysisKeys:
			print "Only run one of these at a time"
			assert 0
				
		if annual:
			av['TotalOMZVolume']['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_ptrc_T.nc"))
			av['TotalOMZVolume']['dataFile'] 		=  WOAFolder+'woa13_all_o00_01.nc'
		else:
			print "OMZ volume not implemented for monthly data"
			assert 0
			
		nc = Dataset(paths.orcaGridfn,'r')
		try:	
			pvol   = nc.variables['pvol' ][:]
			tmask = nc.variables['tmask'][:]
		except:
			tmask = nc.variables['tmask'][:]			
			area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
			pvol = nc.variables['e3t'][:] *area			
			pvol = np.ma.masked_where(tmask==0,pvol)
		nc.close()			

		if 'TotalOMZVolume' in analysisKeys:	omzthreshold = 20.
		if 'TotalOMZVolume50' in analysisKeys:	omzthreshold = 50.	
			
		def modelTotalOMZvol(nc,keys):
			arr = np.ma.array(nc.variables[keys[0]][:].squeeze())
			return np.ma.masked_where((arr>omzthreshold) + pvol.mask + arr.mask,pvol).sum()
	
			
		def woaTotalOMZvol(nc,keys):
			arr = nc.variables[keys[0]][:].squeeze() *44.661
			#area = np.zeros_like(arr[0])
			pvol = np.zeros_like(arr) 
			#np.ma.masked_wjhere(arr.mask + (arr <0.)+(arr >1E10),np.zeros_like(arr))
			lons = nc.variables['lon'][:]
			lats = nc.variables['lat'][:]			
			#lonbnds = nc.variables['lon_bnds'][:]
			latbnds = nc.variables['lat_bnds'][:]
			zthick  = np.abs(nc.variables['depth_bnds'][:,0] - nc.variables['depth_bnds'][:,1])
			
			for y,lat in enumerate(lats):
				area = ukp.Area([latbnds[y,0],0.],[latbnds[y,1],1.])
				for z,thick in enumerate(zthick):
					pvol[z,y,:] = np.ones_like(lons)*area*thick
					
			return np.ma.masked_where(arr.mask + (arr >omzthreshold)+(arr <0.),pvol).sum()
				
		av['TotalOMZVolume']['modelcoords'] 	= medusaCoords 	
		av['TotalOMZVolume']['datacoords'] 		= woaCoords
	
		av['TotalOMZVolume']['modeldetails'] 	= {'name': 'TotalOMZVolume', 'vars':['OXY',], 'convert': modelTotalOMZvol,'units':'m^3'}
		av['TotalOMZVolume']['datadetails']  	= {'name': 'TotalOMZVolume', 'vars':['o_an',], 'convert': woaTotalOMZvol,'units':'m^3'}
	
		av['TotalOMZVolume']['layers'] 		= ['layerless',] 
		av['TotalOMZVolume']['regions'] 	= ['regionless',]
		av['TotalOMZVolume']['metrics']		= ['metricless', ]

		av['TotalOMZVolume']['datasource'] 		= 'WOA'
		av['TotalOMZVolume']['model']		= 'MEDUSA'

		av['TotalOMZVolume']['modelgrid']		= 'eORCA1'
		av['TotalOMZVolume']['gridFile']		= paths.orcaGridfn	
		av['TotalOMZVolume']['Dimensions']		= 1	
	
	
  	#####
  	# Calling timeseriesAnalysis
	# This is where the above settings is passed to timeseriesAnalysis, for the actual work to begin.
	# We loop over all fiels in the first layer dictionary in the autovificiation, av.
	#	
	# Once the timeseriesAnalysis has completed, we save all the output shelves in a dictionairy.
	# At the moment, this dictioary is not used, but we could for instance open the shelve to highlight specific data,
	#	(ie, andy asked to produce a table showing the final year of data.
	
	shelves = {}
	shelves_insitu={}
	for name in av.keys():
		print "------------------------------------------------------------------"	
		print "analysis-Timeseries.py:\tBeginning to call timeseriesAnalysis for ", name

		if len(av[name]['modelFiles']) == 0:
			print "analysis-Timeseries.py:\tWARNING:\tmodel files are not found:",name,av[name]['modelFiles']
			if strictFileCheck: assert 0

		modelfilesexists = [os.path.exists(f) for f in av[name]['modelFiles']]
		if False in modelfilesexists:
			print "analysis-Timeseries.py:\tWARNING:\tnot model files do not all exist:",av[name]['modelFiles'] 
			for f in av[name]['modelFiles']: 
				if os.path.exists(f):continue
				print f, 'does not exist'
			if strictFileCheck: assert 0
			
			
		if av[name]['dataFile']!='':
		   if not os.path.exists(av[name]['dataFile']):
			print "analysis-Timeseries.py:\tWARNING:\tdata file is not found:",av[name]['dataFile']
			if strictFileCheck: assert 0

	
                #####
                # time series and traffic lights.
                tsa = timeseriesAnalysis(
                        av[name]['modelFiles'],
                        av[name]['dataFile'],
                        dataType        = name,
                        modelcoords     = av[name]['modelcoords'],
                        modeldetails    = av[name]['modeldetails'],
                        datacoords      = av[name]['datacoords'],
                        datadetails     = av[name]['datadetails'],
                        datasource      = av[name]['datasource'],
                        model           = av[name]['model'],
                        jobID           = jobID,
                        layers          = av[name]['layers'],
                        regions         = av[name]['regions'],
                        metrics         = av[name]['metrics'],
                        workingDir      = shelvedir,
                        imageDir        = imagedir,
                        grid            = av[name]['modelgrid'],
                        gridFile        = av[name]['gridFile'],
                        clean           = clean,
                )

		
		#####
		# Profile plots 	
		if av[name]['Dimensions'] == 3:
			profa = profileAnalysis(
				av[name]['modelFiles'], 
				av[name]['dataFile'],
				dataType	= name,
	  			modelcoords 	= av[name]['modelcoords'],
	  			modeldetails 	= av[name]['modeldetails'],
	  			datacoords 	= av[name]['datacoords'],
	  			datadetails 	= av[name]['datadetails'],								
				datasource	= av[name]['datasource'],
				model 		= av[name]['model'],
				jobID		= jobID,
				layers	 	= list(np.arange(102)),	# 102 because that is the number of layers in WOA Oxygen
				regions	 	= av[name]['regions'],			
				metrics	 	= ['mean',],
				workingDir	= shelvedir,
				imageDir	= imagedir,					
				grid		= av[name]['modelgrid'],
				gridFile	= av[name]['gridFile'],
				clean 		= clean,
			)
			#shelves[name] = profa.shelvefn
			#shelves_insitu[name] = profa.shelvefn_insitu
		
		#shelves[name] = tsa.shelvefn
		#shelves_insitu[name] = tsa.shelvefn_insitu



def main():
	try:	jobID = argv[1]
	except:	
		jobID = "u-ab749"

	if 'debug' in argv[1:]:	suite = 'debug'
		
	
	analysis_omz(jobID =jobID, )#clean=1)			
	#if suite == 'all':
	#analysis_timeseries(jobID =jobID,analysisSuite='FullDepth', z_component = 'FullDepth',)#clean=1)  

if __name__=="__main__":
	main()	
                    

		
		
