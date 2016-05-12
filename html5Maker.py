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

from UKESMpython import folder
from pftnames import getLongName
from glob import glob
from sys import argv
import os 

import shutil

from html5 import html5Tools

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def addImageToHtml(fn,imagesfold):
	#####
	# Note that we use three paths here.
	# fn: The original file path relative to here
	# newfn: The location of the new copy of the file relative to here
	# relfn: The location of the new copy relative to the index.html
	
	newfn = imagesfold+os.path.basename(fn)
	if not os.path.exists(newfn):
		shutil.copy2(fn, newfn)
	relfn = newfn.replace(reportdir,'./')	
	return relfn
	            
def html5Maker(
		jobID = 'u-ab749',
		reportdir = '../../html5report',
		
	):

	
	#####
	# Delete old files
	print "Removing old files from:",reportdir
	try:shutil.rmtree(reportdir)
	except: pass

	reportdir = folder(reportdir)
	
	####
	# Copy all necceasiry objects and templates to the report location:
	print "Copying html and js assets to", reportdir
	copytree('html5/html5Assets', reportdir)
	indexhtmlfn 	= reportdir+"index.html"
	try:os.rename(reportdir+'index-template.html', indexhtmlfn)
	except: pass

	imagesfold 	= folder(reportdir+'images/')
	def newImageLocation(fn):
		return imagesfold+os.path.basename(fn)
	
	#####
	#

	html5Tools.writeDescription(
				indexhtmlfn,
				'Validation of the job: '+jobID,
				)

	#####
	# Add time series regional plots:
	#key = 'ignoreInlandSeas'
	fields = ['Alkalinity', 
		  'Nitrate',
		  'Silicate', 
		  'Temperature', 
		  'Salinity', 
		  'Oxygen',
		  'DIC',
		  'Chlorophyll_cci', 
		  'IntegratedPrimaryProduction_OSU', 

		  'ExportRatio', 
		  'MLD',
		  #  'IntegratedPrimaryProduction_1x1' , 
		  #'Chlorophyll_pig' , 
		  #'AirSeaFluxCO2' , ]
		 ]
	regions = ['Global',
		  'ignoreInlandSeas',
		  'Equator10', 
		  'ArcticOcean',
		  'NorthernSubpolarAtlantic',
		  'NorthernSubpolarPacific',
		  'SouthernOcean',
		  'Remainder',]
		
	summarySection = 0
	summarySections = True	
	plotbyregion = 0#True
	plotbyfield = 0#
	plotbyfieldandregion = True
	
	
	
	if summarySection:
		files = {}
		vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*ignoreInlandSeas*.png')
		vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*/*/*ignoreInlandSeas*hist.png'))
		vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*/*/*ignoreInlandSeas*robinquad.png'))			
		vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*/*/*ignoreInlandSeas*scatter.png'))	
						
		print './images/'+jobID+'/timeseries/*/percentiles*ignoreInlandSeas*.png:',vfiles
		if len(vfiles) == 0 : assert 0 
				
		for fn in vfiles:

			relfn = addImageToHtml(fn,imagesfold)
			
			#####
			# Create custom title by removing extra bits.
			title = html5Tools.fnToTitle(relfn).split(' ')
				
			files[relfn] = ', '.join([getLongName(t) for t in title])
			print "Adding ",relfn,"to script"
		print files			
		if len(files.keys()) == 0 : assert 0 
						
		html5Tools.AddSection(indexhtmlfn,'Summary','Summary', Description='A selection of the important metrics of model performance',Files = files)


	
	if summarySections:

		SectionTitle= 'Summary'
		hrefs 	= []
		Titles	= {}
		SidebarTitles = {}
		Descriptions= {}
		FileLists	= {}
		
		region = 'ignoreInlandSeas'
		for key in sorted(fields):
			href = 	key+'-'+region
			hrefs.append(href)
			Titles[href] = 	getLongName(region) +' '+	getLongName(key)
			SidebarTitles[href] = getLongName(key)				
			Descriptions[href] = getLongName(key) +' '+	getLongName(region)
			FileLists[href] = {}
			
			#####
			# Determine the list of files:
			vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*.png')
			vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*hist.png'))
			vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*robinquad.png'))			
			vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*scatter.png'))						
		
			#####
			# Create plot headers for each file.
			for fn in vfiles:
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold)
				####
				# WOA fields that also produce transects, etc.
				if key in ['Nitrate', 'Silicate', 'Temperature', 'Salinity', 'Oxygen',] and fn.lower().find('surface')<0:continue
				
				#####
				# Create custom title by removing extra bits.
				title = html5Tools.fnToTitle(relfn).split(' ')
				#for k in ['percentiles', jobID, key,'percentiles']:
				#	try: title.remove(k)
				#	except:pass
		
	
				FileLists[href][relfn] = ', '.join([getLongName(t) for t in title])
				print "Adding ",relfn,"to script"
			
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists)
					
						
	
	
	
	if plotbyregion:
		for key in regions:
			vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*.png')
			files = {}
			for fn in vfiles:
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold)
				
				#####
				# Create custom title by removing extra bits.
				title = html5Tools.fnToTitle(relfn).split(' ')
				#title.remove('Surface').remove('surface')
				title.remove('percentiles')
				title.remove(jobID)
				title.remove(key)
				
			
				files[relfn] = ' '.join([getLongName(t) for t in title])
				print "Adding ",relfn,"to script"
				
			html5Tools.AddSection(indexhtmlfn,'ts'+key,getLongName(key)+ ' TS', Description=getLongName(key)+' Time Series plots',Files = files)

	if plotbyfield:
		for key in fields:
			vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*.png')
			files = {}
			for fn in vfiles:
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold)
				
				#####
				# Create custom title by removing extra bits.
				title = html5Tools.fnToTitle(relfn).split(' ')
				for k in ['percentiles', jobID, key,'percentiles']:
					try: title.remove(k)
					except:pass
				
			
				files[relfn] = ' '.join([getLongName(t) for t in title])
				print "Adding ",relfn,"to script"
				
			html5Tools.AddSection(indexhtmlfn,'ts'+key,getLongName(key)+ ' TS', Description=getLongName(key)+' Time Series plots',Files = files)

	if plotbyfieldandregion:
		for key in sorted(fields):
			SectionTitle= getLongName(key)
			hrefs 	= []
			Titles	= {}
			SidebarTitles = {}
			Descriptions= {}
			FileLists	= {}
			for region in regions:
				href = 	key+'-'+region
				hrefs.append(href)
				Titles[href] = 	getLongName(region) +' '+	getLongName(key)
				SidebarTitles[href] = getLongName(region)				
				Descriptions[href] = getLongName(key) +' '+	getLongName(region)
				FileLists[href] = {}
				
				#####
				# Determine the list of files:
				vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*.png')
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/2007/*'+key+'*/BGCVal/*'+region+'*'+key+'*hist.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/2007/*'+key+'*/BGCVal/*'+region+'*'+key+'*robinquad.png'))			
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/2007/*'+key+'*/BGCVal/*'+region+'*'+key+'*2007.png'))						
				
				#####
				# Create plot headers for each file.
				for fn in vfiles:
					#####
					# Copy image to image folder and return relative path.
					relfn = addImageToHtml(fn, imagesfold)
				
					#####
					# Create custom title by removing extra bits.
					title = html5Tools.fnToTitle(relfn).split(' ')
					for k in ['percentiles', jobID, key,'percentiles']:
						try: title.remove(k)
						except:pass
				
			
					FileLists[href][relfn] = ' '.join([getLongName(t) for t in title])
					print "Adding ",relfn,"to script"
				
			html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
					SidebarTitles=SidebarTitles,#
					Titles=Titles, 
					Descriptions=Descriptions,
					FileLists=FileLists)

#			html5Tools.AddSection(indexhtmlfn,key+'-'+region,longnames, Description=longnames+' plots',Files = files)




	print "-------------\nSuccess\ntest with:\nfirefox",indexhtmlfn
	print "To zip it up:\ntar cfvz  report-"+jobID+".tar.gz ",reportdir





if __name__=="__main__":	
	try:	jobID = argv[1]
	except:	
		jobID = "u-ab749"
	try: 		reportdir = argv[2]
	except: 	reportdir =folder('../../html5report')
		
	html5Maker(jobID =jobID,reportdir=reportdir)

















