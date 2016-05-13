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
            try:shutil.copytree(s, d, symlinks, ignore)
            except:pass
        else:
            try:shutil.copy2(s, d)
            except: pass

def addImageToHtml(fn,imagesfold):
	#####
	# Note that we use three paths here.
	# fn: The original file path relative to here
	# newfn: The location of the new copy of the file relative to here
	# relfn: The location of the new copy relative to the index.html
	
	newfn = imagesfold+os.path.basename(fn)
	if not os.path.exists(newfn):
		shutil.copy2(fn, newfn)
	else:
		####
		# Check if file is newer than the one in images.
		if os.path.getmtime(newfn) > os.path.getmtime(fn):
			shutil.remove(newfn)
			shutil.copy2(fn, newfn)			
	relfn = newfn.replace(reportdir,'./')	
	return relfn
	            
def html5Maker(
		jobID = 'u-ab749',
		reportdir = '../../html5report',
		clean = False
		
	):

	
	if clean:
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
		  'AirSeaFluxCO2' , 
		 ]
	regions = ['Global',
		  'ignoreInlandSeas',
		  'Equator10', 
		  'ArcticOcean',
		  'NorthernSubpolarAtlantic',
		  'NorthernSubpolarPacific',
		  'SouthernOcean',
		  'Remainder',
		   ]
		
	summarySections = True	
	plotbyfieldandregion = True
	
	
	
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
			Descriptions[href] = '' #getLongName(key) +' '+	getLongName(region)
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
				if key in ['Nitrate', 'Silicate', 'Temperature', 'Salinity', 'Oxygen','DIC','Alkalinity'] and fn.lower().find('surface')<0:continue
				
				#####
				# Create custom title by removing extra bits.
				title = html5Tools.fnToTitle(relfn).split(' ')
		
	
				FileLists[href][relfn] = ', '.join([getLongName(t) for t in title])
				print "Adding ",relfn,"to script"
			
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists)
					
						
	
	
	


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
				Descriptions[href] = '' #getLongName(key) +' '+	getLongName(region)
				FileLists[href] = {}
				
				#####
				# Determine the list of files:
				vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*.png')
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*hist.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*robinquad.png'))			
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*scatter.png'))							
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*Transect/*/*'+region+'*'+key+'*hov.png'))
				#####
				# Create plot headers for each file.
				for fn in vfiles:
					#####
					# Copy image to image folder and return relative path.
					relfn = addImageToHtml(fn, imagesfold)
				
					#####
					# Create custom title by removing extra bits.
					title = html5Tools.fnToTitle(relfn).split(' ')
					for k in ['percentiles', jobID, key,'percentiles',key+'vs'+key]:
						try: title.remove(k)
						except:pass
					for i,t in enumerate(title):
						if t[:len(key)] == key: title[i] = t[len(key):]
				
			
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

















