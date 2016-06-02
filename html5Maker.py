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
		print "cp ",newfn,fn	
		shutil.copy2(fn, newfn)
	else:
		####
		# Check if file is newer than the one in images.
		if os.path.getmtime(newfn) > os.path.getmtime(fn):
			print "removing old file",fn
			shutil.remove(newfn)
			shutil.copy2(fn, newfn)			
			print "cp ",newfn,fn
			
	relfn = newfn.replace(reportdir,'./')	
	return relfn
	            
def html5Maker(
		jobID = 'u-ab749',
		reportdir = '../../html5report',
		year = '*',
		clean = False,
		doZip= True,
	):

	
	if clean:
		#####
		# Delete old files
		print "Removing old files from:",reportdir
		try:shutil.rmtree(reportdir)
		except: pass

	reportdir = folder(reportdir)
	year = str(year)
	
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
	
	descriptionText = 'Validation of the job: '+jobID
	if year != '*':	descriptionText+=', in the year: ' +year
	html5Tools.writeDescription(
				indexhtmlfn,
				descriptionText,
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
		  'TotalIntegratedPrimaryProduction',
		  'ExportRatio', 
		  'LocalExportRatio', 		  
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

	#####
	# Two switches to turn on Summary section, and groups of plots of field and region.		
	summarySections = True	
	plotbyfieldandregion = True
	
	
	
	
	#####
	# A list of caveats linked to specific datasets or regions, or jobs.
	ListofCaveats, ListofCaveats_regions = {}, {}
	ListofCaveats['ExportRatio'] = 'Note that there is no historic data set for this variable.'
	ListofCaveats['MLD'] = 'Note that the Model MLD is calculated with based on a sigma_0 difference of 0.01 with the surface where as data uses as \
			sigma_0 difference of +/- 0.2 degrees from a depth on 10m.'
			
	if jobID == 'u-ad371':
		ListofCaveats['Chlorophyll_cci']= 'Note that the Non-diatom chlorophyll failed in run:'+jobID
		ListofCaveats['IntegratedPrimaryProduction_OSU']= 'Note that the Non-diatom chlorophyll does not contribute to IntPP in this run:'+jobID
		
		
	if summarySections:
		sumfields = [
			  'TotalIntegratedPrimaryProduction',
			  'IntegratedPrimaryProduction_OSU', 
			  'AirSeaFluxCO2' ,
			  'TotalAirSeaFluxCO2' ,			  
			  'Chlorophyll_cci', 			   	
			  'DIC',			  		  
			  'Nitrate',
			  'Silicate', 
			  'Temperature', 
			  'Salinity', 
			  'Oxygen',
			  'ExportRatio', 
			  'LocalExportRatio', 			  
			  'MLD',
			  'Alkalinity', 			  
			 ]
		 
		SectionTitle= 'Summary'
		hrefs 	= []
		Titles	= {}
		SidebarTitles = {}
		Descriptions= {}
		FileLists	= {}
		
		region = 'ignoreInlandSeas'
		for key in sumfields:
			href = 	key+'-'+region
			desc = ''
			if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
			if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'			

			hrefs.append(href)
			Titles[href] = 	getLongName(region) +' '+	getLongName(key)
			SidebarTitles[href] = getLongName(key)	
			Descriptions[href] = desc
			FileLists[href] = {}
			
			#####
			# Determine the list of files:
			vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
                        vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
                        vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))                        
			vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
			vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))
			vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*scatter.png'))
			vfiles.extend(glob('./images/'+jobID+'/Targets/'+year+'/*'+key+'*/BGCVal/*.png'))
			
						
		
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
				for i,t in enumerate(title):
					# remove redundent versus field
					if t.find('vs')>-1:	title[i] = ''		
			
				FileLists[href][relfn] = ', '.join([getLongName(t) for t in title])
				print "Adding ",relfn,"to script"
			
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists)
					
						
	
	
	


	if plotbyfieldandregion:
		for key in sorted(fields):
			#if key not in ['Alkalinity','Nitrate']: continue
			SectionTitle= getLongName(key)
			hrefs 	= []
			Titles	= {}
			SidebarTitles = {}
			Descriptions= {}
			FileLists	= {}
			for region in regions:
				href = 	key+'-'+region
				
				desc = ''
				if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
				if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'		
							
				hrefs.append(href)
				Titles[href] = 	getLongName(region) +' '+	getLongName(key)
				SidebarTitles[href] = getLongName(region)				
				Descriptions[href] = desc
				FileLists[href] = {}
				
				#####
				# Determine the list of files:
				vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
 	                        vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
                   	     	vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))                        
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))			
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*scatter.png'))							
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*Transect/*/*'+region+'*'+key+'*'+year+'*hov.png'))

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
						if t[:len(key)] == key: title[i] = t[len(key):]	# replace redudance value. 
						if t.find('vs')>-1:	title[i] = ''		# remove redundent versus field
				
			
					FileLists[href][relfn] = ' '.join([getLongName(t) for t in title])
					print "Adding ",relfn,"to script"
				
			html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
					SidebarTitles=SidebarTitles,#
					Titles=Titles, 
					Descriptions=Descriptions,
					FileLists=FileLists)

#			html5Tools.AddSection(indexhtmlfn,key+'-'+region,longnames, Description=longnames+' plots',Files = files)


        tar = "tar cfvz  report-"+jobID+".tar.gz "+reportdir

	print "-------------\nSuccess\ntest with:\nfirefox",indexhtmlfn
	print "To zip it up:\n",tar
	if doZip:
		import subprocess
		subprocess.Popen(tar.split())



if __name__=="__main__":	
	try:	jobID = argv[1]
	except:	
		jobID = "u-ab749"
	try: 		reportdir = argv[2]
	except: 	reportdir =folder('reports/'+jobID)
	
	try:	year = int(argv[3])
	except: year = '*'
	
	if 'clean' in argv[1:]:
		clean = True
	else:	clean = False
	html5Maker(jobID =jobID,
		   reportdir=reportdir,
		   year = year,
		   clean=clean,
		   )

















