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

from UKESMpython import folder, shouldIMakeFile
from pftnames import getLongName
from glob import glob
from sys import argv
import os 
import shutil
from html5 import html5Tools


#####
# makeReport.py:
# Usage:
#	provide a job ID, a year to look at, and a folder location.


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

def addImageToHtml(fn,imagesfold,reportdir):
	#####
	# Note that we use three paths here.
	# fn: The original file path relative to here
	# newfn: The location of the new copy of the file relative to here
	# relfn: The location of the new copy relative to the index.html
	
	newfn = imagesfold+os.path.basename(fn)
	relfn = newfn.replace(reportdir,'./')	
		
	if not os.path.exists(newfn):
		print "cp ",newfn,fn	
		shutil.copy2(fn, newfn)
	else:
		####
		# Check if the newer file is the same one from images.
		
		if os.path.getmtime(fn) == os.path.getmtime(newfn): return relfn
		####
		# Check if file is newer than the one in images.		
		if shouldIMakeFile(fn, newfn,):
			print "removing old file",fn
			os.remove(newfn)
			shutil.copy2(fn, newfn)			
			print "cp ",newfn,fn
			

	return relfn

	
         
def html5Maker(
		jobID = 'u-ab749',
		reportdir = 'reports/tmp',
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
		  'SouthernOcean',
		  'NorthernSubpolarAtlantic',
		  'NorthernSubpolarPacific',		  	
		  'Equator10', 
		  'ArcticOcean',
		  'Remainder',
		  'ignoreInlandSeas',		  
		   ]
	Transects = ['Transect','PTransect','SOTransect']
	
	#####
	# Two switches to turn on Summary section, and groups of plots of field and region.		
	Level1 = True
	Level1Regional = True
	Level1Profiles = True	
	level2Horizontal = False
	summarySections = False
	plotbyfieldandregion = False
	
	
	
	
	#####
	# A list of caveats linked to specific datasets or regions, or jobs.
	ListofCaveats, ListofCaveats_regions = {}, {}
	ListofCaveats['ExportRatio'] = 'Note that there is no historic data set for this variable.'
	ListofCaveats['Iron'] = 'Note that there is no suitable historic data set for this variable.'	
	ListofCaveats['MLD'] = 'Note that the Model MLD is calculated with based on a sigma_0 difference of 0.01 with the surface where as data uses as \
			sigma_0 difference of +/- 0.2 degrees from a depth on 10m.'
			
	if jobID == 'u-ad371':
		ListofCaveats['Chlorophyll_cci']= 'Note that the Non-diatom chlorophyll failed in run:'+jobID
		ListofCaveats['IntegratedPrimaryProduction_OSU']= 'Note that the Non-diatom chlorophyll does not contribute to IntPP in this run:'+jobID


	if Level1:
		level1Fields = [
			  'TotalIntegratedPrimaryProduction',
			  'ExportRatio', 			  
			  'AirSeaFluxCO2' ,
			  'Nitrate',
			  'DIC',			  		  
			  'Alkalinity', 
			  'TotalOMZVolume',			  			  
			 ]
		 
		SectionTitle= 'Level 1'
		hrefs 	= []
		Titles	= {}
		SidebarTitles = {}
		Descriptions= {}
		FileLists	= {}
		
		region = 'Global'
		for key in level1Fields:

			#####
			# href is the name used for the html 
			href = 	'L1'+key+'-'+region
			hrefs.append(href)
			
			#####
			# Title is the main header, SidebarTitles is the side bar title.
			Titles[href] = 	getLongName(region) +' '+	getLongName(key)
			SidebarTitles[href] = getLongName(key)	
						
			#####
			# Descriptions is a small sub-header
			desc = ''
			if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
			if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'			
			Descriptions[href] = desc
			

			#####
			# A list of files to put in this group.
			FileLists[href] = {}
			#####
			# Determine the list of files:
			vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
	                #vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
	                #vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))      
	                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+'Global*sum.png'))                                                                  
			#vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
			#vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))
			#vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*scatter.png'))
			#vfiles.extend(glob('./images/'+jobID+'/Targets/'+year+'/*'+key+'*/BGCVal/*.png'))
			
						
		
			#####
			# Create plot headers for each file.
			for fn in vfiles:
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold, reportdir)
				
				####
				# WOA fields that also produce transects, etc.
				if key in ['Nitrate', 'Silicate', 'Temperature', 'Salinity', 'Oxygen','DIC','Alkalinity'] and fn.lower().find('surface')<0:continue
				
				#####
				# Create custom title by removing extra bits.
				#title = filenameToTitle(relfn)
	
				FileLists[href][relfn] = html5Tools.fnToTitle(relfn) 
				print "Adding ",relfn,"to script"
			
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists)
					
						
	

	if Level1Regional:
		l1regions = ['Global',
		  'SouthernOcean',
		  'NorthernSubpolarAtlantic',
		  'NorthernSubpolarPacific',		  	
		  'Equator10', 
		  'ArcticOcean',
		  'Remainder',
		  'ignoreInlandSeas',		  
		   ]
		   	
		regionalFields = [
			  'Nitrate',
			  'Silicate', 
			  'Iron',
			  'IntegratedPrimaryProduction_OSU',
			]
		SectionTitle= 'Level 1 - regional'
		hrefs 		= []
		Titles		= {}
		SidebarTitles 	= {}
		Descriptions	= {}
		FileLists	= {}
		FileOrder 	= {}		
		for key in regionalFields:
			#if key not in ['Alkalinity','Nitrate']: continue

			href = 	'L1region'+key#+'-'+region
			
			desc = ''
			if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
			#if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'		
						
			hrefs.append(href)
			Titles[href] = 		getLongName(key)
			SidebarTitles[href] = getLongName(key)				
			Descriptions[href] = desc
			FileLists[href] = {}
			FileOrder[href] = {}
			#####
			# Determine the list of files:
			vfiles = []
			for region in l1regions:				
				vfiles.extend(glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png'))

			#####
			# Create plot headers for each file.
			count=0
			for fn in vfiles:
				#####
				# Skip transects, they'll be added below.
				if fn.find('Transect') >-1: continue
				if fn.lower().find('surface')<0:continue
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold, reportdir)
			
				#####
				# Create custom title by removing extra bits.
				title = html5Tools.fnToTitle(relfn)
		
				FileLists[href][relfn] = title
				FileOrder[href][count] = relfn
				count+=1
				print "Adding ",relfn,"to script"

				
					
				
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists,
				FileOrder=FileOrder)



	if Level1Profiles:
	    for plottype in ['profile','profilehov']:
		l1regions = ['Global',
		  'SouthernOcean',
		  'NorthernSubpolarAtlantic',
		  'NorthernSubpolarPacific',		  	
		  'Equator10', 
		  'ArcticOcean',
		  'Remainder',
		  'ignoreInlandSeas',		  
		   ]
		   	
		regionalFields = [
			  'Nitrate',
			  'Silicate', 
			  'Iron',
			  'DIC',
			  'Alkalinity',
			  'Oxygen',
			]
		if plottype == 'profile':	SectionTitle= 'Level 1 - Profiles'
		if plottype == 'profilehov':	SectionTitle= 'Level 1 - Hovmoeller plots'		
		hrefs 		= []
		Titles		= {}
		SidebarTitles 	= {}
		Descriptions	= {}
		FileLists	= {}
		FileOrder 	= {}		
		for key in regionalFields:
			#if key not in ['Alkalinity','Nitrate']: continue

			href = 	'L1'+plottype+'-'+key#+'-'+region
			
			desc = ''
			if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
			#if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'		
						
			hrefs.append(href)
			Titles[href] = 		getLongName(key)
			SidebarTitles[href] = getLongName(key)				
			Descriptions[href] = desc
			FileLists[href] = {}
			FileOrder[href] = {}
			#####
			# Determine the list of files:
			vfiles = []
			for region in l1regions:
				#vfiles.extend(glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png'))
		                if plottype == 'profile':	vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile_*'+key+'*'+region+'*mean.png'))
		                if plottype == 'profilehov':	vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profilehov_*'+key+'*'+region+'*mean.png'))		                
			#####
			# Create plot headers for each file.
			count=0
			for fn in vfiles:
				#####
				# Skip transects, they'll be added below.
				if fn.find('Transect') >-1: continue
				#if fn.lower().find('surface')<0:continue
				
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold, reportdir)
			
				#####
				# Create custom title by removing extra bits.
				title = html5Tools.fnToTitle(relfn)
		
				FileLists[href][relfn] = title
				FileOrder[href][count] = relfn
				count+=1
				print "Adding ",relfn,"to script"
				
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists,
				FileOrder=FileOrder)
				
				
				

		
	if level2Horizontal:
		l2Fields = [
			  'Nitrate',
			  'Silicate', 		  
			  'DIC',			  		  
			  'Alkalinity',
			  'Chlorophyll_cci', 			   	
			  #'TotalIntegratedPrimaryProduction',
			  'IntegratedPrimaryProduction_OSU', 
			  #'AirSeaFluxCO2' ,
			  #'TotalOMZVolume',
			  #'TotalAirSeaFluxCO2' ,
			  			  
			  #'Temperature', 
			  #'Salinity', 
			  #'Oxygen',
			  #'ExportRatio', 
			  #'LocalExportRatio', 			  
			  #'MLD',
			 ]
		 
		SectionTitle= 'Level 2'
		hrefs 	= []
		Titles	= {}
		SidebarTitles = {}
		Descriptions= {}
		FileLists	= {}
		
		region = 'ignoreInlandSeas'
		for key in sumfields:

			#####
			# href is the name used for the html 
			href = 	'L2'+key+'-'+region
			hrefs.append(href)
			
			#####
			# Title is the main header, SidebarTitles is the side bar title.
			if key == 'SummaryTargets':
				Titles[href] = 	getLongName(key)			
			else:	Titles[href] = 	getLongName(region) +' '+	getLongName(key)
			SidebarTitles[href] = getLongName(key)	
						
			#####
			# Descriptions is a small sub-header
			desc = ''
			if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
			if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'			
			Descriptions[href] = desc
			


			#####
			# A list of files to put in this group.
			FileLists[href] = {}
			if key == 'SummaryTargets':
				vfiles = glob('./images/'+jobID+'/Targets/'+year+'/Summary/*'+region+'*.png')			
			else:
				#####
				# Determine the list of files:
				vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
		                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
		                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))      
		                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+'Global*sum.png'))                                                                  
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*scatter.png'))
				vfiles.extend(glob('./images/'+jobID+'/Targets/'+year+'/*'+key+'*/BGCVal/*.png'))
			
						
		
			#####
			# Create plot headers for each file.
			for fn in vfiles:
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold, reportdir)
				
				####
				# WOA fields that also produce transects, etc.
				if key in ['Nitrate', 'Silicate', 'Temperature', 'Salinity', 'Oxygen','DIC','Alkalinity'] and fn.lower().find('surface')<0:continue
				
				#####
				# Create custom title by removing extra bits.
				#title = filenameToTitle(relfn)
	
				FileLists[href][relfn] = html5Tools.fnToTitle(relfn) 
				print "Adding ",relfn,"to script"
			
		html5Tools.AddSubSections(indexhtmlfn,hrefs,SectionTitle,
				SidebarTitles=SidebarTitles,#
				Titles=Titles, 
				Descriptions=Descriptions,
				FileLists=FileLists)
					
						
						
						
						

			
		
	if summarySections:
		sumfields = [
			  'SummaryTargets',
			  'TotalIntegratedPrimaryProduction',
			  'IntegratedPrimaryProduction_OSU', 
			  'AirSeaFluxCO2' ,
			  'TotalOMZVolume',
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

			#####
			# href is the name used for the html 
			href = 	key+'-'+region
			hrefs.append(href)
			
			#####
			# Title is the main header, SidebarTitles is the side bar title.
			if key == 'SummaryTargets':
				Titles[href] = 	getLongName(key)			
			else:	Titles[href] = 	getLongName(region) +' '+	getLongName(key)
			SidebarTitles[href] = getLongName(key)	
						
			#####
			# Descriptions is a small sub-header
			desc = ''
			if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
			if region in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'			
			Descriptions[href] = desc
			


			#####
			# A list of files to put in this group.
			FileLists[href] = {}
			if key == 'SummaryTargets':
				vfiles = glob('./images/'+jobID+'/Targets/'+year+'/Summary/*'+region+'*.png')			
			else:
				#####
				# Determine the list of files:
				vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
		                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
		                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))      
		                vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+'Global*sum.png'))                                                                  
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*hist.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*robinquad.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*/*/*'+region+'*'+key+'*'+year+'*scatter.png'))
				vfiles.extend(glob('./images/'+jobID+'/Targets/'+year+'/*'+key+'*/BGCVal/*.png'))
			
						
		
			#####
			# Create plot headers for each file.
			for fn in vfiles:
				#####
				# Copy image to image folder and return relative path.
				relfn = addImageToHtml(fn, imagesfold, reportdir)
				
				####
				# WOA fields that also produce transects, etc.
				if key in ['Nitrate', 'Silicate', 'Temperature', 'Salinity', 'Oxygen','DIC','Alkalinity'] and fn.lower().find('surface')<0:continue
				
				#####
				# Create custom title by removing extra bits.
				#title = filenameToTitle(relfn)
	
				FileLists[href][relfn] = html5Tools.fnToTitle(relfn) 
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
				#vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+'*Transect/*/*'+region+'*'+key+'*'+year+'*hov.png'))

				#####
				# Create plot headers for each file.
				for fn in vfiles:
					#####
					# Skip transects, they'll be added below.
					if fn.find('Transect') >-1: continue
					
					#####
					# Copy image to image folder and return relative path.
					relfn = addImageToHtml(fn, imagesfold, reportdir)
				
					#####
					# Create custom title by removing extra bits.
					title = html5Tools.fnToTitle(relfn)
			
					FileLists[href][relfn] = title
					print "Adding ",relfn,"to script"
			for transect in Transects:
				href = 	key+'-'+transect
				
				desc = ''
				if key in ListofCaveats.keys():			desc +=ListofCaveats[key]+'\n'
				if transect in ListofCaveats_regions.keys():	desc +=ListofCaveats_regions[key]+'\n'		
							
				hrefs.append(href)
				Titles[href] = 	getLongName(transect) +' '+	getLongName(key)
				SidebarTitles[href] = getLongName(transect)				
				Descriptions[href] = desc
				FileLists[href] = {}
				
				#####
				# Determine the list of files:
				region = 'Global'
				vfiles = []
				#vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*'+region+'*10-90pc.png')
 	                        #vfiles.extend(glob('./images/'+jobID+'/timeseries/*/profile*'+key+'*'+region+'*median.png'))
                   	     	#vfiles.extend(glob('./images/'+jobID+'/timeseries/*/Sum*'+key+'*'+region+'*sum.png'))                        
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+transect+'/*/'+key+transect+'*'+region+'*'+key+'*'+year+'*hist.png'))
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+transect+'/*/'+key+transect+'*'+region+'*'+key+'*'+year+'*robinquad.png'))			
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+transect+'/*/'+key+transect+'*'+region+'*'+key+'*'+year+'*scatter.png'))		
				vfiles.extend(glob('./images/'+jobID+'/P2Pplots/*/*'+key+transect+'/*/'+key+transect+'*'+region+'*'+key+'*'+year+'*hov.png'))

				#####
				# Create plot headers for each file.
				for fn in vfiles:
					#####
					# Copy image to image folder and return relative path.
					relfn = addImageToHtml(fn, imagesfold, reportdir)
				
					#####
					# Create custom title by removing extra bits.
					title = html5Tools.fnToTitle(relfn)
			
					FileLists[href][relfn] = title
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



def main():
	try:	jobID = argv[1]
	except:	
		print "Please provide a jobID next time"
		exit()
	
	#defaults:
	clean = False
	year = '*'
	reportdir =folder('reports/'+jobID)
	
	for i,arg in enumerate(argv):
		if i <= 1: continue

		try:	
			y = int(arg)
			year = y
			continue
		except: pass

		if arg == 'clean':
			 clean = True
			 continue
		
		reportdir = arg
			

	

	html5Maker(jobID =jobID,
		   reportdir=reportdir,
		   year = year,
		   clean=clean,
		   )
		   	
if __name__=="__main__":	
	main()

















