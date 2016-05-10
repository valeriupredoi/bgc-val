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


	            
def html5Maker(
		jobID = 'u-ab749',
		reportdir = folder('../../html5report'),
	):

	#####
	# Delete old files
	print "Removing old files from:",reportdir
	shutil.rmtree(reportdir)

	####
	# Copy all necceasiry objects and templates to the report location:
	print " Copying html and js assets to", reportdir
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

	plotbyregion = True
	if plotbyregion:
		for key in ['Global','ignoreInlandSeas','Equator10', 'ArcticOcean','NorthernSubpolarAtlantic','NorthernSubpolarPacific','SouthernOcean','Remainder',]:
			vfiles = glob('./images/'+jobID+'/timeseries/*/percentiles*'+key+'*.png')
			files = {}
			for fn in vfiles:
				#####
				# Note that we use three paths here.
				# fn: The original file path
				# newfn: The location of the new copy of the file
				# relfn: The location of the new copy relative to the index.html
				newfn = newImageLocation(fn)	
				if not os.path.exists(newfn):
					shutil.copy2(fn, newfn)
				relfn = newfn.replace(reportdir,'./')
				
				files.append(relativePath)
				print "Adding ",newfn,"to script"
	
			html5Tools.AddSection(indexhtmlfn,'ts'+key,getLongName(key), Description=getLongName(key)+'Time Series plots',Files = files)





	print "-------------\nSuccess\ntest with:\nfirefox",indexhtmlfn





if __name__=="__main__":	
	try:	jobID = argv[1]
	except:	
		jobID = "u-ab749"
	html5Maker(jobID =jobID,)

















