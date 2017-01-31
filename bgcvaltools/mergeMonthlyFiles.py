#!/usr/bin/ipython
#
# Copyright 2014, Plymouth Marine Laboratory
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
.. module:: mergeMonthlyFiles
   :platform: Unix
   :synopsis: A tool for stiching together multiple months of model data into one annual file. 
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import os
from re import findall

import UKESMpython as ukp 
from mergeNC import mergeNC


def getYearFromFile(fn):
	""" 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
	a = findall(r'\d\d\d\d\d\d\d\d',fn)
	mns = [ukp.mnStr(m) for m in range(1,13)]
	for i in a:
	    if i[-2:] in ['28','29','30','31','01'] and i[-4:-2] in mns: 
	    	yr = i[:4]
	      	return yr
	    			      
	return False
	
def getAnnualFilename(files, outfolder,year):
	#####
	# determinig filename on the way out
	if outfolder=='':	outfolder = ukp.folder(os.path.dirname(files[0])+'/Annual')
	
	mintime = ''
	maxtime = ''
	basefile = os.path.basename(files[0])


	for f in files:
		f = os.path.basename(f)
		a = findall(r'\d\d\d\d\d\d\d\d',f)
		if mintime==maxtime=='':	mintime,maxtime = a,a
		
		if a < mintime: mintime = a
		if a < maxtime: maxtime = a
		basefile = basefile.replace(str(a),str(year))
				
	basefile = basefile.replace('--', '-').replace('__','_')
	filenameOut = outfolder+basefile
	return 	filenameOut	

def mergeMonthlyFiles(files,outfolder='',cal='360_day'):
	#####
	# This assuemd that the files have already been split up using the moo filter tool
	# done in the the bgcvalTools/downloadFromMass.py
	
	filesOut=[]
	years = {}

	#####
	# Load file
	for fn in files:
		yr = getYearFromFile(fn)
		try:	years[yr].append(fn)
		except: years[yr] = [fn,]
	
	#####
	# 
	for yr in sorted(years.keys()):
		yearFiles= sorted(years[yr])
		if len(yearFiles)!=12:
			print "Not enough files in ",yr, len(years[yr])
			continue

		filenameOut = getAnnualFilename(years[yr], outfolder,yr)
		
		if  ukp.shouldIMakeFile(years[yr],filenameOut): 
			m = mergeNC( years[yr], filenameOut, [], timeAverage=False,debug=True,calendar=cal)
		
		filesOut.append(filenameOut)
	return filesOut

	












