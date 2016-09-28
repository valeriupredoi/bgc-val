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
.. module:: csvFromShelves
   :platform: Unix
   :synopsis: Creates a CSV file from a shelve.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from shelve import open as shOpen
import os
from glob import glob
from math import copysign

from UKESMpython import folder

#####
# Takes a shelve or a list of shelve filenames and an output file
# This class takes a shelve file(s) from the output of the p2p analysis, and converts it into a csv file for a ValNote.
# Doesn't do any math here.

class csvFromShelves:
  def __init__(self,shelves,outFile, csvOptions = ['check',]):
  	# csv options = 'append', 'overwrite', 'check'
  	print "csvFromShelves:",outFile, csvOptions

	#####
	# Build a list of shelve files to read
  	if type(shelves) == type(['list','of','stings',]):
  		shelves_ = []
  		for s in shelves:shelves_.extend(glob(s)) 
  		shelves = shelves_
  	  	
  	if type(shelves) == type('string'):
  		shelves = glob(shelves)
  	
  	#####
  	# Does the out file exist? If so, load it, unless you're over writing it.
  	if os.path.exists(outFile):
	  	f = open(outFile,'r')
	  	txt_start = f.read()
	  	f.close()
	else:
		
		txt_start = ''
  	  
  	if len(np.intersect1d(['overwrite', 'o'], csvOptions)):
  		if os.path.exists(outFile): 
  			print "csvFromShelves:\tOVERWRITING",outFile
  		else:	print "csvFromShelves:\tCreating",outFile
  		txt_start = ''

  	if len(np.intersect1d(['check','c'], csvOptions)):check = True
  	else: check = False
  	
  	
	#####
	# Prepare to iterate over the  glob-ed shelve files.
  	txt = ''
  	for fn in shelves:
  		print "Loading", fn
  		try:
	  		sh = shOpen(fn)
	  		year 		= sh['year']
	  		name		= sh['name']
	  		region		= sh['region']
	  		newSlice 	= sh['newSlice'] 
	  		gamma	 	= sh['Taylor.gamma']
	  		E0	 	= sh['Taylor.E0']
	  		correlation 	= sh['Taylor.R']
			sh.close()
		except:
			print "unable to read",fn
			continue
		#####
		# Create some new indices:
		sig		= gamma > 1 and 1 or -1  
		inverseCorr = copysign(correlation, (np.ma.abs(1./correlation)-1.))
  		distance 	= np.sqrt(gamma**2 + E0**2 + inverseCorr**2)

		#####
		# Order to data into a string of columns
	  	columns = [year, name, region, newSlice, sig*gamma, E0, correlation, distance]  	
	  	line = ', '.join( [str(c) for c in columns])# + '\n'

		#####
		# Check to see if the line already exists in the outFile
	  	if check and txt_start.find(line) > -1:continue	  	
  		txt += line + '\n'

	#####
	# Write the outFile, according to the option flags above:
	if len(txt)==0:
		print "No changes made", csvOptions
		return

  	head = 'year, name, region, newSlice, gamma, E0, correlation, distance\n'
	if 'overwrite' in csvOptions:
		f = open(outFile,'w')
		f.write(head)

  	if len(np.intersect1d(['append', 'o','check',], csvOptions)):
  		f = open(outFile, 'w')
  		if txt_start.find(head)!=0: f.write(head)
  		f.write(txt_start)
  	
	f.write(txt)
		
	print "Saving",outFile
	
	f.close()
	
if __name__=="__main__":

	####
	# Testing script:
	fn = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-1998/chl/chl_*.shelve'
	
	outfn = folder('CSV/testing')+'test1.csv'
	#csvFromShelves(fn, outfn,csvOptions=['append',])	
	#csvFromShelves(fn, outfn,csvOptions=['overwrite',])
	csvFromShelves(fn, outfn,csvOptions=['overwrite',])				

	outfn = folder('CSV/testing')+'test_all.csv'
	csvFromShelves('/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-1998/*/*.shelve',outfn, csvOptions=['check',])		
	csvFromShelves('/data/euryale7/scratch/ledm/ukesm_postProcessed/NEMO-1998/*/*.shelve',  outfn, csvOptions=['check',])
	 
	outfn = folder('CSV/testing')+'test_all_2.csv'
	lists = ['/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-1998/*/*.shelve','/data/euryale7/scratch/ledm/ukesm_postProcessed/NEMO-1998/*/*.shelve']
	csvFromShelves(lists,outfn, csvOptions=['check',])		
	
	
	
	
			
