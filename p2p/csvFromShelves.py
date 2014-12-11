#!/usr/bin/ipython 
#
# Copyright 2014 Plymouth Marine Laboratory
#
# This file is part of the ukesm-validation library.
#
# ukesm-validation is free software: you can redistribute it and/or modify it
# under the terms of the Lesser GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version. 

# ukesm-validation is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU
# General Public License for more details.
# You should have received a copy of the Lesser GNU General
# Public License along with ukesm-validation. If not, see <http://www.gnu.org/licenses/>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.ukesm
#


import numpy as np
from shelve import open as shOpen
import os
from glob import glob
from math import copysign

from UKESMpython import folder

#####
# Takes a shelve or a list of shelve filenames and an output file
# Doesn't do any math here.

class csvFromShelves:
  def __init__(self,shelves,outFile, csvOptions = ['overwrite',]):
  	# csv options = 'append', 'overwrite', 'check'
  	print "csvFromShelves:",outFile, csvOptions
  	if type(shelves) == type('string'):
  		shelves = glob(shelves)
  	
  	if os.path.exists(outFile):
	  	f = open(outFile,'r')
	  	txt_start = f.read()
	  	f.close()
	else:
		
		txt_start = ''

  	if len(np.intersect1d(['check','c'], csvOptions)):check = True
  	else: check = False
  	  	
  	if len(np.intersect1d(['overwrite', 'o'], csvOptions)):
  		if os.path.exists(outFile): 
  			print "csvFromShelves:\tOVERWRITING",outFile
  		else:	print "csvFromShelves:\tCreating",outFile
  		txt_start = ''
  		
  	txt = ''

  	head = 'year, name, region, newSlice, gamma, E0, correlation, distance\n'
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
		
		sig		= gamma > 1 and 1 or -1  
		
		inverseCorr = copysign(correlation, (np.ma.abs(1./correlation)-1.))
  		distance 	= np.sqrt(gamma**2 + E0**2 + inverseCorr**2)

	  	columns = [year, name, region, newSlice, sig*gamma, E0, correlation, distance]
	  	
	  	line = ', '.join( [str(c) for c in columns])# + '\n'

	  	#print "line:", line, str(line)
	  	if check and txt_start.find(line) > -1:continue
	  	
  		txt += line + '\n'

 
	if 'overwrite' in csvOptions:
		f = open(outFile,'w')
		f.write(head)

  	if len(np.intersect1d(['append', 'o','check',], csvOptions)):
  		f = open(outFile, 'w')
  		if txt_start.find(head)!=0: f.write(head)
  		f.write(txt_start)
  	
	if len(txt)==0:
		print "No changes made", csvOptions
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
	csvFromShelves(fn, outfn,csvOptions=['check',])				

	outfn = folder('CSV/testing')+'test_all.csv'
	csvFromShelves('/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-1998/*/*.shelve',outfn, csvOptions=['check',])		
	csvFromShelves('/data/euryale7/scratch/ledm/ukesm_postProcessed/NEMO-1998/*/*.shelve',  outfn, csvOptions=['check',])
	 
	
	
	
	
	
			
