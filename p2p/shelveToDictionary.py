#
# Copyright 2016, Plymouth Marine Laboratory
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

from UKESMpython import reducesShelves,listShelvesContents
from shelve import open as shOpen


def shelveToDictionary(shelvesAV):
	"""
		In this code, we take the shelvesAV dictionary from the analysis-JASMIN.py,
		and returns a dictionary of the metrics.
	"""
	

	allshelves = listShelvesContents(shelvesAV)
	models 	= allshelves.models
	names  	= allshelves.names
	dls     = allshelves.depthLevels
	years	= allshelves.years
	slices 	= allshelves.sliceslist
	
	outDict={}
	openned = 0
	for model in models:
	  for name in names:
	    for dl in dls:			  
	      for year in years:
	      	for sl in slices:

			shelves = reducesShelves(shelvesAV,  models =[model,],names = [name,],depthLevels = [dl,], years=[year,], sliceslist =[sl,])
			if len(shelves)<1: 
			#	print metricname, 'found nothing'
				continue
			if len(shelves)>1: 
			#	print metricname, ': too many', shelves
				continue
											
			#if len(shelves)!=1: continue
			

			print "--------------\noutPutForJASMIN",[model,name,dl,sl],':',shelves
			
			s = shOpen(shelves[0])
			openned+=1
			for key in [	'N', 'b0', 'stdErr', 'b1', 'pValue', 'rValue',
					'Taylor.E0', 'Taylor.E', 'Taylor.gamma',  'Taylor.R', 'Taylor.p',
					'robust.E0', 'robust.E', 'robust.gamma',  'robust.R', 'robust.p',
					'MNAFE', 'MNFB', 'NMAEF', 'NMBF', ]:
				metricname = '-'.join([model,name,dl,sl,key])	
				if len(	years)>1: metricname = year+'-'+metricname	
				outDict[metricname] = s[key]
				print metricname,':\t',outDict[metricname]
	
	print "And the output dictionary is:",outDict, "\n(still probably too big)"
	
	return outDict
