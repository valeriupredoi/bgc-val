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
import numpy as np


def build_ensemble(timesD, arrD, ensembles={}):

	#####
	# Return no change if no ensembles requested
	if ensembles == {}: return timesD, arrD

	#####
	# Create some ensembles 
	#ensembles = {}
	#ensembles['PiControl'] = ['',]
	#ensembles['NewEmissions'] = ['u-az021', 'u-az417', 'u-az418', 'u-az513', 'u-az515', 'u-az524']
	
	#####
	# Determine time range
	timeRange = {}
	for jobID in sorted(timesD.keys()):
		try:
			minT = 	np.min(timesD[jobID])
			maxT = 	np.max(timesD[jobID])
		except: continue
		timeRange[minT] = True
		timeRange[maxT] = True		
	timeRange = [np.min(timeRange.keys()), np.max(timeRange.keys())]
	allyears = np.arange(int(timeRange[0]), int(timeRange[1])+1)
	
	#####
	# Make empty dicts:
	newTimes = {}
	newArr = {}
        outTimes = {}
        outArr = {}
        for name, ensemble in ensembles.items():
   	    	newTimes[name] = {yr:[] for yr in allyears}
   	    	newArr[name]   = {yr:[] for yr in allyears}
                outTimes[name] = []
                outArr[name]   = []

	#####
	# Load the data
	for name, ensemble in ensembles.items():
	    for jobID in sorted(timesD.keys()):	
		if jobID not in ensemble: continue
			
		for t, d in zip(timesD[jobID], arrD[jobID]):
			yr = int(t)
			newTimes[name][yr].append(t)
			newArr[name][yr].append(d)			

	#####
	# Take the mean of the ensemble	
        for name, ensemble in ensembles.items():
	    for yr in sorted(allyears):
		if len(newTimes[name][yr]) ==0: continue
                if len(newArr[name][yr]) ==0: continue

		outTimes[name].append(np.ma.mean(newTimes[name][yr]))
		outArr[name].append(np.ma.mean(newArr[name][yr]))
		
	return outTimes, outArr
	
	
	
	
	
	
