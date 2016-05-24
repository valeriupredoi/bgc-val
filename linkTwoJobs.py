#!/usr/bin/python 

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

from sys import argv 
import subprocess
from socket import gethostname
import os



def linkTwoJobs(jobID1,jobID2):
	"""
		command:
		linkTwoJobs jobID1 jobID2.
			job 1 is the source job
			job 2 is the new job, where the links will be put.
		
	"""
	
	if '' in [jobID1,jobID2,]: return
	
	
	machine = gethostname()
	knownmachine = False	
	if machine.find('ceda')>-1:
		knownmachine = True
		outputFold1 = "/group_workspaces/jasmin2/ukesm/BGC_data/"+jobID1
		if not os.path.exists(outputFold1):
			print "This folder doesn't exist. Do you have the right jobID 1?", jobID1
			return
		outputFold2 = "/group_workspaces/jasmin2/ukesm/BGC_data/"+jobID2
		if not os.path.exists(outputFold2):
			print "This folder doesn't exist. Do you have the right jobID 2?", jobID2
			return						
			
	if machine.find('monsoon')>-1:
		knownmachine = True
                outputFold1 = "/projects/ukesm/ldmora/UKESM/"+jobID1
		if not os.path.exists(outputFold1):
			print "This folder doesn't exist. Do you have the right jobID 1?", jobID1
			return
		outputFold2 ="/projects/ukesm/ldmora/UKESM/"+jobID2
		if not os.path.exists(outputFold2):
			print "This folder doesn't exist. Do you have the right jobID 2?", jobID2
			return	

		
	if not knownmachine :
		print "Are you running this on the correct machine?"
		print "\tYou should be on mass-cli1.ceda.ac.uk at jasmin or on monsoon at the MO"
		print "\tBut you're at",machine
		return


	linkNetcdfs = True

	########
	# Two aspects of this: 
	#	link the netcdf files
	#	copy the shelve files.
	
	if linkNetcdfs:
		
		for fn1 in sorted(glob(outputFold1+'/*.nc')):
			fn2 = fn1.replace(jobID1,jobID2)

			if os.path.exists(fn2):
				print "Already exists:\t",fn2
				continue
			try:
				os.symlink(fn1, fn2)
				print "linking:\t",fn1, '--->',fn2
			except:
				print "linking:\t",fn1, '--->',fn2, 'FAILED'
				
	


if __name__=="__main__":	
	
	try:	
		jobID1 = argv[1]
		jobID2 = argv[2]		
	except:	
		print "Please provide two jobIDs"
		jobID1 = ''
		jobID2 = ''
		
				
	linkTwoJobs(jobID1,jobID2)
		
	
	
	
