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



def downloadMass(jobID,):
	"""
		command:
		downloadMass jobID options.
		options:
			anymachine : skip check.
	"""
	if jobID == '': return
	
	machine = gethostname()
	
	if machine.find('mass')>-1:
		outputFold = "/group_workspaces/jasmin/esmeval/data/"+jobID
		if not os.path.exists(outputFold):
			print "Making ",outputFold
    			os.makedirs(outputFold)

		
	if machine.find('mass') <0 and "anymachine" not in argv:
		print "Are you running this on the correct machine?"
		print "\tYou should be on mass-cli1.ceda.ac.uk at jasmin"
		print "\tBut you're at",machine
		print "\tTo skip this warning, use the \"anymachine\" option at the command line"
		return
	
	

	
	
	print "Looking at the following files:"
	
	bashCommand = "moo ls moose:/crum/"+jobID+"/ony.nc.file/*.nc "
	print "running the command:",bashCommand
		
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	output = process.communicate()[0]



	print "Downloading at the following files:"
	
	bashCommand = "moo get --fill-gaps moose:/crum/"+jobID+"/ony.nc.file/*.nc "+outputFold
	print "running the command:",bashCommand
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	output = process.communicate()[0]





if __name__=="__main__":	
	
	try:	jobID = argv[1]
	except:	
		print "Please provide a jobID"
		jobID = ''
	downloadMass(jobID)
	
	
	
	
	
	
	
	
