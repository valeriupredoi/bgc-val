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
"""
.. module:: downloadFromMass
   :platform: Unix
   :synopsis: A set of tools to download the UKESM model run data from MASS.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from sys import argv 
import subprocess
from socket import gethostname
import os
from glob import glob
from re import findall


"""
This module includes a series of tools to download the UKESM model run data from MASS.

When run as a script, the command is::

	./downloadFromMass.py jobID

This tool will only work on machines that have mass enabled.
 
"""

def getYearFromFile(fn):
	""" 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
	a = findall(r'\d\d\d\d\d\d\d\d',fn)
	for i in a:
	    if i[-4:] == '1130': 
	    	yr = i[:4]
	      	return yr
			      
	return False
	
def rebaseSymlinks(fn,dryrun=True,debug=False):
	""" 
	:param fn: A full path to a filename. It should be a symbolic link.
	:param dryrun: A boolean switch to do a trial run of this function.

	This function reduces a chain of symbolic links down to one. It takes a full path, 
	checks whether it is a sym link, then checks whether the real path is the  target of the link.
	If not, it replaces the target with the real path.
	
	"""
	
        #####
        # fn is a link file
        #if not os.path.exists(fn):
        #       print "rebaseSymlinks:\tfile does not exist.",fn
        #       return
        if not os.path.islink(fn):
                if debug:print "rebaseSymlinks:\tfile is not a symlink.",fn
                return
        #####
        # Real path and first target:
        realpath = os.path.realpath(fn)         # The final end of the link chin
        linkpath = os.readlink(fn)              # The first target in the symlink chain

        if realpath == linkpath: return

        print "rebaseSymlinks:\tdeleting and re-linking ",fn,'-->', realpath
        if dryrun:      return
        os.remove(fn)
        os.symlink(realpath,fn)

	
def findLastFinishedYear(jobID,dividby=1,numberfiles=6):
	"""
	:param jobID: The job ID, as elsewhere.
	:param 	dividby: Outputs every "dividby" years.
	:param numberfiles: The expected number of files per a year. (usually 6, but sometimes 4)

	This tool find the best year to have a close look at the model, by searching through the files
	and guessing which years are finished.
	
	"""
	if jobID == '': return
	
	machine = gethostname()
	if machine.find('ceda')>-1:
		outputFold = "/group_workspaces/jasmin2/ukesm/BGC_data/"+jobID+'/'

	if machine.find('monsoon')>-1:
                outputFold = "/projects/ukesm/ldmora/UKESM/"+jobID+'/'
                        
	if gethostname().find('pmpc')>-1:	
                outputFold = "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+'/'		
		                        
	fnDict = {}	
	files = sorted(glob(outputFold+jobID+'o_1y_????1201_????1130_????_?.nc'))
	#suffixes = ['diad_T.nc', 'grid_T.nc','grid_U.nc','grid_V.nc','grid_W.nc','ptrc_T.nc']
	for fn in files:
		yr = getYearFromFile(fn)
		print fn, yr
		try: 	fnDict[yr]+=1
		except:	fnDict[yr] =1
	
	years = sorted(fnDict.keys())
	years.reverse()
	
	print years, fnDict

	if len(years) == 0:
		print "findLastFinishedYear:\tNo files found.\t"
		return False
			
	if len(years)< dividby:
		print "findLastFinishedYear:\tLess than",dividby,"years of model run, returning first year:",years[-1]
		return years[0]
	
	for y in years:
		if int(y)%dividby != 0: continue
		
		print y,':', fnDict[y]
		if fnDict[y] >= numberfiles: return y
		
	print "No correct year, there's probably a problem here findLastFinishedYear(",jobID,")"
	print "Machine", machine
	print "outputFold:", outputFold
	return False
	#assert 0	


def downloadMass(jobID,):
	"""
	:param jobID: The job ID
	
	This tool takes the jobID, and using the known structure of universally similar MASS and the local filesystem structure
	from paths.py, downloads the jobID data to the local file structure.
	
	This tool will only work on machines that have mass enabled.
	
	"""
	if jobID == '': return
	
	machine = gethostname()
	knownmachine = False	
	if machine.find('ceda')>-1:
		knownmachine = True
		outputFold = "/group_workspaces/jasmin2/ukesm/BGC_data/"+jobID
		if not os.path.exists(outputFold):
			print "Making ",outputFold
    			os.makedirs(outputFold)

	if machine.find('monsoon')>-1:
		knownmachine = True
                outputFold = "/projects/ukesm/ldmora/UKESM/"+jobID
                if not os.path.exists(outputFold):
                        print "Making ",outputFold
                        os.makedirs(outputFold)

		
	if not knownmachine :
		print "Are you running this on the correct machine?"
		print "\tYou should be on mass-cli1.ceda.ac.uk at jasmin or on monsoon at the MO"
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


	#####
	# The coupled model looses the first two characters of the name in the netcdf file.
	fns = glob(outputFold+"/"+jobID[2:]+"*.nc")
	print "Looking for", outputFold+"/"+jobID[2:]+"*.nc"
	fns.extend(glob(outputFold+'/MetOffice*'))	# Because ocean assess might use the lisence?	
	for fn in sorted(fns):
		#####
        	correctfn = fn.replace('/'+jobID[2:], '/'+jobID)
	        if os.path.exists(correctfn):
	        	print "correct path exists.",correctfn
	        	continue
                print "Fixing file prefix",
        	os.symlink(fn,correctfn)
	        print correctfn
	
	#####
	# This code looks at symoblic links and points them at their ultimate source, removing the long link chains.
	for fn in glob(outputFold+'/*'): rebaseSymlinks(fn,dryrun=False)

if __name__=="__main__":	
	
	try:	jobID = argv[1]
	except:	
		print "Please provide a jobID"
		jobID = ''
	downloadMass(jobID)
	
	
	
	
	
	
	
	
