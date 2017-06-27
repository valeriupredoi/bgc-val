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
import paths 

"""
This module includes a series of tools to download the UKESM model run data from MASS.

When run as a script, the command is::

	./downloadFromMass.py jobID

This tool will only work on machines that have mass enabled.
 
"""
def folder(name):
        """ This snippet takes a string, makes the folder and the string.
            It also accepts lists of strings.
        """
        if type(name) == type(['a','b','c']):
                name='/'.join(name,)
        if name[-1] != '/':
                name = name+'/'
        if os.path.exists(name) is False:
                os.makedirs(name)
                print 'makedirs ', name
	####
	# ensure that permissions are : drwxrwsr-x+
	#os.chmod(name , 02775)
        return name


def getYearFromFile(fn):
	""" 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
	a = findall(r'\d\d\d\d\d\d\d\d',fn)
	for i in a:
	    if i[-4:] in ['1130','1201']: 
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
	files = sorted(glob(outputFold+jobID+'o_1y_*_????_?.nc'))
	
	for fn in files:
		yr = getYearFromFile(fn)
		print "getYearFromFile:",fn, yr
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

def downloadField(jobID, keys, extension='grid[-_]T', timeslice='m',name='',dryrun=False):
	"""
	:param jobID: The job ID
	:param keys: a list of fields as they are saved in the Netcdf. (can also be a single string)
	:param timeslice: The time granularity (monthly or yearly)
	:param dryrun: does not download files, just prints.
	:param extension: Nemo style file extension
	:param name: Name of the analysis group, used for the folder.
	
	This tool takes the jobID, the field name, and using the known structure of universally similar MASS and the local filesystem structure
	from paths.py, downloads the monthly jobID data for the field requested to the local file structure.
	
	This tool will only work on machines that have mass enabled.
	
	"""
	
	if jobID == '': return
	
	#####
	# verify time granularity
	timeslice = str(timeslice).lower()
	if timeslice in ['monthly','month','1m','m']:ts = 'm'
	elif timeslice in ['yearly','year','1y','y']:ts = 'y'
	else: ts = 'y'

	#####
	# verify keys
	if type(keys) == type('string'): keys = [keys,]
	if len(keys)==1 and name=='':
		name = keys[0]

	#####
	# Verify output folder:		
	outputFold = folder(paths.ModelFolder_pref+"/"+jobID+"/"+name)
        ####
        # ensure that permissions are : drwxrwsr-x+
        os.chmod(paths.ModelFolder_pref+"/"+jobID , 02775)
        os.chmod(outputFold , 02775)

	
	print "downloadField:",name, jobID,keys,timeslice,'being saved to:',outputFold

	#####
	# make query file:
	querytxt = '-v '+' '.join(keys)
	queryfile = folder('queryfiles/')+name+'.txt'
	qf = open(queryfile,'w')
	qf.write(querytxt)
	qf.close()
	print "downloadField:\tquery text",querytxt
	
	#####
	# moose file path:
	massfp = "moose:/crum/"+jobID+"/on"+ts+".nc.file/*_1"+ts+"_*"+extension+".nc"
	print "downloadField:\tmoose path:",massfp
		
	######
	# print files
	bashCommand = "moo ls "+massfp
	print "running the command:",bashCommand
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	output = process.communicate()[0]
	print "output",output	
	
	#####
	# create a bash command to download the files:
	shell = '/bin/bash'
	bashCommand ='#!'+shell+' \n'
	for l,line in enumerate(output.split('\n')):
		if line in ['', ' ',]:continue 
		outfn = outputFold+os.path.basename(line)
		if os.path.exists(outfn):continue
                bashCommand +="\nmoo filter "+queryfile+" "+line+" "+outfn 
        bashCommand+="\necho \"The End of "+jobID +' '+name+"\"\n"
        print "running the command:\n######\n",bashCommand
        bashfile = folder('queryfiles/')+jobID+'-'+name+'.sh'
	print "Saved script as ",bashfile
        qf = open(bashfile,'w')
        qf.write(bashCommand)
        qf.close()

	
	#####
	# Run the bash command to download the files:        
        if not dryrun:
		script = shell +" "+bashfile
		print script
                process = subprocess.Popen(script.split(), stdout=subprocess.PIPE,)#shell=True, executable=shell)
                output = process.communicate()
		print "bash out",output
		
	fixFilePaths(outputFold,jobID)

######
# Some spefici wrappers for the downloadField
def nemoMonthlyIce(jobID):
	downloadField(jobID, ['soicecov',], extension='grid[-_]T', timeslice='m',name = 'monthlyIce',dryrun=False)
	
def nemoMonthlyMLD(jobID):
	downloadField(jobID, ['somxl010',], extension='grid[-_]T', timeslice='m',name = 'monthlyMLD',dryrun=False)	


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
	        ####
        	# ensure that permissions are : drwxrwsr-x+
	        os.chmod(outputFold , 02775)


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

	
        fixFilePaths(outputFold,jobID)
        deleteBadLinksAndZeroSize(outputFold,jobID)
	
	doLs = False
	doDL = True
	
	if doLs:	
		print "Looking at the following files:"
	
		bashCommand = "moo ls moose:/crum/"+jobID+"/ony.nc.file/*.nc "
		print "running the command:",bashCommand
		
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
		output = process.communicate()[0]


	if doDL:
		print "Downloading at the following files:"
		bashCommand = "moo get --fill-gaps moose:/crum/"+jobID+"/ony.nc.file/*.nc "+outputFold
		print "running the command:",bashCommand
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
		output = process.communicate()[0]

	fixFilePaths(outputFold,jobID)
	deleteBadLinksAndZeroSize(outputFold,jobID)

def fixFilePaths(outputFold,jobID):
	#####
	# The coupled model looses the first two characters of the name in the netcdf file.
	fns = glob(outputFold+"/*"+jobID[2:]+"*.nc")
	print "downloadFromMass:\tfixFilePaths:\tLooking for", outputFold+"/"+jobID[2:]+"*.nc"
	fns.extend(glob(outputFold+'/MetOffice*'))	# Because ocean assess might use the lisence?	
	for fn in sorted(fns):
		#####
        	correctfn = fn.replace(jobID[2:], jobID)
		correctfn = correctfn.replace('u-u-','u-')

	        if os.path.exists(correctfn):
	        	print "downloadFromMass:\tfixFilePaths:\tcorrect path exists.",correctfn
	        	continue
                print "downloadFromMass:\tfixFilePaths:\tFixing file prefix",
        	os.symlink(fn,correctfn)
	        print "downloadFromMass:\tfixFilePaths:\t", correctfn


        #####
        # Some runs have nemo/medusa as a preface to the file name.
	for pref in ['nemo_','medusa_']:
		#nemo_u-ai886o_1y_26291201-26301201_grid-V.nc
	        fns = glob(outputFold+"/"+pref+jobID+"*.nc")
        	print "downloadFromMass:\tfixFilePaths:\tLooking for new prefix:",pref, outputFold+"/"+pref+jobID+"*.nc"
	        for fn in sorted(fns):
        	        #####
                	correctfn = os.path.dirname(fn) +'/'+ os.path.basename(fn).replace(pref,'')
	                if os.path.exists(correctfn):
        	                print "downloadFromMass:\tfixFilePaths:\tcorrect path exists.",correctfn
                	        continue
	                print "downloadFromMass:\tfixFilePaths:\tFixing file prefix",pref,
	                os.symlink(fn,correctfn)
        	        print "downloadFromMass:\tfixFilePaths:\t", correctfn


        #####
        # Some runs have nemo/medusa as a preface to the file name.
        suffDict= {'grid-T':'grid_T',
                    'grid-U':'grid_U',
                    'grid-V': 'grid_V',
                    'grid-W':'grid_W',
                    'diad-T':'diad_T',
                    'ptrc-T':'ptrc_T',
                    #'diaptr':'diad_T', # diaptr is not the same as BGC diad-T
                    }
	for badsuff, suff in suffDict.items():
                #nemo_u-ai886o_1y_26291201-26301201_grid-V.nc
                fns = glob(outputFold+"/"+jobID+"*"+badsuff+".nc")
                print "downloadFromMass:\tfixFilePaths:\tLooking for new suff:",badsuff, outputFold+"/"+jobID+"*"+badsuff+".nc"
                for fn in sorted(fns):
                        #####
                        correctfn = os.path.dirname(fn) +'/'+ os.path.basename(fn).replace(badsuff,suff)
                        if os.path.exists(correctfn):
                                print "downloadFromMass:\tfixFilePaths:\tcorrect path exists.",correctfn
                                continue
                        print "downloadFromMass:\tfixFilePaths:\tFixing file suffix",badsuff,'->',suff,
                        os.symlink(fn,correctfn)
                        print "downloadFromMass:\tfixFilePaths:\t", correctfn

	#####
	# This code looks at symoblic links and points them at their ultimate source, removing the long link chains.
	for fn in glob(outputFold+'/*'): rebaseSymlinks(fn,dryrun=False)

def deleteBadLinksAndZeroSize(outputFold,jobID):
	

	bashCommand1 = "find "+outputFold+"/. -size 0 -print -delete"
        bashCommand2 = "find -L "+outputFold+"/. -type l -delete  -print"

	print "deleteBadLinksAndZeroSize:\t",bashCommand1
	
	process1= subprocess.Popen(bashCommand1.split(), stdout=subprocess.PIPE)
        output1= process1.communicate()[0]

        print "deleteBadLinksAndZeroSize:\t",bashCommand2


        process2= subprocess.Popen(bashCommand2.split(), stdout=subprocess.PIPE)
        output2= process2.communicate()[0]


		

if __name__=="__main__":	
	
	try:	jobID = argv[1]
	except:	
		print "Please provide a jobID"
		jobID = ''
	try:	
		keys = argv[2:]
	except:	keys = []
	
	#####
	# All yearly files
	if keys == []:	
		downloadMass(jobID) 
	#####
	# Monthly Ice files
	elif keys in [['ice',], ['soicecov',],]:
		nemoMonthlyIce(jobID)
	#####
	# Monthly MLD
	elif keys in [['mld',], ['MLD',],]:
		nemoMonthlyMLD(jobID)		
	#####
	# Other specific monthly files.
	else:
		downloadField(jobID,keys, timeslice='m',dryrun=0)
	
	
	
	
	
	
	
	
