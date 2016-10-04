#!/usr/bin/ipython 
#
# Copyright 2014 Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Lesser GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, 
# or (at your option) any later version. 

# bgc-val is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU
# General Public License for more details.
# You should have received a copy of the Lesser GNU General
# Public License along with bgc-val. If not, see <http://www.gnu.org/licenses/>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.ukesm
#

#Standard Python modules:
from sys import argv
from os.path import exists
from calendar import month_name

#Specific local code:
from UKESMpython import folder,getFileList, AutoVivification, NestedDict,AutoVivToYaml,YamlToDict
from p2p import matchDataAndModel,makePlots,makeTargets

from bgcvaltools.pftnames import MaredatTypes,WOATypes,Ocean_names

#	
#
# Here is the plan:
#	
# 	you want to find the one point where the lat/lon match, then extract that point from your data.
#		is it better to have one 
# 	then m


# setp one#
jobID = 'xhonp'
station = 'HOTS'
	
coords = AutoVivification()
coords['HOTS']['lat'] = 22.75 	# 22deg 45'N,
coords['HOTS']['lon'] = -158. 	# 158deg 00'W

coords['BATS']['lat'] = 32.34 	# 32deg 19' 59" N / 
coords['BATS']['lon'] = -64.75 	# 64deg 45' 0" W

	
	
	
if __name__=="__main__":
	# Can use command line arguments to choose a model.
	models 		= []
	years 		= []
	ERSEMjobIDs 	= []
	
	#####
	# Determine command line arguments
	for a in argv[1:]:	
		try:	
			y = int(a)
			years.append(a)
			continue			
		except:pass
		
		if str(a).upper() in ['MEDUSA','ERSEM','NEMO']:
			models.append(str(a).upper())
			continue			
		if a[:4] in ['xhon','xjez']:
			ERSEMjobIDs.append(a)
			if 'ERSEM' not in models:models.append('ERSEM')
			continue
			
		print "Command line argument not understood:",a
	

	#####
	#Set Defaults:
	if not len(years): 	years = ['1998',]
	if not len(models): 	models = ['MEDUSA','ERSEM','NEMO']
	if not len(ERSEMjobIDs):ERSEMjobIDs = ['xhonp',]	

	print "#############################"
	print "__main__ arguments: "
	print "models:        ",models
	print "year:          ",years
	print "ERSEM jobID:   ",ERSEMjobIDs
	print "#############################"
	#sleep(20)


			
	for year in years:
		#multiERSEMtargets(['yaml/shelvesAVERSEM'+year+e+'.yaml' for e in ERSEMjobIDs],ERSEMjobIDs)	
		#continue
		testsuite_p2p(models = models,	year=year,ERSEMjobID=ERSEMjobIDs[0] ) 
		if len(ERSEMjobIDs)==1:continue
		for e in ERSEMjobIDs[1:]:
			testsuite_p2p(models = ['ERSEM',],year=year,ERSEMjobID=e ) 
	
	print 'The end.'
	






