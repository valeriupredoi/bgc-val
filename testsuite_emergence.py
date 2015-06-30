#!/usr/bin/ipython 
#
# Copyright 2014, Plymouth Marine Laboratory
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

#Standard Python modules:
from sys import argv

#Specific local code:
from UKESMpython import folder,getFileList,machineName
#from cchl import cchl
#from cchlvsIrradiance import cchlvsIrradiance
#from communityfit import communityfit
from emergence import cchl,cchlvsIrradiance, communityfit,primaryproduction


"""
	This is the hold all for the analyses. 
	It can be used to run each of the analyses in series, for any number of input files.
	As more analsyes are added to the package, please copy the template and add more.
	/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc
	/group_workspaces/jasmin/esmeval/bgc/MEDUSA/medusa_bio_1998.nc

"""




def main():
	try:
		filesIn = getFileList(argv[1:])
		print "Using command line arguments:", filesIn	
		
	except:
		print "testsuite:\tERROR:\tSomething wrong with file specified", argv[1:]
		return
	
	if not len(filesIn):
		print "testsuite:\tERROR:\tNo files specified, try:"
		if machineName() == 'PML':
			print "./testsuite.py /data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc"
		else:
			print "./testsuite.py /group_workspaces/jasmin/esmeval/bgc/MEDUSA/medusa_bio_1998.nc"
		
	a = primaryproduction.primaryproduction(filesIn)	
	#assert False
	for fn in filesIn:
		print "testsuite:\tINFO:\t",fn	
		a = cchl.cchl(fn)
		a = communityfit.communityfit(fn)		
		a = cchlvsIrradiance.cchlvsIrradiance(fn)	# C:Chl vs light Will only work if MEDUSA irradiance exists.








	
if __name__=="__main__":
	main() 
	print 'The end.'
