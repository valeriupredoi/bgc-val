#!/usr/bin/ipython 
#Standard Python modules:
from sys import argv

#Specific local code:
from UKESMpython import folder,getFileList
#from cchl import cchl
#from cchlvsIrradiance import cchlvsIrradiance
#from communityfit import communityfit
from emergence import cchl,cchlvsIrradiance, communityfit,primaryproduction


"""
	This is the hold all for the analyses. 
	It can be used to run each of the analyses in series, for any number of input files.
	As more analsyes are added to the package, please copy the template and add more.
	/data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc

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
		print "./testsuite.py /data/euryale7/scratch/ledm/UKESM/MEDUSA/medusa_bio_1998.nc"

	a = primaryproduction.primaryproduction(filesIn)	

	for fn in filesIn:
		print "testsuite:\tINFO:\t",fn	
		a = cchl.cchl(fn)
		a = communityfit.communityfit(fn)		
		a = cchlvsIrradiance.cchlvsIrradiance(fn)	# C:Chl vs light Will only work if MEDUSA irradiance exists.








	
if __name__=="__main__":
	main() 
	print 'The end.'
