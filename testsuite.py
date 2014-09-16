#!/usr/bin/ipython 
#Standard Python modules:
from sys import argv
#Specific local code:
from UKESMpython import folder,getFileList
#from analysisTemplate import analysistemplate

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
		print "analysistemplate:\tERROR:\tNo file specified"
		return
		
	
	for fn in filesIn:
		#a = analysistemplate(fn)











	
if __name__=="__main__":
	main() 
	print 'The end.'
