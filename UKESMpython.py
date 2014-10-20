from sys import argv
from string import join
from os.path  import exists,getmtime
from os import mkdir, makedirs
from glob import glob
import numpy as np
"""	This is a catch all toolkit for the python methods and shorthands used in this code.
"""






def folder(name):
	""" This snippet takes a string, makes the folder and the string.
	    It also accepts lists of strings.
	"""
	if type(name) == type(['a','b','c']):
		name=join(name,'/')
	if name[-1] != '/':
		name = name+'/'
	if exists(name) is False:
		makedirs(name)
		print 'makedirs ', name
	return name


def getCommandJobIDandTime():
	jobID = argv[1]	
	timestamp = argv[2]
	return jobID,timestamp
	
def getFileList(fin):
	if type(fin)==type('abc') and fin.find('*')<0 and fin.find('?')<0: # fin is a string file:
		return [fin,]
	if type(fin)==type('abc') and (fin.find('*')>-1 or fin.find('?')>-1 or fin.find('[')>-1): # fin is a string file:
		return glob(fin)
	if type(fin) == type(['a','b','c',]): # fin is many files:
		filesout = []
		for f in fin:
			filesout.extend(glob(f))
		return filesout


def makeThisSafe(arr,log=False,debug = True, key='',noSqueeze=False):
	if noSqueeze:pass
	else: arr=np.ma.array(arr).squeeze()
	
	ma,mi = arr.max(), arr.min()
	
	if ma > 9.E36:	
		if debug: 	print "makeThisSafe: \tMasked values greater than 9.E36",key
		arr = np.ma.masked_greater(arr, 9.E36)
	
	if np.isinf(ma ) or np.isnan(ma ):
		if debug: print "makeThisSafe: \tMasking infs and Nans",key	
		arr = np.ma.array(arr)
		arr = np.ma.masked_where(np.isnan(arr)+arr.mask +np.isinf(arr) , arr)	
		
	return arr
	
def sliceA(arr,region):
	# you should already have removed time by now.
	#assume time 0.
	#if region == 'DeepestWetCells':return getDeepestWetCells(arr)
	#if region == 'Lutz': return applyLutzMask(arr)
#		region = 'Lutz':			
#		lutzmask = 	
	arr = makeThisSafe(arr,noSqueeze=True)
	

	if arr.shape[-2:] != (292,362):
		print "sliceA:\tERROR:\tThis was not designed for anything except the ORCA1 grid.", arr.shape ,"won't work."
		assert False
				 
	if len(arr.shape) == 2:

		
	    	if region == 'SEA': 		return arr[100:200,:100]
	    	elif region == 'Arctic': 	return arr[230:,45:320]
	    	elif region == 'Antarctic': 	return arr[:80,:]	
	    	elif region == 'NWAtlantic': 	return arr[205:250,220:250]
	    	elif region == 'NEPacific': 	return arr[180:243, 130:185]		    	    	
	    	elif region == 'Med': 		return arr[166:243,268:]	    	
		elif region == 'Transect': 	print 'NOT POSSIBLE TO TRANSECT A 2D FIELD'
		elif region == 'Global': 	return arr[:,:]
		elif region == 'Surface': 	return arr[:,:]
		elif region == 'NoArtics': 	return arr[80:240,:]
		elif region == 'SurfaceNoArtics': 	return arr[80:240,:]		
		elif region == 'TopLayers': 	print 'NOT POSSIBLE TO TRANSECT A 2D FIELD'						
		elif region == 'All': 		return arr[:,:]				
		elif region in ['BATS', 'HOT',]:return arr[:,::-1].squeeze().transpose()
		
	    	elif region == 'SouthPacificOcean': 	return arr[80:160,70:215]
	    	elif region == 'NorthAtlanticOcean': 	return arr[185:251,215:286,]
	    	elif region == 'SouthPacificOcean': 	return arr[80:160,70:215]
	    	elif region == 'SouthAtlanticOcean': 	return arr[72:157,238:308] 	    		    	
	    	elif region == 'NorthPacificOcean': 	return arr[170:247,57:187]	    		    	
	    	elif region == 'IndianOcean': 		return arr[88:187,-40:]    		    		    		    		    			    				 

				 
				 		    		
	elif len(arr.shape) == 3:
	    	if region == 'SEA': 		return arr[0,100:200,:100]
		elif region == 'Arctic': 	return arr[0,230:,45:320]
	    	elif region == 'Antarctic': 	return arr[:,:80,:]	
		elif region == 'Atlantic': 	return arr[::-1,:,260]
	    	elif region == 'NWAtlantic': 	return arr[:,205:250,220:250]	
	    	elif region == 'NEPacific': 	return arr[:,180:243, 130:185]	    		
	    	elif region == 'Med': 		return arr[:,166:243,268:]
		elif region == 'Transect': 	return arr[::-1,:,115]
		elif region == 'Global': 	return arr[:,:,:]
		elif region == 'Surface': 	return arr[0,:,:]
		elif region == 'NoArtics': 	return arr[:,80:240,:]
		elif region == 'SurfaceNoArtics': 	return arr[0,80:240,:]					
		elif region == 'TopLayers': 	return arr[:24,:,:]						
		elif region == 'Top200m': 	return arr[:31,:,:]
		elif region == 'Top200mNoArtics': 	return arr[:31,80:240,:]		
		elif region == 'Top40mNoArtics': 	return arr[:17,80:240,:]				
		elif region == 'Top40m': 	return arr[:17,:,:]		
		elif region == 'DeepLayers': 	return arr[24:,:,:]								
		elif region == 'All': 		return arr[:,:,:]				
		elif region in ['BATS', 'HOT',]:return arr[:,::-1].squeeze().transpose()
	    	elif region == 'SouthPacificOcean': 	return arr[:,80:160,70:215]
	    	elif region == 'NorthAtlanticOcean': 	return arr[:,185:251,215:286,]
	    	elif region == 'SouthPacificOcean': 	return arr[:,80:160,70:215]
	    	elif region == 'SouthAtlanticOcean': 	return arr[:,72:157,238:308] 	    		    	
	    	elif region == 'NorthPacificOcean': 	return arr[:,170:247,57:187]	    		    	
	    	elif region == 'IndianOcean': 		return arr[:,88:187,-40:]    
	elif len(arr.shape) == 4:
	    	if region == 'SEA':		return arr[:,0,100:200,:100].squeeze()
		elif region == 'Arctic': 	return arr[:,0,230:,45:320].squeeze()
	    	elif region == 'Antarctic': 	return arr[:,:,:80,:].squeeze()
		elif region == 'Atlantic': 	return arr[:,::-1,:,260].squeeze()
	    	elif region == 'NWAtlantic': 	return arr[:,:,205:250,220:250].squeeze()
	    	elif region == 'NEPacific': 	return arr[:,:,180:243, 130:185]	
	    	elif region == 'Med': 		return arr[:,:,166:243,268:].squeeze()	    			
		elif region == 'Transect': 	return arr[:,::-1,:,115].squeeze()
		elif region == 'Global': 	return arr[:,:,:,:].squeeze()		    			 
		elif region == 'Surface': 	return arr[:,0,:,:].squeeze()
		elif region == 'NoArtics': 	return arr[:,:,80:240,:].squeeze()
		elif region == 'SurfaceNoArtics':return arr[:,0,80:240,:].squeeze()
		elif region == 'TopLayers': 	return arr[:,0:24:,:,:].squeeze()	
		elif region == 'Top200m': 	return arr[:,:31,:,:].squeeze()						
		elif region == 'Top40m': 	return arr[:,:17,:,:].squeeze()	
		elif region == 'Top200mNoArtics': 	return arr[:,:31,80:240,:]		
		elif region == 'Top40mNoArtics': 	return arr[:,:17,80:240,:]									
		elif region == 'DeepLayers': 	return arr[:,24:,:,:].squeeze()					
		elif region == 'All': 		return arr[:,:,:,:].squeeze()		    			    							
		elif region in ['BATS', 'HOT',]:return arr[:,::-1].squeeze().transpose()
	    	elif region == 'SouthPacificOcean': 	return arr[:,:,80:160,70:215]
	    	elif region == 'NorthAtlanticOcean': 	return arr[:,:,185:251,215:286,]
	    	elif region == 'SouthPacificOcean': 	return arr[:,:,80:160,70:215]
	    	elif region == 'SouthAtlanticOcean': 	return arr[:,:,72:157,238:308] 	    		    	
	    	elif region == 'NorthPacificOcean': 	return arr[:,:,170:247,57:187]	    		    	
	    	elif region == 'IndianOcean': 		return arr[:,:,88:187,-40:]    
	else:
	    	print 'WARNINIG: PLOT SHAPE IS ODD:',arr.shape, region
	    	assert False
	print 'Could not slice', region, arr.shape
	assert False
	return arr

def shouldIMakeFile(fin,fout,debug = True):
	""" the idea is to take the file: returns:
		 True: make the file
		 False: Don't make the file.
		 
	"""
	if not exists(fout): 
		if debug: print 'shouldIMakeFile: out file doesn\'t exit and should be made.'
		return True	

	if type(fin)==type('abc') and fin.find('*')<0: # fin is a string file:
		if not exists(fin): 
			if debug: print 'Warning: ',fin ,'does not exist'
			return False 
	
		if getmtime(fin) > getmtime(fout):
			if debug: print 'shouldIMakeFile: out-file is younger than in-file, you should make it.'
			return True #
		if debug: print 'shouldIMakeFile: out-file is older than in-file, you shouldn\'t make it.'		 
		return False
	if type(fin)==type('abc') and fin.find('*')>0:
		if debug: print 'shouldIMakeFile: in-file contains *, assuming it is a wildcard: ',fin
		fin = glob(fin)
		if debug: print 'shouldIMakeFile: files : ', fin
		
	if type(fin) == type(['a','b','c',]): # fin is many files:
		for f in fin:
			if not exists(f): 
				if debug: print 'Warning: ',f ,'does not exist'
				return False 
			if getmtime(f) > getmtime(fout):
				if debug: print	'shouldIMakeFile: ',f,' is younger than an ',fout,', you should make it'
				return True
		if debug: print 'shouldIMakeFile: no new files in the list. Don\'t make it.'
		return False
	if debug:
		print	'shouldIMakeFile: got to the end somehow:'
		print type(fin), fin, fout
	return False


		
