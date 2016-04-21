
# This code removes a specific mask from the shelve.
# Untested!

from shelve import open as shopen
from glob import glob

def removeFromShelves(fn,removeRegions):
	print 'removing:',removeRegions, 'from', fn
	sh = shopen(fn)
	modeldata = sh['modeldata']

	for [(r,l,m)] in modeldata.keys():
		if r in removeRegions: print modeldata(modeldata[(r,l,m)]) , 'will be deleted'

	sh['modeldata'] = modeldata
	sh.close()

removeRegions = ['ignoreInlandSeas',]

for fn in glob('shelves/timeseries/u-ab749/u-ab749_*'):
	if fn.find('insitu')>-1:continue
	removeFromShelves(fn,removeRegions)
