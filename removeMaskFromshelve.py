
# This code removes a specific mask from the shelve.
# Untested!

from shelve import open as shopen
sh = shopen(fn)

modeldata = sh['modeldata']

removeRegions = ['ignoreInlandSeas',]

for [(r,l,m)] in modeldata.keys():
	if r in removeRegions: print modeldata(modeldata[(r,l,m)]) , 'will be deleted'
