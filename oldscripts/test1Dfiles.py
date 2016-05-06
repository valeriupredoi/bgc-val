# this is a temporary script to figure out the problem with the 1D files.
from ncdfView import ncdfView
import numpy as np
from matplotlib import pyplot

#fnm = "/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/CMIP5_HadGEM2-ES-CMIP5-2005-annualtophalf/oxygen/Model_oxygen__CMIP5_HadGEM2-ES-CMIP5-2005-annual_1D.nc"
#fnd = "/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/CMIP5_HadGEM2-ES-CMIP5-2005-annualtophalf/oxygen/Data_oxygen__CMIP5_HadGEM2-ES-CMIP5-2005-annual_1D.nc"
fnm = "/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/CMIP5_HadGEM2-ES-CMIP5-2005-annuallowhalf/oxygen/Model_oxygen__CMIP5_HadGEM2-ES-CMIP5-2005-annual_1D.nc"
fnd = "/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/CMIP5_HadGEM2-ES-CMIP5-2005-annuallowhalf/oxygen/Data_oxygen__CMIP5_HadGEM2-ES-CMIP5-2005-annual_1D.nc"
#fnm = "/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/CMIP5_HadGEM2-ES-CMIP5-2005-annual/oxygen/Model_oxygen__CMIP5_HadGEM2-ES-CMIP5-2005-annual_1D.nc"
#fnd = "/data/euryale7/scratch/ledm/CMIP5_postProcessed/matchedData/CMIP5_HadGEM2-ES-CMIP5-2005-annual/oxygen/Data_oxygen__CMIP5_HadGEM2-ES-CMIP5-2005-annual_1D.nc"

ncM = ncdfView(fnm,Quiet=True)
ncD = ncdfView(fnd,Quiet=True)

c = 0
notmatches=[]
doesMatch = []
model = {l:ncM(l)[:] for l in ['lat','lon','lev','index_z']}
data  = {l:ncD(l)[:] for l in ['lat','lon','depth','index_z']}
difflat = np.abs(model['lat'] - data['lat'])
difflon = np.abs(model['lon'] - data['lon'])
diffz   = np.abs(model['lev'] - data['depth'])
 
for i, (dla,dlo,dz) in enumerate(zip(difflat,difflon,diffz)):
	if dla >2. and dlo>2. and dz >500.:
		notmatches.append(i)
		continue
		print i,  m,d
		print i,'\tModel:\t',ncM('lev')[i],  '\t',ncM('lat')[i],'\t',ncM('lon')[i]
		print '\tData:\t',ncD('depth')[i],'\t',ncD('lat')[i],'\t',ncD('lon')[i]
		c+=1
		if c>12:assert False
	else:doesMatch.append(i)	
print 'notmatches:',len(notmatches), min(notmatches),max(notmatches)







for ptype in ['m35','m39','d90']:
   for whichindices in ['nomatch','doesmatch']:
	latm,lonm = [],[]
	latd,lond = [],[]
	depd,depm = [],[]
	outmatcjes=[]
	if whichindices=='nomatch': matches = notmatches
	if whichindices=='doesmatch': matches = doesMatch	
	for i in matches:
		if ptype == 'm35' and model['index_z'][i]<35:continue	
		if ptype == 'm39' and model['index_z'][i]!=39:continue
		if ptype == 'd90' and data['index_z'][i]<90:continue		
		#if data['index_z'][i]!=100:continue	
		outmatcjes.append(i)
		latm.append(model['lat'][i])
		lonm.append(model['lon'][i])	
		depm.append(model['lev'][i])
		latd.append(data['lat'][i])
		lond.append(data['lon'][i])
		depd.append(data['depth'][i])


	print "Did matching",whichindices
		
	fig = pyplot.figure()
	fig.set_size_inches(14,8)
	ax = pyplot.subplot(111,)
	pyplot.title('Model ' + whichindices + ' '+ptype)

	print 'making plot'#,lonm[:1000],latm[:1000]
	pyplot.scatter(lonm,latm,c=outmatcjes,lw=0)
	pyplot.colorbar()
	print "Saving Model:",ptype,whichindices,'model'
	pyplot.savefig('images/testData/test1d-'+ptype+whichindices+'-model.png')


	fig = pyplot.figure()
	fig.set_size_inches(14,8)
	ax = pyplot.subplot(111,)
	pyplot.title('Data '+whichindices+ ' '+ptype)
	pyplot.scatter(lond,latd,c=outmatcjes,lw=0)
	pyplot.colorbar()
	print "Saving Model:",ptype,whichindices,'data'
	pyplot.savefig('images/testData/test1d-'+ptype+whichindices+'-data.png')		


