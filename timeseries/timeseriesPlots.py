#!/usr/bin/ipython 

#
# Copyright 2015, Plymouth Marine Laboratory
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

from matplotlib import pyplot
import matplotlib.patches as mpatches

import numpy as np 

import timeseriesTools as tst

def trafficlights(ax,xlims, bands,labels=[]):
	if len(bands) != 6: print "wrong number of bands to traffic light."

	ax.fill_between(xlims,bands[0], bands[1] ,color='r', alpha = 0.2)		
	ax.fill_between(xlims,bands[1] ,bands[2] ,color='DarkOrange', alpha = 0.2)				
	ax.fill_between(xlims,bands[2] ,bands[3] ,color='g', alpha = 0.2)
	ax.fill_between(xlims,bands[3] ,bands[4] ,color='DarkOrange', alpha = 0.2)
	ax.fill_between(xlims,bands[4] ,bands[5] ,color='r', alpha = 0.2)	
	
	if len(labels)==5:
		patch4 = mpatches.Patch(color='r', 		alpha = 0.2,label=labels[4])
		patch3 = mpatches.Patch(color='DarkOrange', alpha = 0.2,label=labels[3])
		patch2 = mpatches.Patch(color='g', 		alpha = 0.2,label=labels[2])
		patch1 = mpatches.Patch(color='DarkOrange', alpha = 0.2,label=labels[1])
		patch0 = mpatches.Patch(color='r', 		alpha = 0.2,label=labels[0])		
		handles = [patch4,patch3,patch2,patch1,patch0,]
		pyplot.legend(handles=handles)
	return ax
	
	
def trafficlightsPlots(
		times, 			# model times (in floats)
		arr,			# model time series
		dataslice,		# in situ data distribution
		title 	='',
		filename='',
	):

	xlims= [times[0],times[-1]]
	
	fig = pyplot.figure()
	
	ax = fig.add_subplot(211)	
	pyplot.plot(times,arr)
	pyplot.xlim(xlims)	
	

	pc1 = np.array([np.percentile(dataslice,25.) for i in xlims]) 
	pc2 = np.array([np.percentile(dataslice,35.) for i in xlims])
	pc3 = np.array([np.percentile(dataslice,45.) for i in xlims])
	pc4 = np.array([np.percentile(dataslice,55.)  for i in xlims])		
	pc5 = np.array([np.percentile(dataslice,65.) for i in xlims])
	pc6 = np.array([np.percentile(dataslice,75.) for i in xlims]) 
	labels = ['25-35 pc','35-45 pc','45-55 pc','55-65 pc','65-75 pc',] 
	ax = trafficlights(ax,xlims, [pc1,pc2,pc3,pc4,pc5,pc6],labels=labels)

	
	pyplot.title(title)
	

	
	
	ax = fig.add_subplot(212)
	newt,cusum = tst.calcCuSum(times,arr)
	
	pyplot.plot(newt,cusum/np.mean(arr))
	pyplot.title('Cumulative sum / in situ mean')
	pyplot.xlim(xlims)	
	pyplot.axhline(y=0.,c='k',ls='-',lw=2,alpha=0.5)
	

	m25 = np.array([-0.25 for i in xlims]) 
	m15 = np.array([-0.15 for i in xlims])
	m05 = np.array([-0.05 for i in xlims])
	p05 = np.array([0.05  for i in xlims])		
	p15 = np.array([0.15 for i in xlims])
	p25 = np.array([0.25 for i in xlims]) 

	ax = trafficlights(ax,xlims, [m25,m15,m05,p05,p15,p25])
	
	
	
		
	print "UKESMpython:\tscatterPlot:\tSaving:" , filename
	pyplot.savefig(filename )
	pyplot.close()	
	
		
