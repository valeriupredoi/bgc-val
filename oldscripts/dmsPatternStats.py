#!/usr/bin/ipython 
#
# Copyright 2014, Plymouth Marine Laboratory
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
#

from matplotlib import pyplot
from calendar import month_name
from glob import glob
from shelve import open as shOpen
from bgcvaltools.StatsDiagram import rmsds
from UKESMpython import folder
from itertools import product
import os
import numpy as np
from ncdfView import ncdfView
from pyproj import Proj
from shapely.geometry import shape

purple = [125./256., 38./256., 205./256.]

oceans = [ 'AntarcticOcean', 'IndianOcean', 'SouthPacificOcean', 'NorthPacificOcean',
	   'SouthAtlanticOcean','NorthAtlanticOcean',
	   'ArcticOcean',] 
months = [month_name[i] for i in xrange(1,13)]	  
hemispheres	=['NorthHemisphere','SouthHemisphere',]
seasons = ['JFM','AMJ','JAS','OND',]

#ncGfn = 'data/Flat1deg-Areas.nc'
#ncG = ncdfView(ncGfn,Quiet=True)
#area = ncG('area')[:]


# To run this from testsuite you'll need:
#	move it to p2p
#	call it with:
#	list of shelves
#	grid (ie flat1deg.)# only for volume field
#	x-axis choice (ie months, oceans, etc)
#	color choice ie, what to actually plot on each subplot

def calculateArea():
	# makes a netcdf with the areas of each pixel in the flat map.
	nc = ncdfView('data/Flat1deg.nc',Quiet=True)
	lat = nc('lat')[:]
	#lon = nc('lon')[:]	 
	areasDict = {}
	print 'calculating area: '	
	for i,l in enumerate(lat): #(u'lat', u'lon')
		co = {"type": "Polygon", "coordinates": [
		    [(1., l+0.5), #('lon', 'lat')
		     (1., l-0.5),
		     (0., l-0.5),
		     (0., l+0.5)]]}
		clon, clat = zip(*co['coordinates'][0])

		pa = Proj("+proj=aea +lat_1="+str(l-0.5)+" +lat_2="+str(l+0.5)+"+lat_0="+str(l)+" +lon_0=0.5")		
		x, y = pa(clon, clat)
		cop = {"type": "Polygon", "coordinates": [zip(x, y)]}

		area = shape(cop).area  

		areasDict[l] = area
	nc.close()
	return areasDict				
areaDict = calculateArea()

def getArea(sh, key = 'datax'):
	areasum = 0.
	for x,dms in zip(sh['x_lat'], sh[key]):
		areasum += areaDict[x]*dms
	return areasum

def run(key_dmsmodels, key_xkeys,):
	
	if key_dmsmodels =='dmsmetrics': 	dmsmodels = ['dms_and','dms_ara','dms_hal','dms_sim','Lana',]
	if key_dmsmodels =='dmspmetrics': 	dmsmodels = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim','Lana',]	
	
	xkeys = []
	if key_xkeys in ['months','oceans',]:	
		shelvefold = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-iMarNet-2526/'		
	else:	shelvefold = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-xkrum-2526/'	
		
	if key_xkeys == 'months':	xkeys 	= months
	if key_xkeys == 'oceans':	xkeys	= oceans
	if key_xkeys == 'seasons':	xkeys	= seasons
	if key_xkeys == 'hemispheres':	xkeys	= hemispheres	

	if key_xkeys == 'OceanMonths':		
		xkeys	= sorted([i[0]+i[1] for i in product(oceans,months)])

	if key_xkeys == 'OceanSeasons':		
		xkeys	= sorted([i[0]+i[1] for i in product(oceans,seasons)])
		
	if key_xkeys == 'HemispheresMonths':		
		xkeys	= sorted([i[0]+i[1] for i in product(hemispheres,months)])


		
	if key_xkeys not in ['months', 'OceanMonths'] and key_xkeys.lower().find('months')>0:
		i = key_xkeys.lower().find('months')	
		xkeys = [key_xkeys[:i]+' '+m for m in months]

	if key_xkeys not in ['seasons', 'OceanSeasons'] and key_xkeys.lower().find('seasons')>0:
		i = key_xkeys.lower().find('seasons')	
		xkeys = [key_xkeys[:i]+' '+s for s in seasons]
		
	modellongname = {m:i for m,i in zip(dmsmodels,['Anderson','Aranami','Halloran','Simo','Lana'])}
	modelmarkers = {m:i for m,i in zip(dmsmodels,['s','<','o','D','v'])}
	modelcolours = {m:i for m,i in zip(dmsmodels,['r','g','b',purple,'k'])}

	bias = {}
	stds = {}
	rcor = {}
	corr = {}
	numb = {}
	beta0 = {}
	beta1 = {}
	

	metricslist =  { 'rgam':'Scale Ratio (Data/Reference)',
			 'rE':'Norm. Diff. Scale',
			 'rE0':'Norm. Bias',
			 'rR':'Spearman Corr.',
			 'rP':'p-value (robust)',
			 #'rsig':'sign',
			 'tgam':'STD Ratio (Data/Reference)',
			 'tE':'Norm. Unbiased RMSD',
			 'tE0':'Norm. Bias',
			 'tR':'Pearson  Corr.',
			 'tP':'p-value (Taylor)',
			 #'tsig':'sign',			 
			 'N':'Number',
			 'b0':'Intersect',
			 'b1':'Slope',
			 'TotalDMS':'Total DMS  Mmol / m',
			 'DMS:Lana':'DMS / Lana ratio',			 
			}
						
	metrics = {me:{} for me in metricslist.keys()}

		
	i = 0

	for o in xkeys:
		for me in metricslist.keys(): metrics[me][o] = {}

		for m in dmsmodels:
			if m == 'Lana':	
				for met in metricslist.keys():	metrics[met][o][m] = np.ma.masked
				gfn = shelvefold+dmsmodels[0]+'Surface/'+dmsmodels[0]+'_'+o.replace(' ','')+'_*.shelve'
				fn = glob(gfn)
				print 'globbing', gfn
				if len(fn)==1:fn = fn[0]
				print i, o,m,'opening',fn
				if not fn:
					for me in metricslist.keys(): metrics[me][o][m] = np.ma.masked
					i+=1
					continue
				sh = shOpen(fn)
				metrics['TotalDMS'][o][m] = getArea(sh,key='datay')/1.e12
				metrics['DMS:Lana'][o][m] = 1.		
				continue			
			gfn = shelvefold+m+'Surface/'+m+'_'+o.replace(' ','')+'_*.shelve'
			fn = glob(gfn)
			print 'globbing', gfn			
			if len(fn)==1:fn = fn[0]
			print i, o,m,'opening',fn
			if not fn:
				for me in metricslist.keys(): metrics[me][o][m] = np.ma.masked
				i+=1
				continue
			sh = shOpen(fn)
			#print sh.keys()
			#assert False
			sig =  sh['robust.gamma']>1 and 1 or -1
			metrics['rgam'][o][m] 	= sh['robust.gamma']
			metrics['rE0'][o][m]  	= sh['robust.E0'] 
			metrics['rE'][o][m]  	= sh['robust.E'] *sig
			metrics['rR'][o][m]	= sh['robust.R'] 			
			metrics['rP'][o][m]	= sh['robust.p'] 						

			sig =  sh['Taylor.gamma']>1 and 1 or -1			
			metrics['tgam'][o][m] 	= sh['Taylor.gamma']
			metrics['tE0'][o][m]  	= sh['Taylor.E0'] 
			metrics['tE'][o][m]  	= sh['Taylor.E'] * sig			
			metrics['tR'][o][m]   	= sh['Taylor.R'] 
			metrics['tP'][o][m]   	= sh['Taylor.p'] 
			
			
			for k in ['N','b1','b0']:
				metrics[k][o][m] = sh[k]
			
			metrics['TotalDMS'][o][m] = getArea(sh)	/1.e12	
			metrics['DMS:Lana'][o][m] = metrics['TotalDMS'][o][m] / (getArea(sh,key='datay')/1.e12)
			sh.close()
			i+=1			
	title = ''#key_dmsmodels
	if 'dms_p_and' in dmsmodels: 
		title += 'DMS (pixels)'
		imagefold = folder('images/MEDUSA-xkrum/DMS_Patterns/pixels')
	if 'dms_and' in dmsmodels: 
		title += 'DMS (extrapolated)'	
		imagefold = folder('images/MEDUSA-xkrum/DMS_Patterns/extrapolated')		

	#if 'SouthPacificOcean' in xkeys: title += ' Oceans'
	#if 'January' in xkeys: 		 title += ' Months'
	
	


				
	plotTypes = {	'Total N':   ['TotalDMS','N',],
		#	'Total':        ['TotalDMS',],
			'Total':   ['TotalDMS','DMS:Lana'],
						
			'Robust':['rE','rE0','rR',],
			'Taylor':['tE','tE0','tR',],
			#'Robust v Taylor Correlation':['rR','tR',],
			#'Robust v Taylor Correlation and N':['rR','tR','N',],
			#'Robust v Taylor Gamma':['rgam','tgam',],
			#'Robust v Taylor E':['rE','tE',],
			#'Robust v Taylor E0':['rE0','tE0',],
			'Linear Regression':['b1','b0','tR',],
			}
	
	for plotType,metrickeys in plotTypes.items():
		plotStyle = 'Lines'	
		#if plotType == 'TotalDMS_N': plotStyle = 'Lines'
		#if plotType == 'TotalDMS':   plotStyle = 'Lines'		
		#if plotType == 'TotalRatio': plotStyle = 'Lines'
		fn = imagefold+key_dmsmodels+'_'+ key_xkeys+'_'+plotType+'.png'
		fn = fn.replace(' ', '')
		#if os.path.exists(fn):continue
				
		fig = pyplot.figure()
		for d,metric in enumerate(metrickeys):
			ax = fig.add_subplot(len(metrickeys),1,d+1)
			xticks = []
			xticklabels=[]

			linesDict = {'x':[]}
			for m in dmsmodels:linesDict[m] = []
			for o,slice in enumerate(xkeys):
				xticks.append(o)
				linesDict['x'].append(o)
				xticklabels.append(slice.replace('Ocean','').replace(' ','\n').replace('North','North ').replace('South','South '))
				for m, val in sorted(metrics[metric][slice].items()):
					if plotStyle=='Lines':
						linesDict[m].append(val)
					else:
						if np.ma.is_masked(val): continue
						pyplot.scatter(o,val,marker=modelmarkers[m],c=modelcolours[m], s = 40,lw=0.)

			if plotStyle=='Lines':
			    for m in dmsmodels: pyplot.plot(linesDict['x'],linesDict[m],c=modelcolours[m],lw=2)# s = 40,lw=0.)#marker=modelmarkers[m]
			
			pyplot.ylabel(metricslist[metric])		
						
			
			if metric ==metrickeys[-1]: 	
				pyplot.xticks(xticks, xticklabels,rotation = 'vertical',)
			else:	pyplot.xticks(xticks, ['',])
			
			if d+1 == 1:	ax.set_title(plotType+' '+title)
		
			if metric in ['b0','rE0','tE0','tE','rO',]:
				pyplot.axhline(y=0.,c='k',ls='--')
			if metric in ['b1','rgam','tgam','rR','tR',]:
				pyplot.axhline(y=1.,c='k',ls='--')
							
			#else:	ax.set_xticklabels([])
	
			for m in dmsmodels:
			    if plotStyle=='Lines':
				pyplot.plot([],[],label=modellongname[m],c=modelcolours[m],lw=2)#s = 40,lw=0.marker=modelmarkers[m],)			    
			    else:
				pyplot.scatter([],[],marker=modelmarkers[m],label=modellongname[m],c=modelcolours[m],s = 40,lw=0.)
		pyplot.subplots_adjust(bottom=0.20)

		legend = pyplot.legend(loc='lower center', ncol=5, borderaxespad=0., numpoints = 1, scatterpoints=1, prop={'size':8},) 	
		legend.draw_frame(False) 
		legend.get_frame().set_alpha(0.) 
		
	

		print 'saving',fn
		pyplot.savefig(fn,dpi=300,)

def main():
	kys=[]
	kys.extend([o + 'Months' for o in hemispheres])
	#kys.extend([o + 'Months' for o in oceans])
	#kys.extend([o + 'Seasons' for o in oceans])
	#kys.extend([o + 'Seasons' for o in oceans])
	kys.extend(['months', 'oceans',])#'OceansSeasons','OceanMonths','hemispheres','seasons',
	for xkeys in kys: #[,]:
	    for dmsmodels in ['dmsmetrics', 'dmspmetrics',]: #
		run(dmsmodels, xkeys,)

if __name__=="__main__":
	main()		
		
	print 'The end.'




