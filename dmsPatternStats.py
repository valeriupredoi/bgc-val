from matplotlib import pyplot
from calendar import month_name
from glob import glob
from shelve import open as shOpen
from bgcvaltools.StatsDiagram import rmsds
from UKESMpython import folder
from itertools import product
import os
import numpy as np

purple = [125./256., 38./256., 205./256.]

oceans = ['SouthPacificOcean',  'ArcticOcean',
	  'AntarcticOcean',   'NorthAtlanticOcean','SouthAtlanticOcean',
	  'NorthPacificOcean','IndianOcean',] 
months = [month_name[i] for i in xrange(1,13)]	  

seasons = ['JFM','AMJ','JAS','OND',]

def run(key_dmsmodels, key_xkeys,):
	
	if key_dmsmodels =='dmsmetrics': 	dmsmodels = ['dms_and','dms_ara','dms_hal','dms_sim']
	if key_dmsmodels =='dmspmetrics': 	dmsmodels = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']	
	
	xkeys = []
	if key_xkeys in ['months','oceans']:	
		shelvefold = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-iMarNet-2526/'		
	else:	shelvefold = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-xkrum-2526/'	
		
	if key_xkeys == 'months':	xkeys 	= months
	if key_xkeys == 'oceans':	xkeys	= oceans
	if key_xkeys == 'seasons':	xkeys	= seasons

	if key_xkeys == 'OceanMonths':		
		xkeys	= sorted([i[0]+i[1] for i in product(oceans,months)])

	if key_xkeys == 'OceanSeasons':		
		xkeys	= sorted([i[0]+i[1] for i in product(oceans,seasons)])
		
		
	if key_xkeys not in ['months', 'OceanMonths'] and key_xkeys.lower().find('months')>0:
		i = key_xkeys.lower().find('months')	
		xkeys = [key_xkeys[:i]+' '+m for m in months]

	if key_xkeys not in ['seasons', 'OceanSeasons'] and key_xkeys.lower().find('seasons')>0:
		i = key_xkeys.lower().find('seasons')	
		xkeys = [key_xkeys[:i]+' '+s for s in seasons]
		
	modellongname = {m:i for m,i in zip(dmsmodels,['Anderson','Aranami','Halloran','Simo'])}
	modelmarkers = {m:i for m,i in zip(dmsmodels,['s','<','o','D'])}
	modelcolours = {m:i for m,i in zip(dmsmodels,['r','g','b',purple])}

	bias = {}
	stds = {}
	rcor = {}
	corr = {}
	numb = {}
	beta0 = {}
	beta1 = {}
	

	metricslist =  { 'rgam':'Scale Ratio (Data/Reference)',
			 'rE':'Normalised Difference Scale',
			 'rE0':'Normalised Bias',
			 'rR':'Spearman Correlation',
			 'rP':'p-value (robust)',
			 #'rsig':'sign',
			 'tgam':'STD Ratio (Data/Reference)',
			 'tE':'Normalised Unbiased RMSD',
			 'tE0':'Normalised Bias',
			 'tR':'Pearson Correlation',
			 'tP':'p-value (Taylor)',
			 #'tsig':'sign',			 
			 'N':'Number',
			 'b0':'Intersect',
			 'b1':'Slope',
			}
						
	metrics = {me:{} for me in metricslist.keys()}

		
	i = 0

	for o in xkeys:
		for me in metricslist.keys(): metrics[me][o] = {}

		for m in dmsmodels:
			gfn = shelvefold+m+'Surface/'+m+'_'+o.replace(' ','')+'*.shelve'
			fn = glob(gfn)
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
			sh.close()
			i+=1			
	title = key_dmsmodels
	if 'dms_p_and' in dmsmodels: 
		title += 'DMS (pixels)'
		imagefold = folder('images/MEDUSA-xkrum/DMS_Patterns/pixels')
	if 'dms_and' in dmsmodels: 
		title += 'DMS (extrapolated)'	
		imagefold = folder('images/MEDUSA-xkrum/DMS_Patterns/extrapolated')		

	#if 'SouthPacificOcean' in xkeys: title += ' Oceans'
	#if 'January' in xkeys: 		 title += ' Months'
	
	


				
	plotTypes = {	'Robust':['rE','rE0','rR',],
			'Taylor':['tE','tE0','tR',],
			#'Robust v Taylor Correlation':['rR','tR',],
			#'Robust v Taylor Correlation and N':['rR','tR','N',],
			#'Robust v Taylor Gamma':['rgam','tgam',],
			#'Robust v Taylor E':['rE','tE',],
			#'Robust v Taylor E0':['rE0','tE0',],
			'LinearRegression':['b1','b0','tR',],
			}
	
	for plotType,metrickeys in plotTypes.items():
		fn = imagefold+key_dmsmodels+'_'+ key_xkeys+'_'+plotType+'.png'
		#if os.path.exists(fn):continue
				
		fig = pyplot.figure()
		for d,metric in enumerate(metrickeys):
			ax = fig.add_subplot(len(metrickeys),1,d+1)
			xticks = []
			xticklabels=[]
			
			for o,slice in enumerate(xkeys):
				xticks.append(o)
				xticklabels.append(slice.replace('Ocean','').replace(' ','\n').replace('North','North ').replace('South','South '))
				for m, val in sorted(metrics[metric][slice].items()):
					if np.ma.is_masked(val): continue
					pyplot.scatter(o,val,marker=modelmarkers[m],c=modelcolours[m], s = 40,lw=0.)
			
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
				pyplot.scatter([],[],marker=modelmarkers[m],label=modellongname[m],c=modelcolours[m],s = 40,lw=0.)
		pyplot.subplots_adjust(bottom=0.20)

		legend = pyplot.legend(loc='lower right', ncol=4, borderaxespad=0., numpoints = 1, scatterpoints=1, prop={'size':8},) 		
		legend.draw_frame(False) 
		legend.get_frame().set_alpha(0.) 
		
	

		print 'saving',fn
		pyplot.savefig(fn,dpi=300,)

def main():
	#kys = [o + 'Months' for o in oceans]
	#kys.extend([o + 'Seasons' for o in oceans])
	kys = [o + 'Seasons' for o in oceans]

	kys.extend(['seasons','OceansSeasons',])# 'OceanMonths','months', 'oceans',
	for xkeys in kys: #[,]:
	    for dmsmodels in ['dmspmetrics','dmsmetrics', ]:
		run(dmsmodels, xkeys,)

if __name__=="__main__":
	main()		
		
	print 'The end.'




