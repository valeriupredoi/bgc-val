from matplotlib import pyplot
from calendar import month_name
from glob import glob
from shelve import open as shOpen
from bgcvaltools.StatsDiagram import rmsds

purple = [125./256., 38./256., 205./256.]

def run(key_dmsmodels, key_xkeys,robust = False):
	
	if key_dmsmodels =='dmsmetrics': 	dmsmodels = ['dms_and','dms_ara','dms_hal','dms_sim']
	if key_dmsmodels =='dmspmetrics': 	dmsmodels = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim']	
	
	if key_xkeys == 'months':
		xkeys = [month_name[i] for i in xrange(1,13)]
	if key_xkeys == 'oceans':		
		xkeys	=['SouthPacificOcean',  'ArcticOcean',
			  'AntarcticOcean',   'NorthAtlanticOcean','SouthAtlanticOcean',
		 	  'NorthPacificOcean','IndianOcean',] 

	
	shelvefold = '/data/euryale7/scratch/ledm/ukesm_postProcessed/MEDUSA-iMarNet-2526/'
	

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
	i = 0

	for o in xkeys:
		bias[o] = {}
		stds[o] = {}
		corr[o] = {}
		numb[o] = {}
		rcor[o] = {}
		beta0[o] = {}
		beta1[o] = {}
		for m in dmsmodels:
			gfn = shelvefold+m+'Surface/'+m+'_'+o+'*.shelve'
			fn = glob(gfn)
			#print gfn
			if len(fn)==1:fn = fn[0]
			print i, o,m,'opening',fn
			sh = shOpen(fn)
			#print sh.keys()
			#assert False
			if robust:
				gam = sh['robust.gamma']
				E0  = sh['robust.E0'] 
			else:
				gam = sh['Taylor.gamma']
				E0  = sh['Taylor.E0'] 

			corr[o][m] = sh['Taylor.R']			
			rcor[o][m] = sh['robust.R']
			numb[o][m] = sh['N']
			beta0[o][m] = sh['b0']
			beta1[o][m] = sh['b1']	
			
			sig= gam>1 and 1 or -1
			E=rmsds(gam,R)
			
			bias[o][m] = sig*E
			stds[o][m] = E0

			
			i+=1
			sh.close()
	title = ''
	if 'dms_p_and' in dmsmodels: title += 'DMS (pixels)'
	if 'dms_and' in dmsmodels: title += 'DMS (extrapolated)'	

	if 'SouthPacificOcean' in xkeys: title += ' Oceans'
	if 'January' in xkeys: 		 title += ' Months'
	
	numbers = False
	if numbers:	subplots = [411,412,413,414]
	else:		subplots = [311,312,313]
	
	fig = pyplot.figure()
	for sbp,data in zip(subplots,[bias,stds,corr, numb]):
		ax = fig.add_subplot(sbp)
		xticks = []
		xticklabels=[]
		for o,slice in enumerate(xkeys):
			xticks.append(o)
			xticklabels.append(slice.replace('Ocean','').replace('North','North\n').replace('South','South\n'))
			for m, val in sorted(data[slice].items()):
				#print sbp, m,modellongname[m],o,val,modelmarkers[m],modelcolours[m]
				pyplot.scatter(o,val,marker=modelmarkers[m],c=modelcolours[m], s = 40,lw=0.)
		if data == bias: pyplot.ylabel('Bias')
		if data == stds: pyplot.ylabel('Std')
		if data == corr: pyplot.ylabel('Corr')
		if data == numb: pyplot.ylabel('Number')
	
		print xticks,xticklabels
		if sbp in [414,313,]: 	pyplot.xticks(xticks, xticklabels,rotation = 'vertical',)
		else:			pyplot.xticks(xticks, ['',])
		if sbp in [411,311]:	ax.set_title(title)
		
		#else:	ax.set_xticklabels([])
	
		for m in dmsmodels:
			pyplot.scatter([],[],marker=modelmarkers[m],label=modellongname[m],c=modelcolours[m],s = 40,lw=0.)
	pyplot.subplots_adjust(bottom=0.20)

#	pyplot.legend()
	legend = pyplot.legend(loc='lower right', ncol=4, borderaxespad=0., numpoints = 1, scatterpoints=1, prop={'size':8},) 		
		
	legend.draw_frame(False) 
	legend.get_frame().set_alpha(0.) 
		
	

	
	

	
	fn = 'images/MEDUSA-iMarNet/DMS_Patterns/'+key_dmsmodels+'_'+ key_xkeys+'.png'
	if robust:	fn = fn.replace('.png','_robust.png')
	print 'saving',fn
	pyplot.savefig(fn,dpi=300,)
	#pyplot.show()

def main():
	for xkeys in ['months', 'oceans']:
	    for dmsmodels in ['dmsmetrics', 'dmspmetrics',]:
	      for robust in (1,0):
		run(dmsmodels, xkeys,robust = robust)

if __name__=="__main__":
	main()		
		
	print 'The end.'




