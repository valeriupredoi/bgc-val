



from paths import orcaGridfn
from netCDF4 import Dataset
import numpy as np
import UKESMpython as ukp
from matplotlib import pyplot
import cartopy.crs as ccrs
from bgcvaltools.pftnames import getLongName


regionList	= [#'Global', 'ignoreInlandSeas',
  		'SouthernOcean','Remainder',
		'Equator10', 
		'NorthernSubpolarAtlantic','NorthernSubpolarPacific','ArcticOcean',
		]

def robinPlotCustom(lons,lats,data,filename,title, zrange=[-100,100],drawCbar=True,cbarlabel='',doLog=False,dpi=100,cmapname='default',crude=True):
	####
	# Based on robinplotSingle
	
	fig = pyplot.figure()
	fig.set_size_inches(10,5)

	lons = np.array(lons)
	lats = np.array(lats)
	data = np.ma.array(data)
	
	rbmi = min([data.min(),])
	rbma = max([data.max(),])
	
	if rbmi * rbma >0. and rbma/rbmi > 100.: doLog=True

	print lons.shape,lats.shape,data.shape
	lon0 = 0.#lons.mean()
	if crude:
		ax = pyplot.subplot(111)#,projection=ccrs.PlateCarree(central_longitude=lon0, ))
		im = pyplot.scatter(lats,lons,c=data,lw=0.,s=3,cmap='viridis',vmin=rbmi,vmax=rbma,)
		#pyplot.colorbar(im)
		#title, zrange=[rbmi,rbma],lon0=lon0,drawCbar=False,cbarlabel=cbarlabel,doLog=doLog,cmap = cmapname)	
	else:
		ax = pyplot.subplot(111,projection=ccrs.PlateCarree(central_longitude=lon0, ))
		fig,ax,im = ukp.makemapplot(fig,ax,lons,lats,data,title, zrange=[rbmi,rbma],lon0=lon0,drawCbar=False,cbarlabel=cbarlabel,doLog=doLog,cmap = cmapname)


          
          
	cmap = im.get_cmap()
	#for i in [0,10,100,1000,10000,1000000,100000]:
	#	print i, cmap(i), data.min(),data.max()
		
	for i,r in enumerate(regionList):
		c = cmap((i)/5.)
		print i,r,c
		pyplot.hist([],color=c,label=getLongName(r))

	# Shrink current axis's height by 10% on the bottom
	box = ax.get_position()
	ax.set_position([box.x0, box.y0 + box.height * 0.1,
		         box.width, box.height * 0.9])

	# Put a legend below current axis
	leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
		  ncol=3,prop={'size':9})
		  	
	#leg = pyplot.legend(loc='lower center',ncol=3, )
	leg.draw_frame(False) 
	leg.get_frame().set_alpha(0.)		
	
	print "robinPlotSingle.py:\tSaving:" , filename
	pyplot.savefig(filename ,dpi=dpi)		
	pyplot.close()





def makeRegionMap():

	plotAll = 0#True	# make plots for all regions
	imageFold = ukp.folder('images/maps')
	#####
	# Load data.
	nc = Dataset(orcaGridfn,'r')
	bathy = nc.variables['mbathy'][:]
	xy = np.ma.masked_where(bathy==0,nc.variables['nav_lat'][:]).compressed()
	xx = np.ma.masked_where(bathy==0,nc.variables['nav_lon'][:]).compressed()
	nc.close()
	
	cbathy = np.ma.masked_where(bathy==0,bathy).compressed()
	xt = np.ones_like(cbathy)
	xz = np.ones_like(cbathy)

	####
	# Calculate masks, based on lat/lon.
	masks = {}
	for r in regionList:
		masks[r] = ~ukp.makeMask('',r, xt,xz,xy,xx,cbathy,debug=True)
	
	#####
	# Turn mask into one field.
	data = np.zeros_like(cbathy)
	for i,r in enumerate(regionList):
		data += (i+1)* masks[r]
		if plotAll:
			fn = imageFold+'Region_Legend_'+r+'.png'		
			ukp.robinPlotSingle(xy, xx, masks[r],fn,r,drawCbar=True,cbarlabel='',doLog=False,dpi=100,)
	data = np.ma.masked_where(data==0,data)
	
	#####
	# Send it to the plotting tool.
	colourmaps = ['default',]#'rainbow','jet','gist_earth','terrain','ocean','hsv','gist_rainbow','nipy_spectral',]
	for c in colourmaps:
		fn = imageFold+'Region_Legend.png'
		robinPlotCustom(xy, xx, data,fn,'',drawCbar=False,cbarlabel='',doLog=False,dpi=200,cmapname = c)
		
	

if __name__=="__main__":
	makeRegionMap()		
	
