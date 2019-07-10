#!/usr/bin/ipython

import os
import numpy as np
from matplotlib import pyplot
from shelve import open as shopen


 
"""
Plot for the CMIP6 implementation paper (section describing how we selected historical initial conditions)..

This page lists the dates of the historical ICS: https://code.metoffice.gov.uk/trac/ukcmip6/ticket/98 they span 1960 to 2815.

 
Global & Southern Ocean mean SST
6 yaer running mean (annual, runing)

LHS: 1960 to 2960 (renormalized so that 1960 = 1850)
RHS: 2550 to 2850

Mark Historical ICs

I’d suggest a plot covering 1960 to 2960 of u-aw310 (with all dates renormalized
 so that 1960 = 1850). And plot both global and S. Ocean SST (with time meaning 
 tested as above) and somehow include an indication of the years of the historical ICs.

 

Then in addition, using basically the same figure, include a blow-up of the period 
2550 to 2850 to highlight the 4 ICs for historical members r16-r19 as these 4 were 
deliberately chosen to sample the centennial timescale variability in the S.Ocean.

 

I am not sure whether we can get away with just using global SST or just S. Ocean
SST or both, so if we can look at both first and then decide, along with a 
suitable time meaning that would be great.

 
"""
job_dicts = {
	'r1i1p1f2': ['u-bc179', '2250', ],
	'r2i1p1f2': ['u-bc292', '2165', ],
	'r3i1p1f2': ['u-bc370', '2120', ],
	'r4i1p1f2': ['u-bb075', '1960', ],
	'r5i1p1f3': ['u-az513', '2020', ],
	'r6i1p1f3': ['u-az515', '2050', ],
	'r7i1p1f3': ['u-az524', '1995', ],
	'r8i1p1f2': ['u-bb277', '2395', ],
	'r9i1p1f2': ['u-bc470', '2285', ],
	'r10i1p1f2': ['u-bd288', '2340', ],
	'r11i1p1f2': ['u-bd416', '2460', ],
	'r12i1p1f2': ['u-bd483', '2200', ],
	'r13i1p1f2': ['u-bf935', '2565', ],
	'r14i1p1f2': ['u-bh100', '2685', ],
	'r15i1p1f2': ['u-bh101', '2745', ],
	'r16i1p1f2': ['u-bf647', '2629', ],
	'r17i1p1f2': ['u-bf656', '2716', ],
	'r18i1p1f2': ['u-bf703', '2760', ],
	'r19i1p1f2': ['u-bh162', '2815', ],}



def movingaverage_DT(data, times, window_len=10.,window_units='years'):
        window_units = window_units.lower()
        if window_units not in ['days','months','years']:
                raise ValueError("movingaverage_DT: window_units not recognised"+str(window_units))

        data = np.ma.array(data)
        times= np.ma.array(times)

        if len(data) != len(times):
                raise ValueError("movingaverage_DT: Data and times are different lengths.")

        #####
        # Assuming time 
        if window_units in ['years',]:  window = float(window_len)/2.
        if window_units in ['months',]: window = float(window_len)/(2.*12.)
        if window_units in ['days',]:   window = float(window_len)/(2.*365.25)

        output = []#np.ma.zeros(data.shape)
        for i,t in enumerate(times):

                tmin = t-window
                tmax = t+window
                arr = np.ma.masked_where((times < tmin) + (times > tmax), data)

                #print [i,t],[tmin,tmax],[t,data[i]], arr.mean(), 'mask:',arr.mask.sum()
                output.append(arr.mean())

        return np.array(output)
        

def get_data(j, field='AMOC'):
        index = ('regionless', 'layerless', 'metricless')

	if field == 'AMOC':
        	fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/'+j+'/'+j+'_AMOC_26N.shelve'
        if field == 'Drake':
                fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/'+j+'/'+j+'_DrakePassageTransport.shelve'
        if field == 'GVT':
                fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/'+j+'/'+j+'_GlobalMeanTemperature.shelve'
        if field == 'GMT':
                fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/'+j+'/'+j+'_Temperature.shelve'
	        index = ('Global', 'Surface', 'mean')
        if field == 'SOMT':
                fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/'+j+'/'+j+'_Temperature.shelve'
	        index = ('SouthernOcean', 'Surface', 'mean')
        if field == 'AirSeaFluxCO2':
                fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/'+j+'/'+j+'_TotalAirSeaFluxCO2.shelve'

	if not os.path.exists(fn):
		print 'Does not exist:', fn	
        shelve = shopen(fn)
	data = shelve['modeldata'][index]
        shelve.close() 
	return data
	

def getClosestPoint(year, times, data):
	times = [int(t) for t in times]
	return data[times.index(year)]
	
	

def fig1():
	data1 = get_data('u-aw310', field='GMT', window_len = 10)

        times1 = sorted(data1.keys())
        data1 = [data1[t] for t in times1]

	
	#pyplot.plot(times, amoc,'k',lw=0.3)
	newd1 = movingaverage_DT(data1, times1,window_len=window_len)
	years = [yr[1] for cmip, jr in job_dicts.items()]

        pyplot.plot(times1, newd1,'k',lw=1.5)	
	
	for yr in years:
		value = getClosestPoint(year, times, newd1)
		pyplot.plot(year, value, ms='o')
	#f,(ax,ax2) = plt.subplots(1,2,sharey=True, facecolor='w')

        
        #ax2.plot(times1, newd1,'k',lw=1.5)             
        
        ax.set_xlim(1960., 2960.)
	#ax2.set_xlim(2550., 2850.)
	
	pyplot.title('Sea Surface Temperature - 10 year moving average')
        pyplot.ylabel('Celsius')
	pyplot.savefig('custom_plots/'+field+'_fig_'+str(window_len)+'.png', dpi=300)
	pyplot.close()
	
fig1()	




