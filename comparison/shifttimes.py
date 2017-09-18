#
# Copyright 2015, Plymouth Marine Laboratory
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

def shifttimes(mdata, jobID,year0=False):


	times,datas=[],[]
	try:	t0 = float(sorted(mdata.keys())[0])
	except:	return times, datas       
				
	if year0=='juggling2':	
               for t in sorted(mdata.keys()):
                        if jobID == 'u-an869': t1 = t - (2061-56)       #-1862
                        if jobID == 'u-ao586': t1 = t - (2561 - 556)    #-1869
	                if t1<540:continue
        	        if t1>700:continue
                	times.append(float(t1))
                        datas.append(mdata[t])	
               return times, datas
               
      	if year0=='FullSpinUp':
               for t in sorted(mdata.keys()):
        		if jobID == 'u-ak900': t1 = t - 2106
        	        if jobID == 'u-an869': t1 = t - (2061-56) + (3996-2106)      #-1862
	                if jobID == 'u-ao586': t1 = t - (2561 - 556)+ (3996-2106)    #-1869
                        print jobID, t1,t, [t0, tm1]
                        times.append(float(t1))
                	datas.append(mdata[t])
               return times, datas
               
	if year0=='UKESM0.8':
               for t in sorted(mdata.keys()):
                        if jobID == 'u-am004': t1 = t -1978 - 203 
                        if jobID == 'u-ao365': t1 = t -1978 - 33  
                        if jobID == 'u-ao837': t1 = t -1978
			if t1<0:continue
			
                        print jobID, t1,t  
                        times.append(float(t1))
                        datas.append(mdata[t])
               return times, datas

        if year0=='u-ai567-minus3':  
               t0 = float(sorted(mdata.keys())[0])
               if jobID =='u-ai567': t0 = t0+3
               for t in sorted(mdata.keys()):
                        times.append(float(t)-t0)
                        datas.append(mdata[t])
               return times, datas
                                

       	if year0=='2000-2250':
               for t in sorted(mdata.keys()):       
			if float(t) <2000.:continue
		        if float(t) >2250.:continue
		        times.append(float(t))
		        datas.append(mdata[t])
               return times, datas		        
               
	if year0 =='2000-2600normu-ak900':
               for t in sorted(mdata.keys()):       	
			if jobID == 'u-ak900': t = t - 1937.
		        if float(t) <2000.:continue
		        if float(t) >2600.:continue
		        times.append(float(t))
		        datas.append(mdata[t])
               return times, datas		        
            
        #####
        # Set all jobs to start at time zero.
        if year0 in ['True', True,'First100Years',]:
               for t in sorted(mdata.keys()):
                        if year0=='First100Years' and float(t) - t0 > 100.:continue
                        times.append(float(t)-t0)
                        datas.append(mdata[t])
               return times, datas
               
	######
	# No year shift requested, returning sorted arrays.	
	times 	= sorted(mdata.keys())
	datas	= [mdata[t] for t in times]
        return times, datas				
				
