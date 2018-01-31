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

	if year0=='piControl':  
               for t in sorted(mdata.keys()):
                        if jobID == 'u-ar766': t1 = t + (2594 -1850)    #-1869
			else: t1 = t 
                        times.append(float(t1))
                        datas.append(mdata[t])  
               return times, datas

				
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
                        if jobID == 'u-ar538': t1 = t - (2061-56) + (3996-2106)      #-1862
	                if jobID == 'u-ao586': t1 = t - (2561 - 556)+ (3996-2106)    #-1869
                        if jobID  in ['u-ar783',]: t1 = t + (4965. - 2108.)

                        print jobID, t1,t, t0 
                        times.append(float(t1))
                	datas.append(mdata[t])
               return times, datas

        if year0=='Staggered':
               for t in sorted(mdata.keys()):
                        if   jobID in ['u-au362', ]: t1 = t + 40. 
                        elif jobID in ['u-au364', ]: t1 = t + 45.
                        elif jobID in ['u-au365', ]: t1 = t + 50.
			else: t1 = t
                        if t1<2215.:continue

                        times.append(float(t1))
                        datas.append(mdata[t])
               return times, datas

	if year0 in ['OOvFC','OOvFC1','OOvFC2',]:
               for t in sorted(mdata.keys()):
                        if jobID == 'u-ak900': t1 = t - 2106
                        if jobID == 'u-an869': t1 = t - (2061-56) + (3996-2106)      #-1862
                        if jobID == 'u-ar538': t1 = t - (2061-56) + (3996-2106)      #-1862
                        if jobID == 'u-ao586': t1 = t - (2561 - 556)+ (3996-2106)    #-1869
                        if jobID  in ['u-ar783',]: t1 = t + (4965. - 2108.)
			if year0 in ['OOvFC','OOvFC2',]:
				if t1 < 4940: continue 
        	                if t1 > 5110: continue
        	        if year0 in ['OOvFC1',]:
                                if t1 < 2390: continue
                                if t1 > 2860: continue

                        print jobID, t1,t, t0
                        times.append(float(t1))
                        datas.append(mdata[t])
               return times, datas

	if year0=='Drift':
               for t in sorted(mdata.keys()):
#                        if jobID in ['u-as462','u-as643',]:	t1 = t  +3312.
			if jobID == 'u-ar977': t1 =t
			else:      t1 = t  +3312.
			if t1 < 5400. : continue
                        if t1 > 5600. : continue

                        times.append(float(t1))
                        datas.append(mdata[t])
               return times, datas

        if year0=='Drift2':
               for t in sorted(mdata.keys()):
#                        if jobID in ['u-as462','u-as643',]:    t1 = t  +3312.
                        if jobID == 'u-ar977': t1 =t
                        else:      t1 = t  +3312.
                        if t1 < 5400. : continue
                        if t1 > 5550. : continue
                        times.append(float(t1))
                        datas.append(mdata[t])
               return times, datas

        if year0=='Drift3':
               for t in sorted(mdata.keys()):
#                        if jobID in ['u-as462','u-as643',]:    t1 = t  +3312.
                        if jobID == 'u-ar977': t1 =t
                        else:      t1 = t  +3312.
                        if t1 < 5405. : continue
                        if t1 > 5565. : continue
			if jobID == 'u-as462':
				if t1 > 5438: continue

                        if jobID == 'u-as858':
                                if t1 > 5485: continue

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

        if year0=='UKESM_0.9.1':
               for t in sorted(mdata.keys()):
                        if jobID == 'u-ar379': t1 = t + 458.
			else: t1=t
                        times.append(float(t1))
                        datas.append(mdata[t])
               return times, datas

        if year0=='UKESM_0.9.2':
               for t in sorted(mdata.keys()):
   #                     if jobID == 'u-ar379': t1 = t + 458.
   #                     else: t1=t
			
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

        if year0=='ignoreStart':
               t0 = float(sorted(mdata.keys())[0])
               for t in sorted(mdata.keys()):
			if float(t) < 4600.:continue
                        times.append(float(t))
                        datas.append(mdata[t])
               return times, datas

       	if year0=='2000-2250':
               for t in sorted(mdata.keys()):       
			if float(t) <2000.:continue
		        if float(t) >2250.:continue
		        times.append(float(t))
		        datas.append(mdata[t])
               return times, datas		        
        if year0=='4830-5000':
               for t in sorted(mdata.keys()):
                        if float(t) <4830.:continue
                        if float(t) >5000.:continue
                        times.append(float(t))
                        datas.append(mdata[t])
               return times, datas
        if year0 in ['4945-5110','4945-5110i','4945-5110ii','4945-5110iii','4800-5100']:
		if jobID in ['u-an869','u-ar783','u-ar538','u-ar766','u-at629','u-at793','u-at628','u-at760','u-as462','u-as858', 'u-at572', 'u-au027']:
	               	for t in sorted(mdata.keys()):
                                if jobID in ['u-at629','u-at793','u-at628','u-at760','u-as462','u-as858','u-at572', 'u-au027']:  t1 = t + 2838.
				#if jobID in ['u-at629','u-at793','u-at628','u-at760',]:#'u-as462','u-as858',]:
				#				t1 = t + 2794.
                                #elif jobID in ['u-as462','u-as858',]:    t1 = t  +2858.
				elif jobID  in ['u-ar783',]: 	t1 = t + 4965. - 2108.
                                elif jobID  in ['u-ar766',]: 	t1 = t + 3095.
				else:	t1 = t
	
				if year0 in ['4945-5110','4945-5110i','4945-5110ii','4945-5110iii',]:
	       		                if float(t1) <4945.: continue
                		        if float(t1) >5410.: continue
                                if year0 in ['4800-5100',]:
                                        if float(t1) <4780. : continue
                                        if float(t1) >5410. : continue

	                        if jobID == 'u-as462':
					if t > 2126.: continue
	                        if jobID == 'u-as858':
                                        if t > 2175.: continue

        	                #        if t1 > 5485 -3312. + 2794.: continue

                        	times.append(float(t1))
	                        datas.append(mdata[t])

                if jobID in ['u-am064','u-am927','u-aq853','u-am064i','u-am927i','u-aq853i','u-am064ii','u-am927ii','u-aq853ii','u-am064iii','u-am927iii','u-aq853iii',]:
                        for t in sorted(mdata.keys()):
                                if jobID in ['u-am064i','u-am927i','u-aq853i',]:     t1 = t+2835     # 2130-2230
                                if jobID in ['u-am064ii','u-am927ii','u-aq853ii',]:      t1 = t+2675     # 2290-2390
                                if jobID in ['u-am064iii','u-am927iii','u-aq853iii',]:     t1 = t+2515     # 2450-2550

				#if year0 == '4945-5110i':	 t1 = t+2835	# 2130-2230
                                #if year0 == '4945-5110ii':      t1 = t+2675	# 2290-2390
                                #if year0 == '4945-5110iii':     t1 = t+2515	# 2450-2550
				if jobID=='u-am927iii':
					if t1 > 4994: continue
                                #if year0 in ['4945-5110','4945-5110i','4945-5110ii','4945-5110iii',]:
                                if float(t1) <4945.:continue
                                if float(t1) >5110.:continue

                                times.append(float(t1))
                                datas.append(mdata[t])
               	return times, datas

        if year0 =='from2228':
		buff= 0
                for t in sorted(mdata.keys()):
                        if float(t) <2228. - buff :continue
                        times.append(float(t))
                        datas.append(mdata[t])
                return times, datas
        if year0 =='from2265':
                buff= 0
                for t in sorted(mdata.keys()):
                        if float(t) <2265. - buff :continue
                        times.append(float(t))
                        datas.append(mdata[t])
                return times, datas
        if year0 =='from4950':
		buff = 0
                for t in sorted(mdata.keys()):
                        if float(t) <4950.-buff:continue
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

        if year0 =='2500-3000':
               for t in sorted(mdata.keys()):
                        if float(t) <2500.:continue
                        if float(t) >3000.:continue
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
				
