# 1. (Easy) a time series of only UKESM1 u-aw310 for the full length (~1960 to year 3???) applying either a 5 year or 10 year running mean (likely the former).

# 2. A time series (using the same time meaning as in 1...so likely 5 year running mean) that consists of the following:
#    (i) A 250 year segment of the piControl (it probably does not matter what 250 years are chosen).
#    (ii) The ensemble mean (of however many historical runs have completed...think it is now 11), 5 year running mean, of the 165 years of the (11) historical runs
#    (iii) A 5-member ensemble mean for each of the 4 Tier 1 scenarios for the 85 years of the projections.



import numpy as np
from matplotlib import pyplot
from shelve import open as shopen



def fig1():
	fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/u-az515/u-az515_AMOC_26N.shelve'
	shelve = shopen(fn)
        print shelve.keys()
