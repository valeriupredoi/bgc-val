Lee's Instructions for Original `bgc-val`
=========================================

Once it’s all set up the script to run a minimal version of BGC-val will be:

```
./analysis_timeseries.py u-bc179 debug
```

But obviously, it needs to do a lot more than that!

## Contents:
- Installation
- download data from mass.
- timeseries
- point to point
- html report maker
- run single job whole package analysis
- run multi-job comparison analysis
- Run all of this from pml or crontab.

## Install:
    standard python packages: numpy, scipy, matplotlib, netCDF4, pyyaml, pyproj
    Harder to install packages: mpl_toolkits, cartopy, python-mpltoolkits.basemap
    netcdf_manip: a set of tools I wrote to manipuate netcdf in python. One job of this update is to remove this dependency.

        https://github.com/ledm/NetCDF_manip
        Unfortunately, this is also in python2! 

    Load jaspy python2.7 module:

        module load jaspy/2.7

    set up paths.py file:

        One job of this project is to replace this with yml instead of python.
        Do this on jasmin by copying Paths/paths_jasmin.py to paths.py

## Download:

Download data from mass. This needs to be done on the mass-cli1.jasmin.ac.uk VM. So, ssh to that machine and run the script:

In BGC-val/bgcvaltools: ./downloadFromMass.py u-bc179


However, this might be weird as you needv all your mass credentials set up. If it fails, you can use my data, which is:

    /gws/nopw/j04/ukesm/BGC_data/u-bc179/u-bc179o_1y_18621201-18631201_grid_T.nc

If you want to be really fancy, I launch this command remotely from pml using:

    ./bgc-val/RemoteScripts/pmlcron_mass.sh u-bc179

 

Note that this doesn't just download the data. It also creates symbolic links to make the format look like the original UKESM beta data. This should probably be fixed in the future – but it’s a massive ball ache as NEMO constantly changes its file naming convention. It happened so often. Plus it may be different for alternative versions of UKESM like UKESM1.1 ort UKESM-fast.


Single model time series - debug

To run a single model analysis, you run:

```./analysis_timeseries.py u-bc179 debug```

If that works, then you can expand the key from debug to kmf (key metrics first), and then level1
 

Single model point to point - debug

```./analysis_p2p.py u-bc179 2010 debug```

If that works, then you can expand the key from debug to level2

Single model html report maker

```./makeReport.py u-bc179 2010```

Single model whole package

```./theWholePackage.py u-bc179 2010```

This runs the time series analysis twice: kmf then level1, then the point to point analysis over the year 2010, then it creates a single model report.

Multi-model comparison report

The Multi-model comparison report command is:

./ analysis_compare.py

 

This part of the tool sucks! All the comparison suites are hardwired into the python – and it copies a lot of the code from analysis_timeseries.py but doesn’t use it. Nightmare. A lot of work is needed to improve it. For now, lets get the single jobID stuff working first.

 

So that’s the bulk of the python.
