.. BGC-Val documentation master file, created by
   sphinx-quickstart on Fri Sep 23 10:33:38 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Biogeochemistry Validation toolkit
===================================

The Biogeochemistry Validation toolkit (BGC-Val) is a toolkit built to assist with the validation of marine ocean models. 
The tools were originally built for analysis of UKESM annual fields, to help identify 
potential issues during the spin up phase. 

Contents:


Accessing the toolkit
--------------

The code itself is available upon request from the 
`PML-hosted gitlab server. <https://gitlab.ecosystem-modelling.pml.ac.uk/UKESM/bgc-val>`_

Access can be requested via the `Code registration form <http://www.pml.ac.uk/Modelling_at_PML/Access_Code>`_.
The PML-contact is "Lee de Mora", please select "Other" in the Model code/Tutorial, and please put
"Access to BGC-val toolkit" as the main purpose of use.  Note that processing may take a few days.

Once access to the gitlab server has been established, and you have set up your ssh keys, the toolkit can be
downloade with the command::

	git clone git@gitlab.ecosystem-modelling.pml.ac.uk:UKESM/bgc-val.git

 



.. toctree::
   :maxdepth: 2


theWholePackage
--------------------
.. automodule:: theWholePackage
    :members: theWholePackage

analysis_timeseries
--------------------
.. automodule:: analysis_timeseries
    :members: analysis_timeseries

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

