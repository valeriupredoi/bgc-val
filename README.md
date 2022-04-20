Environment and installation
============================

Grab miniconda3, install `mamba` in the `base` environment:

```
conda install -c conda-forge mamba
```

then create the `bgcval` environment and activate it:

```
mamba env create -n bgcval -f environment.yml
conda activate bgcval
```

then install the development dependencies:

```
pip install -e .[develop]
```

Python 2 to Python 3
====================

Install `2to3`:

`pip install 2to3`

Usage: use the 3.9 extension and write to disk option:

`2to3-3.9 script.py -w`

Remove the backup `.py.bak` files or stash them.
