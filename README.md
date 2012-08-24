Null tests library for Planck
=============================

This code provides a library and example scripts to run
null tests on the products of a Planck Data eXchange (DX).

Overview
--------

This repository provides a library for interactive and serial usage:

* `differences.py`

and an example script for LFI DX9 that can run either in parallel or serial:

* `run_lfidx9.py` 

Once the tests have been executed, the user can produce a HTML report
using
[plancknull_generate_html](https://github.com/ziotom78/plancknull_generate_html),
which is a separate program.

Library
-------

Input maps are read by a class that inherits from `reader.BaseMapReader` that implements a `__call__` method that given channel, frequency, survey, healfrings and requested polarization returns a map (typically either "I" or "IQU").
Therefore the first step with a new dataset is to implement a new `reader.BaseMapReader` child class.

The `differences.py` library provides 3 functions that create specific null tests, they read the maps thanks to the input reader object:

 * `halfrings`: halfring differences
 * `surveydiff`: survey differences
 * `chdiff`: either channels or horn differences

Those functions can be used interactively, see their docstrings and the example script for reference.

Serial usage
------------

In order to run all null tests serially, just set `paral=False` and run `python run_lfidx9_ipy.py`, this might take few hours, as about 380 maps are produced, the outputs will be stored in the `dx9/` folder.

Parallelization
---------------

The code is trivially parallel, parallelization is achieved using `ipython` cluster computing features.

See [`ipython` documentation][http://ipython.org/ipython-doc/stable/install/install.html#dependencies-for-ipython-parallel-parallel-computing] for installing it, it requires to build `libzmq` and install `pyzmq`.

Once `ipython` is installed, you can just:
 * run `ipcontroller` on the login node
 * submit a pbs job which spawns a number of `ipengines` (see `ipython` documentation)
 * open an `ipython` session on the login node (possibly inside `screen`) and then launch `run run_lfidx9_ipy.py` after setting `paral=True` in the code.

Known issues
------------

Currently logging does not work when the script is launched in parallel.
