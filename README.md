Batch run of null tests for Planck/LFI
======================================

The code in this repository runs a number of null tests on the
products of a Planck/LFI Data eXchange (DX).

Overview
--------

This repository provides two Python scripts that run the tests:

* `run.py` is the sequential version
* `run_parallel.py` is the parallelized version, to be run on a
  multi-core/multi-node machine.

Once the tests have been executed, the user can produce a HTML report
using
[plancknull_generate_html](https://github.com/ziotom78/plancknull_generate_html),
which is a separate program.

Installing and running the code
-------------------------------

You should have a fairly recent version of Python (currently it has
been tested using Python 2.7.2, but it might be runnable even with
2.6.x).

The code needs the [`scoop`](https://code.google.com/p/scoop/)
library, if you have [`pip`](http://pypi.python.org/pypi/pip/) you can
install it using either

    pip install scoop

or

    pip install --user scoop

depending whether you have admin's right or not.

To run the code, you must set the value of the following environment
variables:

* =DX9_LFI= should point to the path containing the FITS files
  relevant for the DX under study
* =NULLTESTS_ENV= specifies the naming convention used for the FITS
  files. It must be either `LFIDPC` or `NERSC`.

As an example, under Bash you can run the code with this command:

    DX9_LFI=/foo/bar NULLTESTS_ENV=LFIDPC python run.py
