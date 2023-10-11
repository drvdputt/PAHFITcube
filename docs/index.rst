##########
PAHFITcube
##########

``pahfitcube``, stylized as PAHFITcube, is a python package that provides tools
to apply `PAHFIT <https://pahfit.readthedocs.io/en/latest/>` to spectroscopic
data cubes.

The main goals are to provide...

- ... a way to straightforwardly bookkeep a large number of PAHFIT fits,
  typically one fit per spatial element of a datacube. The main fitting method
  can run using multiprocessing, and can easily be restarted and resumed in case
  of interruption.
- ... tools to inspect and analyze the results. This includes creating maps and
  tables from the many individual fits, even if the fits are only partially
  complete.
- ... tools to assist with preparing the input data. This includes reprojecting
  the cubes corresponding to individual spectral orders and merging the
  wavelength axis, as well as ways to save and load these intermediate products.

Main Articles
=============

The goals listed above are explained in the following articles

.. toctree::
   :maxdepth: 2

   Preparing input data <cubes.rst>
   Fitting a cube <fitting.rst>
   Vizualizing and exporting <vizualization.rst>

Installation
============

Currently, only installation from source is supported. Since all dependencies
are specified in the package configuration, the installation can still proceed

with pip.

If you are only going to use ``pahfitcube``, and not make modifications, a
simple pip command should suffice.::

    pip install git+https://github.com/drvdputt/PAHFITcube

However, since this tool is in its early stages, I expect that you might want to
edit some of the code. Therefore, I recommend cloning the repository and
doing an editable pip install.::

    git clone https://github.com/drvdputt/PAHFITcube
    cd PAHFITcube
    pip install -e .


Getting Started
===============

WORK IN PROGRESS. Eventually, there should be a Jupyter notebook included in the
repository with a basic setup.
