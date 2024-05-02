PAHFITcube
===========

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

A collection of tools for applying PAHFIT (https://github.com/PAHFIT/pahfit) to IFU data cubes,
and producing feature maps from the fit results.

Usage
-----

.. note::
   The guide below is somewhat outdated. I am in the process of updating the
   sections below, and moving them to the appropriate pages on the new
   `readthedocs site for PAHFITCube <https://pahfitcube.readthedocs.io>`_.

Using the ``-j <nprocs>`` option will run the fits in parallel with the given number of
processes. Currently there’s also a basic ``--resume`` option, which skips pixels when their
output file is already present. It loads the fit results for those pixels from the output files
instead. To use PAHFITcube, the following things are needed.

1. Input data in the right format. The currently supported data is any data cube that can be
   loaded as a Spectrum1D object (see specutils). If you have a cube that works with
   ``Spectrum1D.read()``, you are good to go. Otherwise, some tips are given in the "Preparing
   Input" section below.
2. Settings for PAHFIT that will work for the features and resolution of the given data, i.e.
   choosing the right instrument pack and science pack.
3. A way to tie everything together. We suggest either:

   a. Your own notebook or script that prepares these data and then uses the Python interface
      (``pahfitcube.cube_model``)
   b. The right input files and options for the command line script (``run_pahfitcube.py``)

Python interface
^^^^^^^^^^^^^^^^

The main interface for the user is the ``CubeModel`` class. To set up a cube model, a regular
PAHFIT model needs to be passed. This model will be used as the initial condition for the
fitting.

.. code-block::

   from pahfitcube.cube_model import CubeModel
   from pahfit.model import Model

   initial_model = Model.from_saved("pahfit_output.ecsv")
   cube_model = CubeModel(initial_model)

We recommend performing a PAHFIT fit to e.g. the average SED over the cube, to obtain a suitable
initial condition. As a reminder, PAHFIT models can be saved using
``model.save("pahfit_output.ecsv")``. See the `PAHFIT documentation
<https://pahfit.readthedocs.io/en/latest/fit_spectrum.html#fitting/>`_ for instruction on
fitting.

To perform the fit, the cube data need to be fed to the CubeModel. The wavelengths need to have
the right units, as they will be converted to micron under the hood. Additionally some metadata
needs to be set so that PAHFIT knows which instrument model to use.

.. code-block::

   # code that prepares your data
   # ...
   # package it in Spectrum1D, and set the right metadata so that it is supported by PAHFIT
   from specutils import Spectrum1D
   from astropy.nddata import StdDevUncertainty
   spec = Spectrum1D(spectral_axis=..., flux=..., uncertainty=StdDevUncertainty(...))
   spec.meta["instrument"] = "jwst.nirspec.g395.high"

The fitting can then be performed by the call below. It is recommended to use the
``checkpoint_prefix`` option, so that each fitted pixel is written to disk. This allows a cube
model to resume fitting if interrupted, or to read a (possibly only partially fitted) model from
storage.

.. code-block::

   # command that will start the fitting
   cube_model = fit(spec, './output_dir/output_prefix', maxiter=10000, j=7)

Multi-processing is supported with the ``j`` option, and can provide a large speedup considering
that a typical PAHFIT job takes anywhere between 0.1 and 100 seconds depending on the complexity
of the model. Currently, all pixels are fit independently, but the cube model could be sped up
further by developing "smarter" fitting algorithms making use of cross-pixel information.

Once the fitting has completed, the fit information will be stored in several places

1. The fitted PAHFIT models for each pixel ``(x, y)`` can be accessed as ``cube_model.models[(x, y)]``
2. If a directory was provided, the model for each pixel was saved to './output_dir/output_prefix_xy_...'
3. The fit results can be written to a big table with one row per pixel, and one column per
   parameter, using ``cube_model.maps.save_as_table("big_table.ecsv", format="ascii.ecsv")``.
   This is ideal for pixel-vs-pixel analysis, e.g. scatter plots comparing feature strength
   correlations.

When these per-pixel results are obtained during the ``fit()`` call or while reading from disk,
the fit parameters are collected into maps. These maps can be accessed via ``cube_model.maps``,
which is an instance of ``MapCollection``. The latter offers a transparent way to deal with
many maps on the same spatial grid, and several utilities for inspecting or saving them.

.. code-block::

   map_collection = cube_model.maps

   # get a map array for a parameter
   map_data = map_collection['PAH_15.9_powerz']

   # plot a single map
   map_collection.plot_map('PAH_15.9_power')

   # plot overview of many maps
   map_collection.plot_map_collage(['H2_O(3)_2.9_power', 'H2_O(4)_power', 'H2_O(5)_power'])

   # make a WCS, and save as multi-extention fits file, which can be displayed in e.g. DS9
   wcs = astropy.wcs.WCS(...)
   map_collection.save(wcs, "maps.fits")


Script
^^^^^^
The script ``run_pahfit_cube.py`` has the same command line arguments as the normal PAHFIT run
script, except a data cube is given at the input instead of a single spectrum.

*Things have changed recently and this script needs to be redesigned*

Other Utitlies
--------------

Preparing input
^^^^^^^^^^^^^^^

If multiple instruments are combined to achieve a larger wavelength coverage, the cubes should
be merged which means
1. Reprojecting everything onto the same spatial grid
2. Fixing any jumps in the flux between the different parts
3. Put the data in the right input format (package it in a Spectrum1D object).

Reprojection
,,,,,,,,,,,,

Because the field of view is often different for each IFU cube of an observation, a generic
script that reprojects data cubes is provided: ``merge_cubes.py``. It uses the ``reproject``
package, and writes out the result as a Python pickle, which contains the necessary ingredients
to build a Spectrum1D object. These pickles can be loaded in by the run script.

Spectral Order Stitching
,,,,,,,,,,,,,,,,,,,,,,,,

Spectroscopic observations typically have jumps in the flux, between spectral orders. These need
to be fixed before giving the spectrum to PAHFIT.

A generic stitching implementation might be included in the future.

License
-------

This project is Copyright (c) Dries Van De Putte and licensed under the terms of the GNU GPL v3+
license. This package is based upon the `Astropy package template
<https://github.com/astropy/package-template>`_ which is licensed under the BSD 3-clause
license. See the licenses folder for more information.


Contributing
------------

We love contributions! pahfitcube is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
pahfitcube based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
