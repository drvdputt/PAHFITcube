Fitting
=======

Inputs and workflow
-------------------

To use PAHFITcube, the following things are needed.

1. Input data in the right format. The currently supported data is any data cube
   that can be loaded as a Spectrum1D object (see specutils). If you have a cube
   that works with ``Spectrum1D.read()``, you are good to go. Otherwise, some
   tips are given in the `input preparation article <cubes.rst>`.
2. Settings for PAHFIT that will work for the features and resolution of the
   given data, i.e. choosing the right instrument pack and science pack.

Then, everything needs to be tied together with a script or notebook, as
explained below.

Using the run script
--------------------

The script ``run_pahfit_cube.py`` has the same command line arguments as the normal PAHFIT run
script, except a data cube is given at the input instead of a single spectrum.

Using the ``-j <nprocs>`` option will run the fits in parallel with the given number of
processes. Currently there’s also a basic ``--resume`` option, which skips pixels when their
output file is already present. It loads the fit results for those pixels from the output files
instead.

See also::

  pahfitcube/scripts/run_pahfitcube.py --help


In the Python runtime
---------------------

.. note::

   To better understand the terms below, it is recommended to familiarize
   oneself with the `PAHFIT <https://pahfit.readthedocs.io>` concepts first.

The basic code to perform a PAHFITcube run from within Python will look
something like this.

.. code::
  from pahfit.model import Model
  from pahfitcube.cube_model import CubeModel

  spec = <your spectrum1D cube here>
  spec['instrument'] = 'your pahfit instrument string'
  initial_model = Model.from_saved("pahfit_output.ecsv")
  cm = CubeModel(initial_model)
  cm.fit(spec, checkpoint_prefix='output/model', j=1, maxiter=10000)

The main interface for the user is the ``CubeModel`` class. To set up a cube model, a regular
PAHFIT model needs to be passed. This model will be used as the initial condition for the
fitting.

We recommend performing a PAHFIT fit to e.g. the average SED over the cube, to obtain a suitable
initial condition. As a reminder, PAHFIT models can be saved using
``model.save("pahfit_output.ecsv")``. See the `PAHFIT documentation
<https://pahfit.readthedocs.io/en/latest/fit_spectrum.html#fitting/>`_ for instructions on
fitting.

To perform the fits over an entire data cube, the cube data need to be passed to
the fit function the the CubeModel. The wavelengths need to have the right
units, as they will be converted to micron under the hood. Additionally some
metadata needs to be set so that PAHFIT knows which instrument model to use.
Both PAHFIT and PAHFITcube use the Spectrum1D class as the input container, and
the spectrum with instrument model metadata can be set up analogously to `the
way it is done for PAHFIT
<https://pahfit.readthedocs.io/en/latest/fit_spectrum.html#python-notebook-and-scripts>`_. An example would be

.. code::

   # code that prepares your data
   # ...
   # package it in Spectrum1D, and set the right metadata so that it is supported by PAHFIT
   from specutils import Spectrum1D
   from astropy.nddata import StdDevUncertainty
   spec = Spectrum1D(spectral_axis=..., flux=..., uncertainty=StdDevUncertainty(...))
   spec.meta["instrument"] = "jwst.nirspec.g395.high"

The fitting can then be performed by the ``cm.fit()`` call. It is recommended to
use the ``checkpoint_prefix`` option, so that each fitted pixel is written to
disk. In the above example, there will be a file for every pixel, sotred in the
``output/`` directory, under the name ``model_xy_<x>_<y>.ecsv``. This allows a cube
model to resume fitting if interrupted, or to read amodel from storage.
Multi-processing is supported with the ``j`` option, and can provide a large
speedup considering that a typical PAHFIT job takes anywhere between 0.1 and 100
seconds depending on the complexity of the model. Currently, all pixels are fit
independently, but the cube model could be sped up further by developing
"smarter" fitting algorithms making use of cross-pixel information.

After the fitting, the results can be retrieved from the CubeModel instance. For
a guide to how to use the results, see the `visualization guide
<visualization.rst>`. Briefly summarized

1. ``cm.models`` is a dictionary, where every key is an integer tuple ``(x,
y)``, and the values are regular PAHFIT ``Model`` objects. E.g. ``pahfit_8_10 =
cm.models[(8, 10)]``.
2. ``cm.maps`` is an instance of ``MapCollection``. The latter offers a
   transparent way to deal with many maps on the same spatial grid, and several
   utilities for inspecting or saving them. The results for the fit parameters
   are gathered into easily accessible 2D arrays, which can be accessed as a
   dictionary. E.g., ``cm.maps["PAH_3.3_power"]`` yields a map of the fitted
   power of the 3.3 feature.
3. The fit results can be written to a big table with one row per pixel, and one
   column per parameter, using ``cube_model.maps.save_as_table("big_table.ecsv",
   format="ascii.ecsv")``. This is ideal for pixel-vs-pixel analysis, e.g.
   scatter plots comparing feature strength correlations.
