PAHFIT-cube
===========

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

A collection of tools for applying PAHFIT (https://github.com/PAHFIT/pahfit) to IFU data cubes, with the goal of
producing feature maps.

Functionality
-------------

In this readme, we provide a quick overview of the supported input data, the main run command,
and the output. The full documentation can be found at readthedocs (not published yet).

Run PAHFIT-cube
^^^^^^^^^^^^^^^

The script ``run_pahfit_cube.py`` has the same command line arguments as the normal PAHFIT run
script, except a data cube is given at the input instead of a single spectrum. Additionally, a
suitable science pack needs to be passed. In the example below, the default Spitzer
extragalactic science pack (part of PAHFIT), is given.

::

   python run_pahfit_cube.py reprojected.fits scipack_ExGal_SpitzerIRSSLLL.ipac --fit_maxiter 50

For debugging the script, it is best to set ``--fit_maxiter`` to a low number, so that
everything goes faster. Using the ``-j <nprocs>`` option will run the fits in parallel with the
given number of processes. Currently there’s also a basic ``--resume`` option, which skips
pixels when their output file is already present. It loads the fit results for those pixels from
the output files instead.

Input
^^^^^

The script supports any spectral data cube that can be loaded as a Spectrum1D object by
specutils, and from which astropy can extract a WCS.

Reprojection
,,,,,,,,,,,,

Because the field of view is often different for each IFU cube of an observation, some scripts
to reproject data cubes are provided. These are specialized per instrument, as data formats can
be different. These scripts use the ``reproject`` package, and write out the result as a Python
pickle, which can also be loaded in by the run script.

Spectral Order Stitching
,,,,,,,,,,,,,,,,,,,,,,,,

The reprojection scripts contain a basic spectral order stitching step. For example, for Spitzer
data cubes (LL1, LL2, SL1, SL2), the spectra are rescaled on a pixel-per-pixel basis, according
to the average values in the wavelength overlap region. This is done in order of decreasing
wavelength: - scale LL2 to match LL1 - scale SL1 to match LL2 - scale SL2 to match SL1
After this step, the slices of the 4 cubes are merged and sorted by wavelength.

Output
^^^^^^

Feature maps
,,,,,,,,,,,,

The main output is a multi-extension fits file ``reprojected_parameter_maps.fits``, with the
same WCS as the input spectral cube. Each extension contains a map of one of the fit parameters.
It can be viewed in DS9 when opened as a multi-extension cube.

Per pixel files and resume option
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

For every pixel, the normal PAHFIT output is written out (fit parameters
and a plot). This is useful for diagnosing problems with individual
pixels, and determining/recalling progress.

Test data
^^^^^^^^^

The following info should be reworked, and moved to readthedocs.io.

We do our experiments and testing with Spitzer IFU cubes from the
SAGE-Spec (citation) program. From this dataset, one HII region was
taken, for which the wavelength range from 35 ~ 5 micron is covered by 4
IFU cubes: LL1, LL2, SL1, and SL2 (in the order from long to short
wavelengths).

If you are on the ERS-PDRs Slack space, you can look for the data there.
Or ask me for it. The list of files should look like this

::

   ls -1 PAHFIT-cube/data/sage-spec_hii1_hii8_4dec08
   hii1_hii8_ll_LL1_cube.fits
   hii1_hii8_ll_LL1_cube_unc.fits
   hii1_hii8_ll_LL2_cube.fits
   hii1_hii8_ll_LL2_cube_unc.fits
   hii1_hii8_ll_LL3_cube.fits
   hii1_hii8_ll_LL3_cube_unc.fits
   hii1_hii8_sl_SL1_cube.fits
   hii1_hii8_sl_SL1_cube_unc.fits
   hii1_hii8_sl_SL2_cube.fits
   hii1_hii8_sl_SL2_cube_unc.fits
   hii1_hii8_sl_SL3_cube.fits
   hii1_hii8_sl_SL3_cube_unc.fits

Because PAHFIT needs the whole wavelength range to work with, the 4
Spitzer cubes need to be merged first. The script
``merge_spitzer_cubes.py`` takes care this. Since this is an early test,
it does not have command line arguments yet, so it only works for the
data above. Simply run

::

   python merge_spitzer_cubes.py


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
