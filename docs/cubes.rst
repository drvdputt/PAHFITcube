PAHFITcube input
================

As you will read in the article about running the `fitting <fitting.rst>`, the
input format of PAHFITcube is centered around the ``Spectrum1D`` class from the
``specutils`` package.

Easiest option: Auto-loadable fits file
---------------------------------------

In short: anything that is loadable by ``Spectrum1D.read`` can be given to the
PAHFITcube run script.

Most general option: Pass a Spectrum1D object in the Python runtime
-------------------------------------------------------------------

PAHFITcube requires that the spectra are on a uniform spatial grid. If one
wishes to fit a wavelength ranges that spans across multiple spectral segments,
then the data needs to be prepared by the user first, into a single data cube.

Best practices
,,,,,,,,,,,,,,

It is alway best to avoid multiple reprojections. For JWST specifically, we
recommend that all data cubes are build on the same spatial grid. This can be
achieved directly in the pipeline, by using the appropriate cube build options.
E.g., to build all cubes (NIRSpec + MIRI) to the same grid as Channel 4 LONG,
one can use the following options.

.. code::

   strun cube_build asn.json \
      --output_dir coarse_wcs \
      --output_type 'band' \
      --scalexy 0.35 \
      --ra_center 83.83535169747073 \
      --dec_center -5.419729392828107 \
      --cube_pa 70.42338204969806 \
      --nspax_x 79 \
      --nspax_y 15

Please consult the jwst pipeline documentation instructions on how to `run the
pipeline
<https://jwst-pipeline.readthedocs.io/en/latest/jwst/user_documentation/running_pipeline_command_line.html>`_
and `cube build
https://jwst-pipeline.readthedocs.io/en/latest/jwst/cube_build/main.html`_ in
particular, and for the format of asn.json, the association file.

Cube merge tool
,,,,,,,,,,,,,,,

If there is no "native" or "best practice" way to merge your cube data (e.g. you
have data from different observatories with incompatible toolkits), a few
general tools are provided by PAHFITcube.

Since for most instruments, IFU data are typically delived as a collection of
spectral segments, this is a task that most users will have to deal with.
Therefore, PAHFITcube provides some generic utilities which can help the user
with preparing the data. Common tasks are
- Reprojecting the segment cubes onto the same spatial grid
- Concatenating the reprojected cubes along the wavelength axis
- Saving the data so that they can be loaded later to perform a fit (instead of
  having to run the reprojection before every fit).

The module ``pahfitcube.merge_cubes`` provides tools for these tasks. Each of
the above steps can be performed individually, to make intermediate adjustments.
But to get started, I recommend using
``pahfitcube.merge_cubes.make_usable_cube``. The code to apply it would look
something like this ::

    from specutils import Spectrum1D
    from pahfitcube.merge_cubes import make_usable_cube

    your_file_names = ['segment_a.fits', 'segment_b.fits']
    list_of_cubes = [Spectrum1D.read(fn) for fn in your_file_names]
    merged_cube, celestial_wcs = make_usable_cube(list_of_cubes, min_wavelength, max_wavelength)

The assumption here is that your fits files are loadable by ``Spectrum1D``, and
that the loading process properly includes the WCS info into
``s.meta['header']``, with ``s`` the ``Spectrum1D`` instance produced by
``Spectrum1D.read()``.
