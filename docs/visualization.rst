Working with PAHFITcube results
===============================

Default visualizations
----------------------

A default vizualization is offered via ``CubeModel.maps.plot_map`` which simply
applies ``matplotlib.pyplot.imshow`` to the map data.

More complex vizualizations can be made by accessing the data stored in
``CubeModel.models`` and ``CubeModel.maps``. The former is a dictionary of
``pahfit.Model`` objects, the latter is a class called ``MapCollection``, in
which the maps can be accessed using their names as dictionary keys.
Below are a few examples

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

The maps.fits file can be loaded in DS9 by using the menu bar > file > open as >
Multiple Extension Cube.

Exporting the results
---------------------

If you would like to analyze your results in another tool, the results can be
exported as either a multi-extension fits file containing all the maps
(``MapCollection.save()``), or as a long table containing one row per pixel and
one column per map (``MapCollection.save_as_table(file_name.hdf5)``).

The former is compatible with DS9 (if loaded as multi-extension cube). The
latter can be easily loaded into Glueviz to quickly make scatter plots and
inspect relations between groups of spaxels.
