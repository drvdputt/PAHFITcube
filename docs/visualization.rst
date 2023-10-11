Working with PAHFITcube results
===============================

Default visualizations
----------------------

A default vizualization is offered via ``CubeModel.maps.plot_map`` .

More complex vizualizations can be made by accessing the data stored in
``CubeModel.models`` and ``CubeModel.maps``. The former is a dictionary of
``pahfit.Model`` objects, the latter is a class called ``MapCollection``, in
which the maps can be accessed using their names as dictionary keys. Examples to
be put here.

Exporting the results
---------------------

If you would like to analyze your results in another tool, the results can be
exported as either a multi-extension fits file containing all the maps
(``MapCollection.save()``), or as a long table containing one row per pixel and
one column per map (``MapCollection.save_as_table(file_name.hdf5)``).

The former is compatible with DS9 (if loaded as multi-extension cube). The
latter can be easily loaded into Glueviz to quickly make scatter plots and
inspect relations between groups of spaxels.
