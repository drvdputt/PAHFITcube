Fitting
=======

Inputs and workflow
-------------------


Using the run script
--------------------

See::

  pahfitcube/scripts/run_pahfitcube.py --help

In the Python runtime
---------------------

The explanation here could be quite elaborate. Make sure you understant the
relevant PAHFIT concepts first! For now TL;DR::

  from pahfit.model import Model
  from pahfitcube.cube_model import CubeModel
  
  spec = <your spectrum1D cube here>
  spec['instrument'] = 'your pahfit instrument string'
  cm = CubeModel(Model.from_saved('pahfit_result.ecsv'))
  cm.fit(spec, checkpoint_prefix='output/model', j=1, maxiter=10000)

During the fitting, the results are stored in ``cm.models``, and on disk at the
given prefix (one PAHFIT model per spaxel).
