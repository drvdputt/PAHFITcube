import numpy as np
from astropy.table import Table
from pahfit.base import PAHFITBase
import os
import pkg_resources


def make_trimmed_model(packfile, obsdata):
    """
    Load a complete science model, and remove any components that are outside of the wavelength range.

    Parameters
    ----------

    packfile:
        Path to the ipac file from which the complete model will be loaded.

    obsdata:
        dict containing wavelengths 'x', flux 'y', and uncertainty 'unc'

    Returns
    -------
    pmodel: PAHFITBase model
        PAHFIT model based on trimmed science pack table

    """
    # If packfile not in current directory, try to find one of the same
    # name that comes with PAHFIT
    if not os.path.isfile(packfile):
        pack_path = pkg_resources.resource_filename("pahfit", "packs/")
        test_packfile = "{}/{}".format(pack_path, packfile)
        if os.path.isfile(test_packfile):
            packfile = test_packfile
        else:
            raise ValueError("Input packfile {} not found".format(packfile))

    # determine wavelength range
    w = obsdata["x"].value
    wmin = np.amin(w)
    wmax = np.amax(w)

    # Different idea: trim the astropy table first, then continue as normal
    t = Table.read(packfile, format="ipac")

    # decide which rows we are going to keep
    keep_row = np.full(len(t), True)

    is_drude_or_gauss = np.logical_or(t["Form"] == "Drude1D", t["Form"] == "Gaussian1D")
    # Only keep drudes and gauss with center within range. It is called
    # x0 for both so we can use the same code.
    feature_x0 = t[is_drude_or_gauss]["x_0"]
    keep_row[is_drude_or_gauss] = np.logical_and(wmin < feature_x0, feature_x0 < wmax)

    # now parse the trimmed table
    param_info = PAHFITBase.parse_table(t[keep_row])

    # and create a new model (and don't pass a file name, so that the
    # current contents of param_info are used)
    trimmed_model = PAHFITBase(
        obsdata["x"].value,
        obsdata["y"].value,
        estimate_start=True,
        param_info=param_info,
    )
    return trimmed_model
