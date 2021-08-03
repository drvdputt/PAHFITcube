from pahfit.helpers import read_spectrum, initialize_model, fit_spectrum
from pahfit.scripts.run_pahfit import initialize_parser
from itertools import product
import numpy as np


def main():
    parser = initialize_parser()
    print("Hack warning: spectrumfile needs to be a cube for this script.")
    # Need to figure out the common parts between this set of arguments
    # and those set in initialize_parser, so that the single spectrum
    # and the cube case can each have there own arguments, without
    # unnecessary duplication.

    args = parser.parse_args()

    cube_spaxel_infos = read_cube(args.spectrumfile)

    # setup the base model. The fitting for each spaxel will start from a
    # copy of this. Later, fitting could be optimized by being smarter,
    # and using information about neighbouring spaxels instead of
    # starting from scratch each time.
    pmodel_original = initialize_model(
        args.packfile, obsdata, not args.no_starting_estimate
    )

    for spaxel_info in cube_spaxel_infos:
        # get pixel indices (tentative)
        x = spaxel_info["x"]
        y = spaxel_info["y"]

        # get obsdata equivalent for this cube
        obsdata = spaxel_info["obsdata"]
        pmodel = pmodel_original.copy()
        obsfit = fit_spectrum(obsdata, pmodel, maxiter=args.fit_maxiter)

        # for now, we will save the results to separate files. But
        # later, we should not have a file per pixel, with all features,
        # but a file per feature, with all pixels.
        outputname = args.spectrumfile.split(".")[0] + "_x{}y{}"
        pmodel.save(obsfit, outputname, args.saveoutput)


def fit_spaxel():
    pass


def read_cube(cubefile):
    """
    Read multiple spectra from a data cube and convert input units to
    the expected internal PAHFIT units.

    Parameters
    ----------
    cubefile : string
        File with the spectral cube to be fit. The cube should be in the
        format desribed in <link needed>.

    Returns
    -------
    spaxel_infos : list of dict
        One entry per pixel. 'x' is x coordinate of pixel (in pixel
        space), 'y' is y coordinate of pixel (in pixel space), 'obsdata'
        is a dict in the same format as the output of read_spectrum.
    """
    # loop over all spaxels of the cube. Still need to write out the cube
    # we will use here in resample_and_merge_spitzer_cubes.py.
    spaxel_infos = []
    nx = 5
    ny = 5
    nwavs = 80
    cube_data = np.zeros((nwavs, ny, nx))
    cube_unc = np.zeros((nwavs, ny, nx))
    cube_wavs = np.zeros(nwavs)
    for x, y in product(range(nx), range(ny)):
        obsdata = {"x": cube_wavs, "y": cube_data[:, y, x], "unc": cube_unc[:, y, x]}
        spaxel_infos.append({"x": x, "y": y, "obsdata": obsdata})
    return spaxel_infos
