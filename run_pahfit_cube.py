#!/usr/bin/env python3
from pahfit.helpers import initialize_model, fit_spectrum
from pahfit.scripts.run_pahfit import initialize_parser
from itertools import product, repeat
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from multiprocessing import Pool


def main():
    parser = initialize_parser()
    print("Hack warning: spectrumfile needs to be a cube for this script.")
    # Need to figure out the common parts between this set of arguments
    # and those set in initialize_parser, so that the single spectrum
    # and the cube case can each have there own arguments, without
    # unnecessary duplication.
    args = parser.parse_args()
    cube_spaxel_infos = read_cube(args.spectrumfile)
    with Pool(7) as p:
        p.starmap(
            fit_spaxel,
            [(cube_spaxel_infos[i], args) for i in range(len(cube_spaxel_infos))],
            1,
        )


def fit_spaxel(spaxel_info, args):
    """
    Fits a single spaxel, write out to separate files
    """
    # get pixel indices (tentative)
    x = spaxel_info["x"]
    y = spaxel_info["y"]
    # get obsdata for this spaxel
    obsdata = spaxel_info["obsdata"]
    # setup the base model. Later, fitting could be optimized by being
    # smarter, and using information about neighbouring spaxels instead
    # of starting from scratch each time.
    pmodel = initialize_model(args.packfile, obsdata, not args.no_starting_estimate)
    obsfit = fit_spectrum(obsdata, pmodel, maxiter=args.fit_maxiter)

    # for now, we will save the results to separate files. But
    # later, we should not have a file per pixel, with all features,
    # but a file per feature, with all pixels.
    outputname = args.spectrumfile.split(".")[0] + f"_x{x}y{y}"
    pmodel.save(obsfit, outputname, args.saveoutput)


def read_cube(cubefile):
    """
    Read multiple spectra from a data cube and convert input units to
    the expected internal PAHFIT units.

    Parameters
    ----------
    cubefile : string
        File with the spectral cube to be fit. The cube should be in the
        format desribed in <link needed>. The cube file (of the same
        shape) containing the uncertainties should be "name_unc.fits".

    Returns
    -------
    spaxel_infos : list of dict
        One entry per pixel. 'x' is x coordinate of pixel (in pixel
        space), 'y' is y coordinate of pixel (in pixel space), 'obsdata'
        is a dict in the same format as the output of read_spectrum.
    """
    cube_hdu = fits.open(cubefile)["PRIMARY"]
    cube_unit = u.Unit(cube_hdu.header["BUNIT"])
    cube_data = cube_hdu.data
    cube_qty = cube_data * cube_unit

    cube_unc_data = fits.getdata(cubefile.replace(".fits", "_unc.fits"))
    cube_unc_qty = cube_unc_data * cube_unit

    wavtable = Table.read(cubefile)
    cube_wavs = wavtable["WAVELENGTH"].quantity

    if cube_data.shape != cube_unc_data.shape:
        print("Uncertainties cube has wrong shape!")
        exit()

    if len(wavtable) != cube_data.shape[0]:
        print("Wavelength table and cube are not compatible")
        exit()

    nwavs, ny, nx = cube_data.shape

    spaxel_infos = []
    for x, y in product(range(nx), range(ny)):
        obsdata = {"x": cube_wavs, "y": cube_qty[:, y, x], "unc": cube_unc_qty[:, y, x]}
        spaxel_infos.append({"x": x, "y": y, "obsdata": obsdata})
    return spaxel_infos


if __name__ == "__main__":
    main()
