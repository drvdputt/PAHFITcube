#!/usr/bin/env python3
from pahfit.helpers import initialize_model, fit_spectrum
from pahfit.scripts.run_pahfit import initialize_parser
from itertools import product
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from multiprocess import Pool
from pathlib import Path
from astropy.wcs import WCS
import numpy as np


def main():
    parser = initialize_parser()
    print("Hack warning: spectrumfile needs to be a cube for this script.")
    # Need to figure out the common parts between this set of arguments
    # and those set in initialize_parser, so that the single spectrum
    # and the cube case can each have there own arguments, without
    # unnecessary duplication.
    parser.add_argument("-j", type=int, default=1, help="Number of parallel processes")
    args = parser.parse_args()
    cube_spaxel_infos, map_info = read_cube(args.spectrumfile)

    num_fits = len(cube_spaxel_infos)
    if args.j > 1:
        with Pool(args.j) as p:
            # do it like this as long as memory for all the pmodels is not
            # an issue. If it becomes an issue, look at parallel iterators
            # obsfits = p.starmap(
            #     fit_spaxel,
            #     [(cube_spaxel_infos[i], args) for i in range(len(cube_spaxel_infos))],
            #     1,
            # )
            parallel_iterator = p.imap(
                fit_spaxel_wrapper,
                ((cube_spaxel_infos[i], args) for i in range(num_fits)),
            )
            obsfits = []
            for i, obsfit in enumerate(parallel_iterator):
                obsfits.append(obsfit)
                print(f"Finished fit {i}/{num_fits}")
    else:
        obsfits = [fit_spaxel(s, args) for s in cube_spaxel_infos]

    # make sure the spaxel infos and the obsfits are in the same order:
    maps_dict = initialize_maps_dict(obsfits[0], shape=(map_info["ny"], map_info["nx"]))
    for spaxel_info, obsfit in zip(cube_spaxel_infos, obsfits):
        x = spaxel_info["x"]
        y = spaxel_info["y"]
        for component in obsfit:
            for i, value in enumerate(component.parameters):
                key = feature_name(component, i)
                maps_dict[key][y, x] = value

    # save the maps
    header = map_info["wcs"].to_header()
    new_hdul = fits.HDUList()
    for key in maps_dict:
        hdu = fits.ImageHDU(data=maps_dict[key], header=header, name=key)
        new_hdul.append(hdu)

    basename = args.spectrumfile.split(".")[0]
    output_dir = Path(".") / basename
    output_dir.mkdir(exist_ok=True)
    filename = str(output_dir / basename) + "_parameter_maps.fits"
    new_hdul.writeto(filename, overwrite=True)


def initialize_maps_dict(obsfit, shape):
    """Initialize every output map using np.zeros

    Parameters
    ----------

    obsfit : fit result for one of the pixels, to provide the feature names

    shape : shape of the array used to represent the maps

    Returns
    -------

    maps_dict : dictionary mapping feature names to
"""
    maps_dict = {}
    for component in obsfit:
        for i in range(len(component.parameters)):
            key = feature_name(component, i)
            maps_dict[key] = np.zeros(shape)
    return maps_dict


def feature_name(component, param_index):
    """Consistent naming scheme for features (= parameter of component)"""
    return f"{component.name}_{component.param_names[param_index]}"


def fit_spaxel_wrapper(x):
    return fit_spaxel(x[0], x[1])


def fit_spaxel(spaxel_info, args):
    """
    Fits a single spaxel, write out to separate files

    Parameters
    ----------
    spaxel_info : dictionary where x = pixel x coordinate, y = pixel y coordinate and obsdata dict

    args : the command line arguments

    Returns
    -------
    pmodel : PAHFITBase model
        PAHFIT model
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
    basename = args.spectrumfile.split(".")[0]
    output_dir = Path(".") / basename
    output_dir.mkdir(exist_ok=True)
    outputname = str(output_dir / basename) + f"_x{x}y{y}"
    pmodel.save(obsfit, outputname, args.saveoutput)
    return obsfit


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
    cube_2dwcs = WCS(cubefile, naxis=2)

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

    map_info = {"nx": nx, "ny": ny, "wcs": cube_2dwcs}
    return spaxel_infos, map_info


if __name__ == "__main__":
    main()
