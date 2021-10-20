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
import matplotlib as mpl
from matplotlib import pyplot as plt


def main():
    parser = initialize_parser()
    # Need to figure out the common parts between this set of arguments
    # and those set in initialize_parser, so that the single spectrum
    # and the cube case can each have there own arguments, without
    # unnecessary duplication.
    parser.add_argument("-j", type=int, default=1, help="Number of parallel processes")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Load pmodel for pixel from file if present",
    )
    args = parser.parse_args()
    cube_spaxel_infos, map_info = read_cube(args.spectrumfile)

    # run everything
    num_fits = len(cube_spaxel_infos)
    if args.j > 1:
        with Pool(args.j) as p:
            parallel_iterator = p.imap(
                fit_spaxel_wrapper,
                ((cube_spaxel_infos[i], args) for i in range(num_fits)),
            )
            pmodels = []
            for i, pmodel in enumerate(parallel_iterator):
                pmodels.append(pmodel)
                print(f"Finished fit {i}/{num_fits}")
    else:
        pmodels = [fit_spaxel(s, args) for s in cube_spaxel_infos]

    # make sure the spaxel infos and the obsfits are in the same order:
    maps_dict = initialize_maps_dict(pmodels[0], shape=(map_info["ny"], map_info["nx"]))
    for spaxel_info, pmodel in zip(cube_spaxel_infos, pmodels):
        x = spaxel_info["x"]
        y = spaxel_info["y"]
        for component in pmodel.model:
            for i, value in enumerate(component.parameters):
                key = feature_name(component, i)
                # initialize_maps_dict has already chosen which
                # parameters we are going to write out
                if key in maps_dict:
                    maps_dict[key][y, x] = value

    # save the maps
    header = map_info["wcs"].to_header()
    new_hdul = fits.HDUList()
    for key in maps_dict:
        hdu = fits.ImageHDU(data=maps_dict[key], header=header, name=key)
        new_hdul.append(hdu)

    basename = args.spectrumfile.split(".")[0]
    output_dir = Path(".") / (basename + "_output_per_pixel")
    output_dir.mkdir(exist_ok=True)
    outputname = str(output_dir / basename) + "_parameter_maps.fits"
    new_hdul.writeto(outputname, overwrite=True)

    # plot everything
    fontsize = 18
    font = {"size": fontsize}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=2)
    mpl.rc("axes", linewidth=2)
    mpl.rc("xtick.major", size=5, width=1)
    mpl.rc("ytick.major", size=5, width=1)
    mpl.rc("xtick.minor", size=3, width=1)
    mpl.rc("ytick.minor", size=3, width=1)
    for spaxel_info, pmodel in zip(cube_spaxel_infos, pmodels):
        obsdata = spaxel_info["obsdata"]
        obsfit = pmodel.model
        x = spaxel_info["x"]
        y = spaxel_info["y"]

        fig, axs = plt.subplots(
            ncols=1,
            nrows=2,
            figsize=(15, 10),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True,
        )

        try:
            pmodel.plot(
                axs,
                obsdata["x"],
                obsdata["y"],
                obsdata["unc"],
                obsfit,
                scalefac_resid=args.scalefac_resid,
            )
        except ValueError as error:
            print(error)
            print(f"Skipping plot x{x}y{y} due to the above error")
            continue

        # use the whitespace better
        fig.subplots_adjust(hspace=0)
        outputname = str(output_dir / basename) + f"_x{x}y{y}"
        fig.savefig("{}.{}".format(outputname, args.savefig))


def initialize_maps_dict(pmodel, shape):
    """Initialize every output map using np.zeros

    Parameters
    ----------

    pmodel : fit result for one of the pixels, to provide the feature names

    shape : shape of the array used to represent the maps

    Returns
    -------

    maps_dict : dictionary mapping feature names to
"""
    maps_dict = {}
    for component in pmodel.model:
        for i, name in enumerate(component.param_names):
            # only write out non-fixed parameters
            if not component.fixed[name]:
                key = feature_name(component, i)
                maps_dict[key] = np.zeros(shape)

    return maps_dict


def feature_name(component, param_index):
    """Consistent naming scheme for features (= parameter of component)"""
    # key will be <component name>_<constant names and values>_<fitted parameter name>
    key = component.name
    # for value, name in zip(component.parameters, component.param_names):
    #     if component.fixed[name]:
    #         key += f"_{name}{value}"
    key += f"_{component.param_names[param_index]}"
    return key


def fit_spaxel_wrapper(x):
    return fit_spaxel(*x)


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
    # get pixel indices and observed data
    x = spaxel_info["x"]
    y = spaxel_info["y"]
    obsdata = spaxel_info["obsdata"]

    # determine file name for this pixel
    basename = args.spectrumfile.split(".")[0]
    output_dir = Path(".") / (basename + "_output_per_pixel")
    output_dir.mkdir(exist_ok=True)
    outputformat = str(output_dir / (basename + f"_x{x}y{y}"))
    outputpath_full = output_dir / (basename + f"_x{x}y{y}_output." + args.saveoutput)

    if args.resume and outputpath_full.exists():
        # load model from previous (possibly partially completed) run
        pmodel = initialize_model(str(outputpath_full), obsdata, estimate_start=False)
    else:
        # setup the base model. Later, fitting could be optimized by being
        # smarter, and using information about neighbouring spaxels instead
        # of starting from scratch each time.
        pmodel = initialize_model(args.packfile, obsdata, not args.no_starting_estimate)
        obsfit = fit_spectrum(obsdata, pmodel, maxiter=args.fit_maxiter)
        # Save each pixel to separate file. Useful for "--continue" option.

        pmodel.save(obsfit, outputformat, args.saveoutput)

    return pmodel


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
    # The wavelengths can be stored in a weird format sometimes.
    # Flattening should make it consistent.
    cube_wavs = wavtable["WAVELENGTH"].quantity.flatten()

    if cube_data.shape != cube_unc_data.shape:
        print("Uncertainties cube has wrong shape!")
        exit()

    if len(cube_wavs) != cube_data.shape[0]:
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
