#!/usr/bin/env python3
from itertools import product
from astropy.io import fits
from multiprocess import Pool
from pathlib import Path
from astropy.wcs import WCS
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from specutils import Spectrum1D
import pickle
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from dataclasses import dataclass

from pahfit.helpers import initialize_model, fit_spectrum
from pahfit.scripts.run_pahfit import initialize_parser
from pahfit.base import PAHFITBase

from pahfit_cube.deprecated.make_trimmed_model import initialize_trimmed_model




@dataclass
class Spaxel:
    """Properties of a single spaxel, pass around frequently.

    x and y: pixel coordinates, to remember which spectrum goes where

    obsdata: dict in the same format as the output of pahfit.helpers.read_spectrum:
        {'x': wavelengths, 'y': flux, 'unc': uncertainty}.
    """

    x: int
    y: int
    obsdata: dict


def make_output_path(spectrumfile, suffix):
    """Creates output path for file

    New file will be <working directory>/<spectrumfile base name>/<spectrumfile base name + suffix>

    spectrumfile: the file on which the suggested file name will be based

    suffix: output name will be the stem of spectrumfile + suffix.

    """
    path = Path(spectrumfile)
    output_dir = Path(".") / (path.stem + "_output_per_pixel")
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / f"{path.stem}_{suffix}"
    return output_path


def run_pahfit_cube(args):
    """Main function, called after parsing the args.

    Can also be run by importing it, and passing the right args. Run
    ./run_pahfit_cube.py --help to see which arguments can be specified.
    The ones currently used

    """
    spaxels, map_info = read_cube(args.spectrumfile)
    # When joint fitting is supported in PAHFIT, allow a list of cubes
    # to be passed, instead of just a single cube. PAHFIT will then fit
    # the cubes (with different wavelength ranges) simultaneously. They
    # will all need to have the same WCS and spatial coverage though.

    # run everything
    num_fits = len(spaxels)
    if args.j > 1:
        with Pool(args.j) as p:
            parallel_iterator = p.imap(
                fit_spaxel_wrapper,
                ((spaxels[i], args) for i in range(num_fits)),
            )
            models = []
            for i, model in enumerate(parallel_iterator):
                models.append(model)
                print(f"Finished fit {i}/{num_fits}")
    else:
        models = [fit_spaxel(s, args) for s in spaxels]

    # make sure the spaxel infos and the obsfits are in the same order:
    maps_dict = initialize_maps_dict(models[0], shape=(map_info["ny"], map_info["nx"]))
    for spaxel, model in zip(spaxels, models):
        for component in model:
            for i, value in enumerate(component.parameters):
                key = feature_name(component, i)
                # initialize_maps_dict has already chosen which
                # parameters we are going to write out
                if key in maps_dict:
                    maps_dict[key][spaxel.y, spaxel.x] = value

    # save the maps
    header = map_info["wcs"].to_header()
    new_hdul = fits.HDUList()
    for key in maps_dict:
        hdu = fits.ImageHDU(data=maps_dict[key], header=header, name=key)
        new_hdul.append(hdu)
        outputname = make_output_path(args.spectrumfile, "parameter_maps.fits")
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
    for spaxel, model in zip(spaxels, models):
        plot_spaxel_result(args, spaxel, model)


def initialize_maps_dict(model, shape):
    """Initialize every output map using np.zeros

    Parameters
    ----------

    model : fit result for one of the pixels, to provide the feature names

    shape : shape of the array used to represent the maps

    Returns
    -------

    maps_dict : dictionary containing feature names as keys, and zeros
    arrays to work with"""
    maps_dict = {}
    for component in model:
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


def fit_spaxel(spaxel, args):
    """
    Fits a single spaxel, write out to separate files

    Parameters
    ----------
    spaxel : Spaxel

    args : the command line arguments

    Returns
    -------
    pmodel : PAHFITBase model
        PAHFIT model
    """
    # determine file name for this pixel
    suffix = f"x{spaxel.x}y{spaxel.y}"
    # pahfit adds stuff to the above format depending on file type
    suffix_ipac = f"{suffix}_output.{args.saveoutput}"
    output_path_format = make_output_path(args.spectrumfile, suffix)
    output_path_ipac = make_output_path(args.spectrumfile, suffix_ipac)

    if args.resume and output_path_ipac.exists():
        # load model from previous (possibly partially completed) run
        pmodel = initialize_model(
            str(output_path_ipac), spaxel.obsdata, estimate_start=False
        )
        print(f"Loaded existing fit results from {output_path_ipac}")
    else:
        # setup the base model. Later, fitting could be optimized by being
        # smarter, and using information about neighbouring spaxels instead
        # of starting from scratch each time.
        # pmodel = initialize_model(
        #     args.packfile, spaxel.obsdata, not args.no_starting_estimate
        # )
        pmodel = initialize_trimmed_model(args.packfile, spaxel.obsdata)
        obsfit = fit_spectrum(spaxel.obsdata, pmodel, maxiter=args.fit_maxiter)
        # Save each pixel to separate file. Useful for "--continue" option.
        pmodel.save(obsfit, str(output_path_format), args.saveoutput)

    return obsfit


def plot_spaxel_result(args, spaxel, model):
    obsdata = spaxel.obsdata
    x = spaxel.x
    y = spaxel.y

    fig, axs = plt.subplots(
        ncols=1,
        nrows=2,
        figsize=(15, 10),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    try:
        PAHFITBase.plot(
            axs,
            obsdata["x"],
            obsdata["y"],
            obsdata["unc"],
            model,
            scalefac_resid=args.scalefac_resid,
        )
    except ValueError as error:
        print(error)
        print(f"Skipping plot x{x}y{y} due to the above error")
        # raise error
        # continue

    # use the whitespace better
    fig.subplots_adjust(hspace=0)

    # save using meaningful file name (contains coordinates + figure file type)
    output_path = make_output_path(args.spectrumfile, f"x{x}y{y}.{args.savefig}")
    fig.savefig(output_path)
    plt.close(fig)


def read_cube(cubefile, only_one=False):
    """
    Read multiple spectra from a data cube and convert input units to
    the expected internal PAHFIT units.

    Spaxels for which the flux is entirely zero will not be loaded.

    Parameters
    ----------
    cubefile : string
        File with the spectral cube to be fit. The cube should be in the
        format desribed in <link needed>. The cube file (of the same
        shape) containing the uncertainties should be "name_unc.fits".

    Returns
    -------
    spaxels : list of Spaxel
    """
    if cubefile.endswith(".fits"):
        # Try to load in this way. Confirmed to work for JWST cubes.
        # Other cubes might need 'format' argument.
        spec = Spectrum1D.read(cubefile)
        # This does not seem to load the spatial wcs properly, so load
        # that in a different way
        cube_2dwcs = WCS(cubefile, naxis=2)
    elif cubefile.endswith(".pickle"):
        # We're dealing with output from the merging script here. See
        # wcshacks.write_merged_cube.
        with open(cubefile, "rb") as f:
            data_dict = pickle.load(f)
        try:
            flux = data_dict["flux"]
            spectral_axis = data_dict["spectral_axis"]
            # mock the uncertainty to 10% for now
            uncertainty = StdDevUncertainty(0.1 * flux)
            spec = Spectrum1D(flux, spectral_axis, uncertainty=uncertainty)
            cube_2dwcs = data_dict["spatial_wcs"]
        except Exception as e:
            print(
                "Something wrong with pickle. Needs to come from wcshacks.write_merged_cube."
            )
            raise e

    # problem: pahfit assumes Jy, not Jy / sr. Should ask about this.
    cube_qty = spec.flux.to(u.MJy / u.sr)
    cube_unc_qty = spec.uncertainty.quantity.to(u.MJy / u.sr)
    cube_wavs = spec.spectral_axis.to(u.micron)
    ny, nx, _ = spec.shape

    spaxels = []
    # option to use only center pixel, useful for testing
    xy_tuples = [(nx // 2, ny // 2)] if only_one else product(range(nx), range(ny))
    for x, y in xy_tuples:
        obsdata = {"x": cube_wavs, "y": cube_qty[y, x], "unc": cube_unc_qty[y, x]}
        if not all(obsdata["y"] == 0):
            spaxels.append(Spaxel(x, y, obsdata))

    print(f"{len(spaxels)}/{nx * ny} spaxels will be fit")

    map_info = {"nx": nx, "ny": ny, "wcs": cube_2dwcs}
    return spaxels, map_info

def main(args_list=None):
    """
    args_list: list of str
        override the command line arguments (useful for calling main
        from another script)
    """
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
    if args_list is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args_list)

    run_pahfit_cube(args)


if __name__ == "__main__":
    main()
