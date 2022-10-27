from specutils import Spectrum1D
from astropy.wcs import WCS
from astropy.nddata import StdDevUncertainty
import pickle
from itertools import product

def read_cube(cubefile, only_one=False):
    """
    Read data cube and determine which pixels are suitable for fitting.

    E.g. Spaxels for which the flux is entirely zero will not be loaded.

    Parameters
    ----------
    cubefile : string
        File with the spectral cube to be fit. The cube should be in the
        format desribed in <link needed>. The cube file (of the same
        shape) containing the uncertainties should be "name_unc.fits".

    Returns
    -------
    spec: Spectrum1D
        The data cube

    xy_tuples: list of (x,y) pairs suitable for fitting

    wcs: the spatial wcs

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
    else:
        raise ValueError("Cube file not compatible.")

    ny, nx, _ = spec.shape
    xy_tuples = [(nx // 2, ny // 2)] if only_one else product(range(nx), range(ny))
    # ^ has option to use only center pixel, useful for testing

    return spec, xy_tuples, cube_2dwcs
