"""Prototype for cube wrapper script around PAHFIT.

For now, it will only work for the cubes of the SAGE-Spec program.

Given a set of data cubes, a set of spectra is extracted by performing
aperture photometry on each slice. The apertures used are rectangular
and adjacent in the pixel space of the map we are trying to create. The
axes of the map align with RA,DEC by construction. 

We use a manually created WCS to turn the pixel grid into a coordinate
grid with equal angular separation along each axis. The coordinates are
used to create a set sky-apertures, which can then be used to perform
aperture photometry on the cube slices, given the right WCS for each
cube.

"""
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from pathlib import Path
from matplotlib import pyplot as plt
from photutils.aperture import SkyRectangularAperture, aperture_photometry
import reproject
from collections import namedtuple
from dataclasses import dataclass


@dataclass
class Cube:
    """Fixed pattern for the information in one cube that we want to pass
    around

    """

    file_handle: fits.HDUList
    data: np.ndarray
    wavelength: np.ndarray
    wcs: WCS


@dataclass
class CubeSet:
    """The four cubes we have for each object

    Constructor is automatically generated"""

    ll1: Cube
    ll2: Cube
    sl1: Cube
    sl2: Cube

    def all_cubes(self):
        """You can add functions to this type of class like this"""
        return [self.ll1, self.ll2, self.sl1, self.sl2]


def read_spitzer_cube(fn):
    """
    Get the data cube, wavelength table and WCS.

    Returns
    -------
    Cube
    """
    hdulist = fits.open(fn)
    wavs = Table.read(fn)["WAVELENGTH"][0].flatten()
    wcs = WCS(fn, naxis=2)
    return Cube(file_handle=hdulist, data=hdulist[0].data, wavelength=wavs, wcs=wcs)


def get_SAGE_cubes(target, uncertainty=False):
    """Get all 4 Spitzer IFU cubes (LL1, LL2, SL1, SL2).

    uncertainty : if True, loads the uncertainty cubes instead

    Returns
    -------
    CubeSet
    """
    dname = f"data/sage-spec_{target}_4dec08"
    d = Path(dname)
    # load them in the right order...
    cube_names = ["ll_LL1", "ll_LL2", "sl_SL1", "sl_SL2"]
    cubes = []
    for cube_name in cube_names:
        fname = f"{target}_{cube_name}_cube.fits"
        # reuse this function
        if uncertainty:
            fname = fname.replace("cube.fits", "cube_unc.fits")
        cubes.append(read_spitzer_cube(str(d / fname)))

    # ... in the future we can unambiguously access them here
    return CubeSet(*cubes)


def quicklook_cubes(cube_set, apertures=None):
    """Display the 4 cubes and a set of apertures in a 2x2 grid."""
    fig = plt.gcf()
    for i, cube in enumerate(cube_set.all_cubes()):
        ax = fig.add_subplot(2, 2, i + 1, projection=cube.wcs)
        ax.imshow(cube.data[-1])
        ax.grid()
        ax.coords[0].set_format_unit(u.degree, decimal=True)
        ax.coords[1].set_format_unit(u.degree, decimal=True)
        if apertures is not None:
            apertures.to_pixel(cube.wcs).plot(axes=ax, color="r")


# the wcs we want our final map to have
def make_ra_dec_wcs(center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec):
    """Make simple WCS where X and Y are aligned with RA and DEC, respectively.

    center_ra, center_dec: determines crval

    pix_angle_delta: physical distance between pixels, in decimal degrees

    npix_ra, npix_dec: number of pixels along each axis. Is needed to
    make sure that (center_ra, center_dec) corresponds to the middle of
    the image. The physical image size will be pix_angle_delta * npix.

    """
    w = WCS(naxis=2)
    # center of each pixel = 1, 2, 3, ...
    # 3 pixels --> center is 2
    # 4 pixels --> center is 2.5 (border between 2 and 3)
    center_x = npix_ra / 2 + 0.5
    center_y = npix_dec / 2 + 0.5
    w.wcs.crpix = [center_x, center_y]
    w.wcs.crval = [center_ra, center_dec]
    w.wcs.cdelt = [-pix_angle_delta, pix_angle_delta]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def reproject_and_merge_cubes(
    cube_set, center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec, filename=None
):
    """
    Reproject all four cubes onto pixel grid wcs, and merge the result.

    Result is sorted by wavelength"""
    output_projection = make_ra_dec_wcs(
        center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec
    )
    output_wavs = np.concatenate([c.wavelength for c in cube_set.all_cubes()]).flatten()
    nw = len(output_wavs)

    ny, nx = npix_dec, npix_ra
    output_cube_array = np.zeros((nw, ny, nx))

    start = 0
    for c in cube_set.all_cubes():
        stop = start + len(c.wavelength)
        output_cube_array[start:stop] = reproject_cube(c, output_projection, ny, nx)
        start = stop

    # stitching here

    order = np.argsort(output_wavs)
    output_wavs = output_wavs[order]
    output_cube_array = output_cube_array[order]

    if filename is not None:
        # multi extension (wav list and cube)
        new_hdul = fits.HDUList()
        # cube as primary hdu
        header = output_projection.to_header()
        header["PC3_3"] = 1
        header["CRPIX3"] = 1
        header["CRVAL3"] = 1
        header["CTYPE3"] = "WAVE-TAB"
        header["CUNIT3"] = "um"
        header["PS3_0"] = "WCS-TAB"
        header["PS3_1"] = "WAVELENGTH"
        header["BUNIT"] = "MJy/sr"
        new_hdul.append(fits.PrimaryHDU(data=output_cube_array, header=header))
        # wavs as bintable hdu
        wav_col = fits.Column(name="WAVELENGTH", array=output_wavs, format="D")
        wavhdu = fits.BinTableHDU.from_columns([wav_col])
        wavhdu.header["EXTNAME"] = "WCS-TAB"
        wavhdu.header["TUNIT1"] = "um"
        new_hdul.append(wavhdu)
        # write
        new_hdul.writeto(filename, overwrite=True)

    return output_wavs, output_cube_array


def reproject_cube(cube, wcs, ny, nx):
    """
    Reproject every slice of cube onto wcs, ny, nx.

    Returns
    -------
    output_array: np.ndarray indexed on wavelength, y, x
    """
    num_wavs = len(cube.wavelength)
    input_array = cube.data
    input_wcs = cube.wcs

    output_array = np.zeros((num_wavs, ny, nx))
    for w in range(num_wavs):
        output_array[w], footprint = reproject.reproject_adaptive(
            input_data=(input_array[w], input_wcs),
            output_projection=wcs,
            shape_out=(ny, nx),
        )
    return output_array


def make_square_aperture_grid(
    center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec
):
    """
    Create sky apertures representing pixels of RA-DEC aligned map.

    Use sky apertures immediately (instead of starting with pixel
    apertures and then converting), to make picking the right size
    easier.
    """
    # all xy pairs
    X, Y = np.mgrid[:npix_ra, :npix_dec]
    x = X.flatten()
    y = Y.flatten()

    w = make_ra_dec_wcs(center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec)
    positions = w.pixel_to_world(x, y)
    size = pix_angle_delta * u.degree
    return SkyRectangularAperture(positions, size, size)


##################################OBSOLETE######################################
def extract_spectra(cube_dicts, sky_apertures):
    wavtables = [d["wavelength"] for d in cube_dicts]
    cubes = [d["cube"] for d in cube_dicts]
    wcses = [d["wcs"] for d in cube_dicts]

    # prepare output table. output will be table with wavelength in
    # first column, and the flux for each aperture in the following
    # columns
    num_wavs = sum((len(w) for w in wavtables))
    num_spectra = len(sky_apertures)
    output = np.zeros((num_wavs, 1 + num_spectra))

    # perform aperture photometry on every slice of every cube
    row_counter = 0
    for i in range(len(cubes)):
        for w in range(len(wavtables[i])):
            output[row_counter, 0] = wavtables[i][w]
            image = cubes[i].data[w]
            photometry_result = aperture_photometry(image, sky_apertures, wcs=wcses[i])
            output[row_counter, 1:] = photometry_result["aperture_sum"]
            row_counter += 1

    return output


##################################OBSOLETE######################################


def main():
    ra_center = 73.03
    dec_center = -66.923
    num_ra_pix = 15
    num_dec_pix = 10
    delt = 0.001
    apr = make_square_aperture_grid(
        ra_center, dec_center, delt, num_ra_pix, num_dec_pix
    )
    # print(apr)
    target = "hii1_hii8"
    cube_dicts = get_SAGE_cubes(target)
    quicklook_cubes(cube_dicts, apr)
    plt.title("A slice from each cube and apertures used")

    output_fn = "reprojected.fits"
    wavs, cube = reproject_and_merge_cubes(
        cube_dicts,
        ra_center,
        dec_center,
        delt,
        num_ra_pix,
        num_dec_pix,
        filename=output_fn,
    )

    # do the same for the uncertainties
    cube_unc_dicts = get_SAGE_cubes(target, uncertainty=True)
    output_fn_unc = output_fn.replace(".fits", "_unc.fits")
    wavs_unc, cube_unc = reproject_and_merge_cubes(
        cube_unc_dicts,
        ra_center,
        dec_center,
        delt,
        num_ra_pix,
        num_dec_pix,
        filename=output_fn_unc,
    )

    plt.figure()
    w = cube.shape[0] // 2
    wval = wavs[w]
    plt.imshow(cube[w])
    plt.title(f"Final cube at {wval} micron")
    plt.figure()
    plt.plot(wavs, cube.reshape(cube.shape[0], cube.shape[1] * cube.shape[2]))
    plt.title("SED for each pixel of final cube")
    plt.show()


main()
