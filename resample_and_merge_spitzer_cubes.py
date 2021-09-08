"""Take multiple spectral cubes, and merged them into a format suitable
for PAHFIT-cube.

For now, it will only work for the cubes of the SAGE-Spec program: LL1,
LL2, SL1, SL2. Given these four data cubes, a merged cube is created by
reprojecting each slice onto the same WCS, and merging the reprojected
slices. A WCS for the output is created manually.

A pixel-by-pixel spectral order stitching step is performed, based on
the average value in the region of wavelength overlap.

"""
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from pathlib import Path
from matplotlib import pyplot as plt
from photutils.aperture import SkyRectangularAperture
import reproject
from itertools import product
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
    titles = ["LL1", "LL2", "SL1", "SL2"]
    for i, cube in enumerate(cube_set.all_cubes()):
        ax = fig.add_subplot(2, 2, i + 1, projection=cube.wcs)
        ax.imshow(cube.data[-1], origin="lower")
        ax.grid()
        ax.coords[0].set_format_unit(u.degree, decimal=True)
        ax.coords[1].set_format_unit(u.degree, decimal=True)
        if apertures is not None:
            apertures.to_pixel(cube.wcs).plot(axes=ax, color="r", alpha=0.5)
        ax.set_title(titles[i])
    return fig


def indicate_overlap(ax):
    """Plot axvspan on ax to indicate each overlap area."""
    ax.axvspan(7.5, 7.6, color="k", alpha=0.1)
    ax.axvspan(14.2, 14.8, color="k", alpha=0.1)
    ax.axvspan(20.4, 21.1, color="k", alpha=0.1)


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


def stitch_long_to_short(rpj_cube_l, rpj_cube_s):
    """
    Stitch two cubes together, pixelwise, by scaling the one that is of
    shorter wavelengths.

    The cubes must have been reprojected onto the same pixel grid.

    Parameters
    ----------
    rpj_cube_l : Cube
        reprojected cube on the long wavelength side

    rpj_cube_s : Cube
        reprojected cube on the short wavelength side

    Returns
    -------
    rpj_cube_s_rescaled : ndarray
        the data of shorter wavelength cube of the two, for which every
        pixel has been rescaled to to match rpj_cube_l in the overlaping
        wavelength region.
    """
    lw = rpj_cube_l.wavelength
    sw = rpj_cube_s.wavelength

    lw_min = min(lw)
    sw_max = max(sw)

    # sum the data along the wavelength window
    ldata_sum = np.average(rpj_cube_l.data[lw < sw_max], axis=0)
    sdata_sum = np.average(rpj_cube_s.data[sw > lw_min], axis=0)

    # 2d array
    pixel_ratios = ldata_sum / sdata_sum

    # rescale the 3d array pixel-wise
    return rpj_cube_s.data * pixel_ratios[np.newaxis]


def reproject_and_merge_cubes(
    cube_set, center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec, filename=None
):
    """
    Reproject all four cubes onto pixel grid wcs, and merge the result.

    Result is sorted by wavelength.

    """
    if not ".fits" in filename:
        print("filename should end in .fits")
        raise

    output_projection = make_ra_dec_wcs(
        center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec
    )
    rpj_cube_set = CubeSet(
        *[
            Cube(
                file_handle=None,
                data=reproject_cube(c, output_projection, npix_dec, npix_ra),
                wavelength=c.wavelength,
                wcs=output_projection,
            )
            for c in cube_set.all_cubes()
        ]
    )

    # write out cube before stitching
    path = Path(filename)
    merge_and_write_cubes(
        rpj_cube_set.all_cubes(),
        path.parent / path.name.replace(".fits", "_nostitch.fits"),
    )

    # pixel wise stitching: match ll2 to ll1, sl1 to ll1, sl2 to sl1, in that order
    rpj_cube_set.ll2.data = stitch_long_to_short(rpj_cube_set.ll1, rpj_cube_set.ll2)
    rpj_cube_set.sl1.data = stitch_long_to_short(rpj_cube_set.ll2, rpj_cube_set.sl1)
    rpj_cube_set.sl2.data = stitch_long_to_short(rpj_cube_set.sl1, rpj_cube_set.sl2)
    merge_and_write_cubes(rpj_cube_set.all_cubes(), path)


def reproject_cube(cube, wcs, ny, nx):
    """
    Reproject every slice of cube onto wcs using ny, nx grid

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


def merge_and_write_cubes(cubes, filename):
    """Merge a set of spectral cubes and write fits file

    They need to be of the same size along the spatial axes. The WCS
    contained in the first cube will be used.

    Parameters
    ----------
    cubes: list of Cube

    """
    # put the cube data in one big array
    output_wavs = np.concatenate([c.wavelength for c in cubes])
    output_cube_array = np.concatenate([c.data for c in cubes], axis=0)

    # sort the slices by wavelength
    order = np.argsort(output_wavs)
    output_wavs = output_wavs[order]
    output_cube_array = output_cube_array[order]

    if filename is not None:
        path = Path(filename)
        # multi extension (wav list and cube)
        new_hdul = fits.HDUList()
        # cube as primary hdu
        header = cubes[0].wcs.to_header()
        header["BUNIT"] = "MJy/sr"
        # write out slice as test
        test_hdul = fits.HDUList()
        test_hdul.append(fits.PrimaryHDU(data=output_cube_array[0], header=header))
        test_hdul.writeto(path.parent / ("slice_test_" + path.name), overwrite=True)

        # manually set these cards, in an attempt to make the wavelength
        # slider work properly in DS9
        header["PC3_3"] = 1
        header["CRPIX3"] = 1
        header["CRVAL3"] = 1
        header["CTYPE3"] = "WAVE-TAB"
        header["CUNIT3"] = "um"
        header["PS3_0"] = "WCS-TAB"
        header["PS3_1"] = "WAVELENGTH"

        new_hdul.append(fits.PrimaryHDU(data=output_cube_array, header=header))
        # wavs as bintable hdu. Try to recreate the format of the
        # Spitzer cubes.
        weird_output_format = np.zeros(
            shape=(1,), dtype=[("WAVELENGTH", ">f4", (len(output_wavs), 1))]
        )
        for i in range(len(output_wavs)):
            weird_output_format["WAVELENGTH"][0][i][0] = output_wavs[i]
        wavhdu = fits.table_to_hdu(Table(data=weird_output_format))
        wavhdu.header["EXTNAME"] = "WCS-TAB"
        wavhdu.header["TUNIT1"] = "um"
        wavhdu.header["TDIM1"] = str((1, len(output_wavs)))

        new_hdul.append(wavhdu)
        # write
        new_hdul.writeto(path, overwrite=True)


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


def plot_cube(filename, name_in_title):
    """Plots some slices and SEDs in a cube"""
    wavs = Table.read(filename)["WAVELENGTH"][0].flatten()
    with fits.open(filename) as hdulist:
        cube = hdulist["PRIMARY"].data
        wcs = WCS(filename, naxis=2)

        fig = plt.figure()
        ax = fig.add_subplot(projection=wcs)
        w = cube.shape[0] // 2
        wval = wavs[w]
        ax.imshow(cube[w], origin="lower")
        ax.set_title(f"{name_in_title} at {wval:.2f} micron")
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")

        plt.figure()
        nw, ny, nx = cube.shape
        pixel_x_choice = (nx // 2, nx // 4, nx // 2 + nx // 4)
        pixel_y_choice = (ny // 2, ny // 4, ny // 2 + ny // 4)
        for (x, y) in product(pixel_x_choice, pixel_y_choice):
            plt.plot(wavs, cube[:, y, x])

        indicate_overlap(plt.gca())

        plt.xlabel("wavelength (micron)")
        plt.ylabel("pixel (MJy / sr)")
        plt.title(f"A few SEDs from {name_in_title}")


def main():
    ra_center = 73.03
    dec_center = -66.923
    npix_ra = 15
    npix_dec = 10
    delt = 0.001
    apr = make_square_aperture_grid(ra_center, dec_center, delt, npix_ra, npix_dec)
    # print(apr)
    target = "hii1_hii8"
    cube_dicts = get_SAGE_cubes(target)
    fig = quicklook_cubes(cube_dicts, apr)
    fig.suptitle("Reprojection grid")

    output_fn = "reprojected.fits"
    reproject_and_merge_cubes(
        cube_dicts, ra_center, dec_center, delt, npix_ra, npix_dec, filename=output_fn,
    )

    # do the same for the uncertainties
    cube_unc_dicts = get_SAGE_cubes(target, uncertainty=True)
    output_fn_unc = output_fn.replace(".fits", "_unc.fits")
    reproject_and_merge_cubes(
        cube_unc_dicts,
        ra_center,
        dec_center,
        delt,
        npix_ra,
        npix_dec,
        filename=output_fn_unc,
    )

    # plot a preview of the merged cube
    plot_cube(output_fn, "stitched cube")
    plot_cube(output_fn.replace(".fits", "_nostitch.fits"), "not stitched")

    plt.show()


main()
