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
from dataclasses import dataclass
from plotting import plot_cube

# local imports
import wcshacks


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

    output_projection = wcshacks.make_ra_dec_wcs(
        center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec
    )
    rpj_cube_set = CubeSet(
        *[
            Cube(
                file_handle=None,
                data=wcshacks.reproject_cube_data(
                    c.data, c.wcs, output_projection, npix_dec, npix_ra
                ),
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
        newwcs = cubes[0].wcs
        wcshacks.write_merged_cube(filename, output_cube_array, output_wavs, newwcs)


def main():
    ra_center = 73.03
    dec_center = -66.923
    npix_ra = 15
    npix_dec = 10
    delt = 0.001
    apr = wcshacks.make_square_aperture_grid(
        ra_center, dec_center, delt, npix_ra, npix_dec
    )
    # print(apr)
    target = "hii1_hii8"
    cube_dicts = get_SAGE_cubes(target)
    fig = quicklook_cubes(cube_dicts, apr)
    fig.suptitle("Reprojection grid")

    output_fn = "reprojected.fits"
    reproject_and_merge_cubes(
        cube_dicts,
        ra_center,
        dec_center,
        delt,
        npix_ra,
        npix_dec,
        filename=output_fn,
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
    indicate_overlap(plt.gca())
    plot_cube(output_fn.replace(".fits", "_nostitch.fits"), "not stitched")
    indicate_overlap(plt.gca())

    plt.show()


if __name__ == "__main__":
    main()
