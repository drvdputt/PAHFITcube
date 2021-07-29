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


def read_spitzer_cube(fn):
    cube = fits.open(fn)[0]
    wavs = Table.read(fn)["WAVELENGTH"][0]
    wcs = WCS(fn, naxis=2)
    return dict(filename=fn, cube=cube, wavelength=wavs, wcs=wcs)


def get_SAGE_cubes(target):
    dname = f"data/sage-spec_{target}_4dec08"
    d = Path(dname)
    cube_names = ["ll_LL1", "ll_LL2", "sl_SL1", "sl_SL2"]
    cubes = []
    for cube_name in cube_names:
        fname = d / f"{target}_{cube_name}_cube.fits"
        cubes.append(read_spitzer_cube(str(fname)))

    return cubes


def quicklook_cubes(cube_dicts, apertures=None):
    fig = plt.gcf()
    for i, cube in enumerate(cube_dicts):
        ax = fig.add_subplot(2, 2, i + 1, projection=cube["wcs"])
        ax.imshow(cube["cube"].data[-1])
        ax.grid()
        ax.coords[0].set_format_unit(u.degree, decimal=True)
        ax.coords[1].set_format_unit(u.degree, decimal=True)
        if apertures is not None:
            apertures.to_pixel(cube["wcs"]).plot(axes=ax, color="r")


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
    center_x = center_ra / 2 + 0.5
    center_y = center_dec / 2 + 0.5
    w.wcs.crpix = [center_x, center_y]
    w.wcs.crval = [center_ra, center_dec]
    w.wcs.cdelt = [pix_angle_delta, pix_angle_delta]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def reproject_and_merge(
    cube_dicts, center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec
):
    """Reproject all four cubes onto pixel grid wcs, and merge the result."""
    output_projection = make_ra_dec_wcs(
        center_ra, center_dec, pix_angle_delta, npix_ra, npix_dec
    )
    wavtables = [d["wavelength"] for d in cube_dicts]
    output_wavs = np.concatenate(wavtables)
    nw = len(output_wavs)

    ny, nx = npix_dec, npix_ra
    output_cube_array = np.zeros((nw, ny, nx))

    start = 0
    for d in cube_dicts:
        stop = start + len(d["wavelength"])
        output_cube_array[start:stop] = reproject_cube(d, output_projection, ny, nx)
        start = stop

    return output_wavs, output_cube_array


def reproject_cube(cube_dict, wcs, ny, nx):
    """Reproject every slice of cube onto wcs, ny, nx.

    Returns
    -------
    output_array: np.ndarray indexed on wavelength, y, x
    """
    num_wavs = len(cube_dict["wavelength"])
    input_array = cube_dict["cube"].data
    input_wcs = cube_dict["wcs"]

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
    """Create sky apertures representing pixels of RA-DEC aligned map.

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


def main():
    ra_center = 73.03
    dec_center = -66.923
    num_ra_pix = 1
    num_dec_pix = 1
    delt = 0.01
    apr = make_square_aperture_grid(
        ra_center, dec_center, delt, num_ra_pix, num_dec_pix
    )
    print(apr)
    c = get_SAGE_cubes("hii1_hii8")
    quicklook_cubes(c)
    # quicklook_cubes(c, apr)
    # result = extract_spectra(c, apr)
    # print(result)
    # plt.figure()
    # for i in range(1, result.shape[1]):
    #     plt.plot(result[:, 0], result[:, i])
    # plt.show()

    # Different approach: try using reproject package
    wavs, cube = reproject_and_merge(
        c, ra_center, dec_center, delt, num_ra_pix, num_dec_pix
    )
    print(wavs)
    print(cube)
    plt.figure()
    plt.imshow(cube[cube.shape[0] // 2])
    plt.figure()
    plt.plot(wavs, cube.reshape(cube.shape[0], cube.shape[1] * cube.shape[2]))
    plt.show()

main()
