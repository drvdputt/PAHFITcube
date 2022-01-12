from astropy.wcs import WCS
import reproject
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
import numpy as np
import astropy.units as u
from photutils import SkyRectangularAperture


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


def reproject_cube_data(cube_data, cube_wcs, wcs, ny, nx):
    """
    Reproject every slice of cube onto wcs using ny, nx grid

    Returns
    -------
    output_array: np.ndarray indexed on wavelength, y, x
    """
    num_wavs = cube_data.shape[0]
    output_array = np.zeros((num_wavs, ny, nx))
    for w in range(num_wavs):
        output_array[w], footprint = reproject.reproject_adaptive(
            input_data=(cube_data[w], cube_wcs),
            output_projection=wcs,
            shape_out=(ny, nx),
        )
    return output_array


def write_merged_cube(fn, data, wavs, spatial_wcs, spectral_axis=None):
    """Use Spectrum1D to do write out cube

    Parameters
    ----------

    fn: str
        file name

    data: array (no unit!)
        flux (in MJy/sr). Last axis must be spectral, or spectral_axis
        should be set, so that it can be moved there.

    wavs: array (no unit!)
        wavelengths in micron

    """

    if isinstance(fn, Path):
        path = fn
    else:
        path = Path(fn)

    new_hdul = fits.HDUList()
    header = spatial_wcs.to_header()
    header["BUNIT"] = "MJy/sr"
    # manually set these cards, but still can't seem to make the
    # wavelength slider work properly in DS9
    header["PC3_3"] = 1
    header["CRPIX3"] = 1
    header["CRVAL3"] = 1
    header["CTYPE3"] = "WAVE-TAB"
    header["CUNIT3"] = "um"
    header["PS3_0"] = "WCS-TAB"
    header["PS3_1"] = "WAVELENGTH"
    # wavs as bintable hdu. Try to recreate the format of the
    # Spitzer cubes.
    new_hdul.append(fits.PrimaryHDU(data=data, header=header))
    weird_output_format = np.zeros(
        shape=(1,), dtype=[("WAVELENGTH", ">f4", (len(wavs), 1))]
    )
    for i in range(len(wavs)):
        weird_output_format["WAVELENGTH"][0][i][0] = wavs[i]
    wavhdu = fits.table_to_hdu(Table(data=weird_output_format))
    wavhdu.header["EXTNAME"] = "WCS-TAB"
    wavhdu.header["TUNIT1"] = "um"
    wavhdu.header["TDIM1"] = str((1, len(wavs)))
    new_hdul.append(wavhdu)
    new_hdul.writeto(path, overwrite=True)
