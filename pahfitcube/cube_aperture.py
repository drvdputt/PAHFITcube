import numpy as np
from astropy.wcs.utils import proj_plane_pixel_area
from specutils import Spectrum1D
from astropy.wcs import WCS
from astropy import units as u


def wcs_from_spec1d(spec1d):
    """Create a WCS from the header stored in the metadata of a Spectrum1D
    object representing a cube"""
    return WCS(spec1d.meta["header"])


def make_cube_image_mask(wcs_2d, shape_2d, sky_aperture):
    # use wcs to make pixel mask
    pixel_mask = sky_aperture.to_pixel(wcs_2d).to_mask(method="exact")

    # use shape to make image mask
    image_mask = pixel_mask.to_image(shape_2d)

    if image_mask is None:
        print("Something wrong with make overlap!")

    return image_mask


def cube_sky_aperture_extraction(cube_spec1d, sky_aperture):
    # make 2D wcs of cube
    wcs_2d = wcs_from_spec1d(cube_spec1d).celestial

    # use wcs and shape to make mask
    image_mask = make_cube_image_mask(wcs_2d, cube_spec1d.shape[:2], sky_aperture)

    # broadcast the mask over the slices and multiply
    masked_cube = image_mask[:, :, None] * cube_spec1d.data

    # collapse the masked cube
    spectrum = np.average(masked_cube, axis=(0, 1))

    # convert from average MJy sr-1 to total flux in MJy
    # average per sr * total sr = average * area * number
    pix_area = (proj_plane_pixel_area(wcs_2d) * u.deg**2).to(u.sr)
    spectrum = spectrum * pix_area * cube_spec1d.shape[0] * cube_spec1d.shape[1]

    # make a spectrum1d object for convenience
    s1d = Spectrum1D(
        spectral_axis=cube_spec1d.spectral_axis, flux=spectrum * cube_spec1d.flux.unit
    )

    return s1d
