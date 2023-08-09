import numpy as np
from astropy.wcs.utils import proj_plane_pixel_area
from specutils import Spectrum1D
from astropy.wcs import WCS
from astropy import units as u
from astropy.nddata import StdDevUncertainty


def wcs_from_spec1d(spec1d):
    """Create a WCS from header in metadata of a Spectrum1D cube

    Returns
    -------

    Celestial wcs

    """
    return WCS(spec1d.meta["header"])


def make_cube_array_mask(wcs_2d, shape_2d, sky_aperture):
    """Make an array representing an aperture discretized on the spatial wcs/grid.

        Parameters
        ----------
        wcs_2d: WCS
            the celestial wcs of a cube
    _
        shape_2d: (int, int)
            the celestial part of the shape of a cube

        Returns
        -------
        array_mask: np.array of shape shape_2d.

    """
    try:
        pixel = sky_aperture.to_pixel(wcs_2d)
    except ValueError as e:
        raise ValueError(
            "Aperture to_pixel failed. Probably because aperture is not inside FOV of WCS"
        ) from e
    pixel_mask = pixel.to_mask(method="exact")

    # watch out here! To_image works in y,x coordinates! Need to provide
    # shape as (y,x), and then convert the image mask to a mask that
    # works for our array, by transposing.
    image_mask = pixel_mask.to_image(shape_2d[::-1])

    if image_mask is None:
        print("Something wrong with make overlap!")
        print("sky aperture = ", sky_aperture)
        print("pixel aperture = ", pixel)
        print("WCS = ", wcs_2d)

    array_mask = image_mask.T
    return array_mask


def cube_sky_aperture_plot(ax1, ax2, cube_spec1d, sky_aperture, wcs_2d=None):
    """Plot aperture on top of cube slice

    For the background image, the middle slice of the cube is used
    (along the wavelength direction)"""
    if wcs_2d is None:
        the_wcs_2d = wcs_from_spec1d(cube_spec1d).celestial
    else:
        the_wcs_2d = wcs_2d

    array_mask = make_cube_array_mask(the_wcs_2d, cube_spec1d.shape[:2], sky_aperture)
    example_image = cube_spec1d.flux.value[:, :, cube_spec1d.shape[2] // 2] * (
        1 + array_mask
    )
    sky_aperture.to_pixel(the_wcs_2d).plot(ax1, color="red")
    ax1.imshow(np.log10(example_image).T, origin="lower")
    ax2.imshow(array_mask.T, origin="lower")


def cube_sky_aperture_extraction(
    cube_spec1d, sky_aperture, average_per_sr=True, wcs_2d=None
):
    """Extract a 1D spectrum from a cube using an aperture over all the slices.

    Parameters
    ----------
    cube_spec1d: Spectrum1D
        Object loaded from cube file using Spectrum1D.read

    sky_aperture: SkyAperture
        The aperture to apply. Is converted to pixel mask suitable for the cube.

    Returns
    -------
    spectrum: Spectrum1D
        The collapsed spectrum.

    """
    # make 2D wcs of cube
    if wcs_2d is None:
        the_wcs_2d = wcs_from_spec1d(cube_spec1d).celestial
    else:
        the_wcs_2d = wcs_2d

    # use wcs and shape to make mask
    array_mask = make_cube_array_mask(the_wcs_2d, cube_spec1d.shape[:2], sky_aperture)

    # broadcast the mask over the slices and multiply (and turn the data
    # in to a masked cube to avoid problems with nan)
    masked_cube = array_mask[:, :, None] * np.ma.masked_invalid(cube_spec1d.data)

    # Add units again, and convert to MJy
    area = (proj_plane_pixel_area(the_wcs_2d) * u.deg**2).to(u.sr)
    masked_cube = masked_cube * cube_spec1d.flux.unit * area

    # collapse
    spectrum = np.sum(masked_cube, axis=(0, 1))

    # do the same for uncertainty
    if cube_spec1d.uncertainty is not None:
        masked_unc = array_mask[:, :, None] * np.ma.masked_invalid(
            cube_spec1d.uncertainty.array
        )
        masked_unc = masked_unc * cube_spec1d.flux.unit * area
        # collapse the uncertainty. Variances = (mask value * sigma) **2
        sigmas = np.sqrt(np.sum(np.square(masked_unc), axis=(0, 1)))
    else:
        sigmas = None

    if average_per_sr:
        # total area of the unmasked pixels. Note that array_mask is
        # generally 0 or 1, but at the edges, there can be fractional
        # values.
        masked_area = area * np.sum(array_mask)
        spectrum /= masked_area
        if sigmas is not None:
            sigmas /= masked_area

    # make a spectrum1d object for convenience
    s1d = Spectrum1D(
        spectral_axis=cube_spec1d.spectral_axis,
        flux=spectrum,
        uncertainty=StdDevUncertainty(sigmas) if sigmas is not None else None,
    )
    return s1d
