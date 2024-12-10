from astropy.wcs import WCS
import reproject
import numpy as np
import astropy.units as u
from photutils.aperture import SkyRectangularAperture
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
from copy import deepcopy


def cube_wcs_extent(wcs, shape):
    """Get center and four corners"""
    # get only the celestial wcs
    wcs2d = wcs.sub((1, 2))

    # the four corners in pixel space
    corners = [(0, 0), (0, shape[1]), (shape[0], 0), (shape[0], shape[1])]
    sky_corners = wcs2d.pixel_to_world(corners)
    # sky_center = wcs2d.pixel_to_world([(shape[0] // 2, shape[1] // 2)])
    return sky_corners


def make_reasonable_wcs(spec1d_cubes, fov, res=None):
    """
    fov = (ra span, dec span) in arcsec
    """
    # extents = [
    #     cube_wcs_extent(WCS(s.meta["header"]), s.shape[:2]) for s in spec1d_cubes
    # ]

    # choose center of first cube
    header = spec1d_cubes[0].meta["header"]
    center_ra = header["CRVAL1"]
    center_dec = header["CRVAL2"]

    # size choice is manual for now
    fov_ra = fov[0] / 3600
    fov_dec = fov[1] / 3600

    # determine number of pixels along each axis
    if res is None:
        # default to .1 arcsec
        pix_degree_delta = 0.1 / 3600
    else:
        pix_degree_delta = res / 3600

    nx = int(fov_ra / pix_degree_delta)
    ny = int(fov_dec / pix_degree_delta)

    return make_ra_dec_wcs(center_ra, center_dec, pix_degree_delta, nx, ny), ny, nx


def make_ra_dec_wcs(center_ra, center_dec, pix_degree_delta, npix_ra, npix_dec):
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
    w.wcs.cdelt = [-pix_degree_delta, pix_degree_delta]
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


def reproject_cube_data(cube_data, cube_wcs, new_wcs, n0, n1):
    """Reproject every slice of cube onto wcs using n0, n1 grid

    This function assumes the spectrum1D convention of putting the
    wavelength index last.

    The two spatial dimensions are axis 0 and axis 1

    Returns
    -------
    output_array: np.ndarray indexed on axis 0, axis 1, wavelength

    """
    # the reprojection is vectorized over the first axes. So if we move
    # wavelength to the first index, we should be able to do all slices
    # at once. At the end, move the wavelength axis back to the last
    # index.
    input_data = np.moveaxis(cube_data, 2, 0)
    output_array, _ = reproject.reproject_adaptive(
        input_data=(input_data, cube_wcs),
        output_projection=new_wcs,
        shape_out=(input_data.shape[0], n0, n1),
    )
    return np.moveaxis(output_array, 0, 2)


def celestial_wcs_from_s3d(s):
    """s: Spectrum1D with 3D wcs in meta['header']"""
    # return WCS(s.meta["header"]).sub((1, 2))
    return WCS(s.meta["header"]).celestial


def reproject_s1d(s3d, new_wcs, nx, ny):
    """Reproject every slice of Spectrum1D cube onto wcs using ny, nx grid

    Reprojects both flux and uncertainty, and creates new Spectrum1D
    object. Metadata is copied over.

    Returns
    -------
    new_s3d: new Spectrum1D object

    """
    old_wcs = celestial_wcs_from_s3d(s3d)
    print("reprojecting from", old_wcs, "to", new_wcs)
    rpj_flux = reproject_cube_data(s3d.flux.value, old_wcs, new_wcs, nx, ny)

    if s3d.uncertainty is not None:
        rpj_unc = StdDevUncertainty(
            np.sqrt(
                reproject_cube_data(
                    np.square(s3d.uncertainty.array), old_wcs, new_wcs, nx, ny
                )
            )
        )
    else:
        rpj_unc = None

    # deepcopies just to be sure. I've been burned by this before.
    new_s3d = Spectrum1D(
        rpj_flux * s3d.flux.unit,
        deepcopy(s3d.spectral_axis),
        uncertainty=rpj_unc,
        meta=deepcopy(s3d.meta),
    )
    add_celestial_wcs_to_s1d(new_s3d, new_wcs)
    return new_s3d


# def reproject_s1d_multi(s3ds, wcs, nx, ny):
#     """Parallel version, since this is quite slow"""
#     args = ((s.flux.value, celestial_wcs_from_s1d(s), ny, nx) for s in s3ds)
#     with Pool(8) as p:
#         rpj_data = p.starmap(
#             reproject_cube_data,
#             [
#                 (s.data, WCS(s.meta["header"]).sub((1, 2)), newwcs, ny, nx)
#                 for s in specs
#             ],
#         )

#     rpj_fluxes = Pool


def add_celestial_wcs_to_s1d(s3d, wcs2d):
    """Add 2D wcs header to Spectrum1D.meta['header']"""
    spatial_header = wcs2d.to_header()
    # deepcopy because if meta comes from another s3d object, both of
    # their WCS will be changed
    new_meta = deepcopy(s3d.meta)
    for key in spatial_header:
        new_meta["header"][key] = spatial_header[key]
    s3d.meta = new_meta
