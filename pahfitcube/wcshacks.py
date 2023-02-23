from astropy.wcs import WCS
import reproject
import numpy as np
import astropy.units as u
from photutils import SkyRectangularAperture
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
import pickle
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
from multiprocess import Pool


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


def reproject_cube_data(cube_data, cube_wcs, wcs, n0, n1):
    """Reproject every slice of cube onto wcs using n0, n1 grid

    This function assumes the spectrum1D convention of putting the
    wavelength index last.

    The two spatial dimensions are axis 0 and axis 1

    Returns
    -------
    output_array: np.ndarray indexed on axis 0, axis 1, wavelength

    """
    num_wavs = cube_data.shape[-1]
    output_array = np.zeros((n0, n1, num_wavs))
    for w in range(num_wavs):
        output_array[..., w], footprint = reproject.reproject_adaptive(
            input_data=(cube_data[..., w], cube_wcs),
            output_projection=wcs,
            shape_out=(n0, n1),
        )
    return output_array


def celestial_wcs_from_s1d(s):
    """s: Spectrum1D with 3D wcs in meta['header']"""
    return WCS(s.meta["header"]).sub((1, 2))


def reproject_s1d(s3d, wcs, nx, ny):
    """Reproject every slice of Spectrum1D cube onto wcs using ny, nx grid

    Reprojects both flux and uncertainty, and creates new Spectrum1D object. Metadata is copied over.

    Returns
    -------
    new_s3d: new Spectrum1D object
    """
    old_wcs = celestial_wcs_from_s1d(s3d)
    rpj_flux = reproject_cube_data(s3d.flux.value, old_wcs, wcs, nx, ny)

    if s3d.uncertainty is not None:
        rpj_unc = StdDevUncertainty(
            reproject_cube_data(s3d.uncertainty.array, old_wcs, wcs, nx, ny)
        )
    else:
        rpj_unc = None

    new_s3d = Spectrum1D(
        rpj_flux * s3d.flux.unit, s3d.spectral_axis, uncertainty=rpj_unc, meta=s3d.meta
    )
    add_celestial_wcs_to_s1d(new_s3d, wcs)
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
    for key in spatial_header:
        s3d.meta["header"][key] = spatial_header[key]


def write_wavetab_cube(fn, cube_data, wave, spatial_wcs, wav_axis_index=-1):
    """Write out cube data with unevenly spaced wavelengths

    Based on some code I found in the JWST package, cube_build/ifu_cube.py

    wav_axis_index default is -1, to conform with spectrum1D.
    """
    if wav_axis_index != 0:
        data = np.moveaxis(cube_data, wav_axis_index, 0)
    else:
        data = cube_data

    num = len(wave)
    header = spatial_wcs.to_header()
    header["BUNIT"] = "MJy/sr"
    # translate these statements from the jwst package
    # ifucube_model.meta.wcsinfo.ctype3 = 'WAVE-TAB'
    #   ifucube_model.meta.wcsinfo.ps3_0 = 'WCS-TABLE'
    #   ifucube_model.meta.wcsinfo.ps3_1 = 'wavelength'
    #   ifucube_model.meta.wcsinfo.crval3 = 1.0
    #   ifucube_model.meta.wcsinfo.crpix3 = 1.0
    #   ifucube_model.meta.wcsinfo.cdelt3 = None
    #   ifucube_model.meta.ifu.roi_wave = np.mean(self.roiw_table)
    #   ifucube_model.wavedim = '(1,{:d})'.format(num)
    wavedim = "(1,{:d})".format(num)
    header["CTYPE3"] = "WAVE-TAB"
    header["PS3_0"] = "WCS-TABLE"
    header["PS3_1"] = "wavelength"
    header["CRVAL3"] = 1.0
    header["CRPIX3"] = 1.0
    # header["CDELT3"] = None
    #   ifucube_model.meta.wcsinfo.cunit3 = 'um'
    #   ifucube_model.meta.wcsinfo.wcsaxes = 3
    #   ifucube_model.meta.wcsinfo.pc1_1 = -1
    #   ifucube_model.meta.wcsinfo.pc1_2 = 0
    #   ifucube_model.meta.wcsinfo.pc1_3 = 0
    header["CUNIT3"] = "um"
    header["WCSAXES"] = 3
    header["PC1_1"] = -1
    header["PC1_2"] = 0
    header["PC1_3"] = 0
    #   ifucube_model.meta.wcsinfo.pc2_1 = 0
    #   ifucube_model.meta.wcsinfo.pc2_2 = 1
    #   ifucube_model.meta.wcsinfo.pc2_3 = 0
    header["PC2_1"] = 0
    header["PC2_2"] = 1
    header["PC2_3"] = 0
    #   ifucube_model.meta.wcsinfo.pc3_1 = 0
    #   ifucube_model.meta.wcsinfo.pc3_2 = 0
    #   ifucube_model.meta.wcsinfo.pc3_3 = 1
    header["PC3_1"] = 0
    header["PC3_2"] = 0
    header["PC3_3"] = 1

    # this header is now ready. It is attached to the main data
    # header["EXTNAME"] = 'SCI'
    new_hdul = fits.HDUList()
    new_hdul.append(fits.ImageHDU(data=data, header=header, name="SCI"))

    # jwst package creates wavelength table like this
    alldata = np.array([(wave[None].T,)], dtype=[("wavelength", "<f4", (num, 1))])
    # but does not show how it is actually put into fits file
    # wavhdu = fits.table_to_hdu(Table(data=weird_output_format))
    wavhdu = fits.BinTableHDU(alldata, name="WCS-TABLE")
    # wavhdu.header["TUNIT1"] = "um"
    wavhdu.header["TDIM1"] = wavedim
    wavhdu.header["TDIM2"] = wavedim
    wavhdu.header["TTYPE1"] = "wavelength"
    new_hdul.append(wavhdu)

    new_hdul.writeto(fn, overwrite=True)


def write_cube(fn, data, wavs, spatial_wcs, spectral_axis=None):
    """Write out cube in fits format (for DS9) and as a pickle (for Spectrum1D)

    The pickle can be given to run_pahfit_cube, which will use its
    contents to make a Spectrum1D and obtain the WCS.

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

    write_wavetab_cube(fn, data, wavs, spatial_wcs, spectral_axis)

    # Ideally, I want to make a Spectrum1D, and save it as a fits file
    # that readable by both DS9 and Spectrum1D.read(). But that doesn't
    # seem supported yet, so I resort to pickling the input data for
    # Spectrum1D below.
    if spectral_axis is not None:
        spec1d_flux_array = np.moveaxis(data, spectral_axis, -1)
        # spec1d_flux_array = np.moveaxis(data, spectral_axis, 0)
    else:
        spec1d_flux_array = data

    # MJy/sr
    spec1d_flux_array_with_units = spec1d_flux_array * 1e6 * u.Jy / u.sr

    pickle_path = path.with_suffix(".pickle")
    with open(pickle_path, "wb") as f:
        # Spectrum1D itself isn't pickleable. To make one, we need flux,
        # wavelengths and spatial wcs. And perhaps uncertainties and
        # data quality too, if we know how to calculate those for the
        # merged cubes.
        obj = {
            "flux": spec1d_flux_array_with_units,
            "spectral_axis": wavs * u.micron,
            "spatial_wcs": spatial_wcs,
        }
        pickle.dump(obj, f)


def write_merged_cube_old(fn, data, wavs, spatial_wcs):
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
    new_hdul.writeto(fn, overwrite=True)
