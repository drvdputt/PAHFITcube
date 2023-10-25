"""Tricks for saving and loading cube files, either pickle or fits."""
from specutils import Spectrum1D
from astropy.wcs import WCS
from astropy.nddata import StdDevUncertainty
import pickle
import numpy as np
from astropy.io import fits
from pathlib import Path
from astropy import units as u
from astropy.table import Table


def _write_wavetab_cube(fn, flux, uncertainty, wave, spatial_wcs, wav_axis_index=-1):
    """Write out cube data with unevenly spaced wavelengths

    Based on some code I found in the JWST package, cube_build/ifu_cube.py

    wav_axis_index default is -1, to conform with spectrum1D.
    """
    if wav_axis_index != 0:
        f = np.moveaxis(flux, wav_axis_index, 0)
        unc = np.moveaxis(uncertainty, wav_axis_index, 0)
    else:
        f = flux
        unc = uncertainty

    f = np.swapaxes(f, 1, 2)
    unc = np.swapaxes(f, 1, 2)

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
    # header["PC1_1"] = -1
    # header["PC1_2"] = 0
    header["PC1_3"] = 0
    #   ifucube_model.meta.wcsinfo.pc2_1 = 0
    #   ifucube_model.meta.wcsinfo.pc2_2 = 1
    #   ifucube_model.meta.wcsinfo.pc2_3 = 0
    # header["PC2_1"] = 0
    # header["PC2_2"] = 1
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
    new_hdul.append(fits.ImageHDU(data=f, header=header, name="SCI"))
    new_hdul.append(fits.ImageHDU(data=unc, header=header, name="ERR"))

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


def write_cube(fn, flux, uncertainty, wavs, spatial_wcs, spectral_axis=None):
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

    uncertainty: array (no unit!)
        uncertainty array to store. Same conventions as flux data

    wavs: array (no unit!)
        wavelengths in micron

    spatial_wcs: wcs for the celestial axes. Typically
        WCS(header).celestial for most cube fits files.

    spectra_axis: int
        if spectral axis is not last axis, indicate the right axis here.

    """
    if isinstance(fn, Path):
        path = fn
    else:
        path = Path(fn)

    _write_wavetab_cube(
        fn,
        flux,
        uncertainty,
        wavs,
        spatial_wcs,
        -1 if spectral_axis is None else spectral_axis,
    )

    # Ideally, I want to make a Spectrum1D, and save it as a fits file
    # that readable by both DS9 and Spectrum1D.read(). But that doesn't
    # seem supported yet, so I resort to pickling the input data for
    # Spectrum1D below.
    if spectral_axis is not None:
        spec1d_flux_array = np.moveaxis(flux, spectral_axis, -1)
        spec1d_unc_array = np.moveaxis(uncertainty, spectral_axis, -1)
        # spec1d_flux_array = np.moveaxis(data, spectral_axis, 0)
    else:
        spec1d_flux_array = flux
        spec1d_unc_array = uncertainty

    # MJy/sr
    spec1d_flux_array_with_units = spec1d_flux_array * u.MJy / u.sr
    spec1d_unc_array_with_units = spec1d_unc_array * u.MJy / u.sr

    pickle_path = path.with_suffix(".pickle")
    with open(pickle_path, "wb") as f:
        # Spectrum1D itself isn't pickleable. To make one, we need flux,
        # wavelengths and spatial wcs. And perhaps uncertainties and
        # data quality too, if we know how to calculate those for the
        # merged cubes.
        obj = {
            "flux": spec1d_flux_array_with_units,
            "uncertainty": spec1d_unc_array_with_units,
            "spectral_axis": wavs * u.micron,
            "spatial_wcs": spatial_wcs,
        }
        pickle.dump(obj, f)


def write_s3d(fn, s: Spectrum1D, spatial_wcs):
    """Utility to write cube from Spectrum1D

    For now, due to hardcoding, this will only produce correct results
    if the flux is in MJy / sr and the spectral axis is in micron.

    will write both a fits file for DS9, and pickle that PAHFITcube can load

    While the spatial WCS could be retrieved from the Spectrum1D meta
    data, were asking for it explicitly here, since I found it not to
    always be reliable.

    """
    write_cube(
        fn, s.flux.value, s.uncertainty.array, s.spectral_axis.value, spatial_wcs
    )


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


def read_cube(cubefile):
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
            uncertainty = StdDevUncertainty(data_dict["uncertainty"])
            spec = Spectrum1D(flux, spectral_axis, uncertainty=uncertainty)
            cube_2dwcs = data_dict["spatial_wcs"]
        except Exception as e:
            print(
                "Something wrong with pickle. Needs to come from iohacks.write_merged_cube."
            )
            raise e
    else:
        raise ValueError("Cube file not compatible.")

    ny, nx, _ = spec.shape
    return spec, cube_2dwcs
