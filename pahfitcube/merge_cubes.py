"""Generic cube combination function. Assumes that all the cubes given
on the command line can be loaded as spec using Spectrum1D, and that
their WCS can be retrieved from spec.meta['header'].

Should probably import and use haute couture by Amelie Canin (if
available somewhere).

"""
import numpy as np
from pathlib import Path
import argparse
from specutils import Spectrum1D
from astropy import units as u

from pahfitcube import wcshacks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fov_ra_dec", help="field of view in arcsec", type=float, nargs=2)
    ap.add_argument("--format", help="format for Spectrum1D.read()")
    ap.add_argument("-o", help="output directory")
    ap.add_argument(
        "--res", type=float, default=0.1, help="pixel scale in arcsec (default 0.1)"
    )
    ap.add_argument("cubes", nargs="+", type=str)
    args = ap.parse_args()
    output_dir = Path(args.o)
    output_dir.mkdir(exist_ok=True)

    # open all the cubes
    specs = [Spectrum1D.read(f) for f in args.cubes]

    # choose wcs
    newwcs, ny, nx = wcshacks.make_reasonable_wcs(specs, args.fov_ra_dec, args.res)

    # reproject the data
    rpj_specs = rpj_all(specs, newwcs, ny, nx)

    # merge into one big cube
    output_wavs, output_cube_array = merge_wav_and_data_arrays(
        [s.spectral_axis.to(u.micron).value for s in rpj_specs],
        [s.flux.value for s in rpj_specs],
    )

    # save individual reprojected cubes
    for i, s in enumerate(rpj_specs):
        rpj_fn = "rpj" + Path(args.cubes[i]).name
        wcshacks.write_cube(
            output_dir / rpj_fn,
            s.flux.value,
            s.spectral_axis.to(u.micron).value,
            newwcs,
            spectral_axis=-1,
        )

    # save the big merged cube
    wcshacks.write_cube(
        output_dir / "reprojected_allcube.fits",
        output_cube_array,
        output_wavs,
        newwcs,
        spectral_axis=-1,
    )


def rpj_all(specs, newwcs, nx, ny):
    """Parallel loop that reprojects cubes onto the same grid

    Parameters
    ----------
    specs: list of Spectrum1D cubes

    newwcs: wcs on which to reproject the cubes

    ny, nx: spatial dimensions of the new flux array

    Returns
    -------

    rpj_specs: list of Spectrum1D
        Reprojected Spectrum1D objects. Flux and uncertainty have been
        changed, Wavelength axis is still the same. Metadata is copied
        over, and the celestial WCS in the metadata has been adjusted.

    """
    rpj_specs = [wcshacks.reproject_s1d(s3d, newwcs, nx, ny) for s3d in specs]
    return rpj_specs


def merge_wav_and_data_arrays(wavs_list, rpj_data_list):
    output_wavs = np.concatenate(wavs_list)
    output_cube_array = np.concatenate([data for data in rpj_data_list], axis=-1)

    # sort the slices by wavelength
    order = np.argsort(output_wavs)
    output_wavs = output_wavs[order]
    output_cube_array = output_cube_array[..., order]

    print(
        f"merged shapes {[d.shape for d in rpj_data_list]} into cube of size {output_cube_array.shape}"
    )

    return output_wavs, output_cube_array


if __name__ == "__main__":
    main()
