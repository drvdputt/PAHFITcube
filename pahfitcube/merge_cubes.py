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
from astropy.wcs import WCS
from pahfitcube import wcshacks, iohacks
from scipy.interpolate import interp1d
from astropy.nddata import StdDevUncertainty


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

    # save individual reprojected cubes
    for i, s in enumerate(rpj_specs):
        rpj_fn = "rpj" + Path(args.cubes[i]).name
        iohacks.write_s3d(output_dir / rpj_fn, s, newwcs)

    # merge the reprojected segments
    output_spec = merge_nd(rpj_specs)
    iohacks.write_s3d(output_dir / "reprojected_allcube.fits", output_spec, newwcs)


def rpj_all(specs, new_wcs, nx, ny):
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
    rpj_specs = [wcshacks.reproject_s1d(s3d, new_wcs, nx, ny) for s3d in specs]
    return rpj_specs


def find_overlap_ranges(ss):
    """Find the wavelength overlap regions of a list of spectra.

    Parameters
    ----------

    ss: list of Spectrum1D
        Assumes that the spectra are already sorted, and that each
        spectral segment overlaps only with the previous and
        the next in the list.

    Returns
    -------

    list of 2-tuples representing the ranges where spectral overlap occurs.
        typically [(min1, max0), (min2, max1), ...]

    """
    wav_overlap_ranges = []
    for i in range(1, len(ss)):
        pr = ss[i - 1]
        cr = ss[i]
        v1 = cr.spectral_axis[0]
        v2 = pr.spectral_axis[-1]
        if v2 > v1:
            wav_overlap_ranges.append((v1, v2))

    wav_channel_span = []
    for i in range(0, len(ss), 3):
        short = ss[i]
        long = ss[i + 2]
        wav_channel_span.append((short.spectral_axis[0], long.spectral_axis[-1]))

    return wav_overlap_ranges


def merge_nd(ss):
    """Merge a list of sorted (by wavelength) Spectrum1D segments

    A newer merging method, copied over from my personal astro package.
    Assumes that no spectral matching correction is needed. At the
    moment there are multiple merging methods in this module, but this
    one should supersede the other eventually.

    Parameters
    ----------
    ss: sorted segments (list of Spectrum1D), of same spatial shape

    Returns
    -------
    Spectrum1D

    """
    # find the wavelength regions where the segments overlap
    overlap_ranges = find_overlap_ranges(ss)

    # merge everything
    new_spectral_axis = np.sort(np.concatenate([s.spectral_axis for s in ss]))

    # replace the points in the overlap ranges by the finest of the two
    # relevant segments
    for ileft, (wmin, wmax) in enumerate(overlap_ranges):
        # keep the points outside the overlap region
        wavs_outside = new_spectral_axis[
            (new_spectral_axis < wmin) | (new_spectral_axis > wmax)
        ]
        # replace the points inside the overlap region with points from the left segment
        left_wavs = ss[ileft].spectral_axis
        wavs_inside = left_wavs[(left_wavs > wmin) & (left_wavs < wmax)]
        new_spectral_axis = np.sort(np.concatenate([wavs_outside, wavs_inside]))

    flux_new = [
        interp1d(
            s.spectral_axis.value,
            s.flux.value,
            axis=-1,  # this is default, but I'm being explicit here for clarity
            bounds_error=False,
            fill_value=np.nan,
        )(new_spectral_axis.value)
        for s in ss
    ]
    unc_new = [
        interp1d(
            s.spectral_axis.value,
            s.uncertainty.array,
            axis=-1,
            bounds_error=False,
            fill_value=np.nan,
        )(new_spectral_axis.value)
        for s in ss
    ]

    # first, naively combine into one big spectrum without caring about the jumps
    # At this point, flux_new has indices [segment, space, space, wavelength]
    flux_merged = np.nanmedian(flux_new, axis=0)

    for ileft, (wmin, wmax) in enumerate(overlap_ranges):
        f_left = flux_new[ileft]
        f_right = flux_new[ileft + 1]
        overlap = np.logical_and(new_spectral_axis > wmin, new_spectral_axis < wmax)

        # f(w) = 0 at wmin, 1 at wmax
        sliding_f = (new_spectral_axis[overlap] - wmin) / (wmax - wmin)
        flux_merged[..., overlap] = (1 - sliding_f) * f_left[
            ..., overlap
        ] + sliding_f * f_right[..., overlap]

    new_uncertainty_sum = np.sqrt(np.nansum([a**2 for a in unc_new], axis=0))
    new_uncertainty_count = np.count_nonzero(np.isfinite([a for a in unc_new]), axis=0)
    new_uncertainty = StdDevUncertainty(new_uncertainty_sum / new_uncertainty_count)

    return Spectrum1D(
        flux_merged * ss[0].flux.unit, new_spectral_axis, uncertainty=new_uncertainty
    )


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


def make_usable_cube(cubes, wmin, wmax):
    """Exemplary workflow of how to deal with a set of cubes, to get data withing a certain wavelength range.
    1. Determine which cubes we need based on wmin/wmax
    2. Reproject the cubes onto spatial grid of the longest wavelength cube
    3. Merge the cubes of the same shape. Need something like merge_1d that also works in general.
    4. Cut out the right wavelength range
    5. Save to disk, so we don't have to run this function in the pahfitcube run script.

    Parameters
    ----------
    cubes: Spectrum1D
        Input cubes of different shapes and wavelength ranges. It is
        assumes that the cubes are already sorted by wavelength.

    wmin, wmax: float
        Wavelength range to extract, in micron. Only cubes with data in
        this range will be chosen from the provided list of cubes.

    Returns
    -------
    Spectrum1D: a single merged cube cut to the desired wavelength range
    """
    use_cubes = [
        c
        for c in cubes
        if np.any(
            (c.spectral_axis > wmin * u.micron) & (c.spectral_axis < wmax * u.micron)
        )
    ]

    # reproject everything to the WCS of the longest relevant cube
    clong = use_cubes[0]
    cwcs = WCS(clong.meta["header"]).celestial
    nx = clong.shape[0]
    ny = clong.shape[1]
    rpj_cubes = rpj_all(cubes, cwcs, nx, ny)

    merged_cube = merge_nd(rpj_cubes)
    return merged_cube[wmin * u.micron : wmax * u.micron], cwcs


if __name__ == "__main__":
    main()
