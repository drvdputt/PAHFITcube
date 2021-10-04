from jwst.datamodels import CubeModel
from pathlib import Path
import wcshacks
from matplotlib import pyplot as plt
import numpy as np
from plotting import plot_cube

# Glossary:
# rpj = reprojected


def read_miri_cube_file(fn):
    return CubeModel(fn)


def get_miri_cubes(dn, pfx):
    """Get all the MIRI cubes.

    Order is from short to long wavelength, indexed on channel and
    grating: Ch1 S, Ch1 M, Ch1 L, Ch2 S, ...

    Parameters
    ----------
    dn : directory name

    pfx : prefix. Files are found as pfx_ch{n}-{short | medium | long}_s3d.fits.

    Returns
    -------
    cubes : list of CubeModel
    """
    names = [
        f"{pfx}_ch{n}-{g}_s3d.fits"
        for n in range(1, 5)
        for g in ["short", "medium", "long"]
    ]
    return [read_miri_cube_file(Path(dn) / fn) for fn in names]


def wavelengths(cube_model):
    """Get list of wavelengths.

    For some reason, cube_model.wavelength is not working. It just
    returns zeros.

    Returns value (instead of astropy quantity) to prevent problems down
    the line.

    """
    return cube_model.wcs.spectral.pixel_to_world(range(0, cube_model.shape[0])).value


def rpj_and_merge(cube_models, newwcs, ny, nx):
    rpj_data = [
        wcshacks.reproject_cube_data(c.data, c.wcs.sub((1, 2)), newwcs, ny, nx)
        for c in cube_models
    ]

    def plot12cubes(cbs):
        fig, axs = plt.subplots(3, 4)
        for ax, dat in zip(axs.flatten(), cbs):
            ax.imshow(dat[len(dat) // 2])
        return fig, axs

    fig, axs = plot12cubes([c.data for c in cube_models])
    for ax, cube_model in zip(axs.flatten(), cube_models):
        title = (
            "ch"
            + cube_model.meta.instrument.channel
            + " "
            + cube_model.meta.instrument.band
        )
        ax.set_title(title)
    plot12cubes(rpj_data)

    output_wavs = np.concatenate([wavelengths(c) for c in cube_models])
    output_cube_array = np.concatenate([data for data in rpj_data], axis=0)

    # sort the slices by wavelength
    order = np.argsort(output_wavs)
    output_wavs = output_wavs[order]
    output_cube_array = output_cube_array[order]

    wcshacks.write_merged_cube(
        "miri_merged.fits", output_cube_array, output_wavs, newwcs
    )
    return output_wavs, output_cube_array


def main():
    cubes = get_miri_cubes("../simulated-data-stsci/MRS2/stage3/", "Level3")
    nx, ny = 5, 5
    newwcs = wcshacks.make_ra_dec_wcs(0, 0, 1 / 3600, nx, ny)
    wavs, cube = rpj_and_merge(cubes, newwcs, ny, nx)

    plot_cube("miri_merged.fits", "Resampled MIRI cube")

    plt.show()


if __name__ == "__main__":
    main()
