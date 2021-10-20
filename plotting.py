from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
from itertools import product
from astropy.wcs import WCS


def plot_cube(filename, name_in_title):
    """Plots some slices and SEDs in a cube"""
    wavs = Table.read(filename)["WAVELENGTH"][0].flatten()
    with fits.open(filename) as hdulist:
        cube = hdulist["PRIMARY"].data
        wcs = WCS(filename, naxis=2)

        fig = plt.figure()
        ax = fig.add_subplot(projection=wcs)
        w = cube.shape[0] // 2
        wval = wavs[w]
        ax.imshow(cube[w], origin="lower")
        ax.set_title(f"{name_in_title} at {wval:.2f} micron")
        ax.set_xlabel("RA")
        ax.set_ylabel("DEC")

        plt.figure()
        nw, ny, nx = cube.shape
        pixel_x_choice = (nx // 2, nx // 4, nx // 2 + nx // 4)
        pixel_y_choice = (ny // 2, ny // 4, ny // 2 + ny // 4)
        for (x, y) in product(pixel_x_choice, pixel_y_choice):
            plt.plot(wavs, cube[:, y, x], lw=1, label=str((x,y)))

        plt.xlabel("wavelength (micron)")
        plt.ylabel("pixel (MJy / sr)")
        plt.title(f"A few SEDs from {name_in_title}")
