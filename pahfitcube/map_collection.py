import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
import itertools as it
from astropy.table import Table


class MapCollection:
    """Object representing a collection of nd named maps with the same wcs."""

    def __init__(self, names, shape):
        """
        names: list of str

        shape: nx, ny
        """
        self.names = names
        self.data = np.zeros((len(names), *shape))
        self.index = {name: i for i, name in enumerate(self.names)}

    def __getitem__(self, name):
        """Access one of the maps using [name]"""
        if name in self.names:
            return self.data[self.index[name]]
        else:
            raise RuntimeError(
                f"Map of {name} not in map collection\n"
                "Available names are" + str([k for k in self.index])
            )

    def plot_map(self, name, log_color=False, axes=None):
        if axes is None:
            plt.figure()
            ax = plt.gca()
        else:
            ax = axes

        array = self[name]
        vmin = np.amin(array)
        vmax = np.percentile(array, 99)
        kwargs = {"origin": "lower"}
        if log_color:
            vmin = np.amin(array[array > 0])
            kwargs["norm"] = LogNorm(vmin=vmin, vmax=vmax)
        else:
            kwargs["vmin"] = vmin
            kwargs["vmax"] = vmax
        # imshow assumes (y, x), so transpose
        cax = ax.imshow(self[name].T, **kwargs)
        ax.figure.colorbar(cax, ax=ax)

    def save(self, wcs, fits_fn):
        """Save to one (or many?) files."""
        header = wcs.to_header()
        new_hdul = fits.HDUList()
        for k in self.index:
            i = self.index[k]
            hdu = fits.ImageHDU(data=self.data[i], header=header, name=k)
            new_hdul.append(hdu)
        new_hdul.writeto(fits_fn, overwrite=True)

    def save_as_table(self, fn, wcs=None, **table_write_kwargs):
        """Save the map as an astropy table.

        This way, it becomes easy to explore the data in glue or topcat
        (e.g. figure out which part of the nebula a certaint grouping in
        a scatter plot belongs to). The table contains one row per
        pixel, has the following columns:

        x, y, <all map values>

        RA and DEC columns optional?

        """
        nx, ny = self.data.shape[1:]
        num_rows = nx * ny
        names = ["x", "y"] + [name for name in self.index]
        table_data = np.zeros((num_rows, len(names)))

        # x and y
        for i, (x, y) in enumerate(it.product(range(nx), range(ny))):
            table_data[i, 0] = x
            table_data[i, 1] = y

        # one column for every map
        for j, map_name in enumerate(names[2:], start=2):
            m = self[map_name]
            for i, (x, y) in enumerate(it.product(range(nx), range(ny))):
                table_data[i, j] = m[x, y]

        # do not write row when all map values are zero
        nonzero = np.any(table_data[:, 2:], axis=1)
        n = np.count_nonzero(nonzero)
        print(f"{n}/{len(nonzero)} rows are nonzero")
        t = Table(data=table_data[nonzero], names=names)
        t.write(fn, **table_write_kwargs)
