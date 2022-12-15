import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
import itertools as it
from astropy.table import Table
from scipy import ndimage


class MapCollection:
    """Object representing a collection of nd named maps with the same wcs."""

    def __init__(self, names, shape):
        """
        names: list of str

        shape: nx, ny
        """
        self.names = names
        self.shape = shape
        self.data = np.zeros((len(names), *shape))
        self.index = {name: i for i, name in enumerate(names)}

    def __getitem__(self, name):
        """Access one of the maps using [name]"""
        if name in self.names:
            return self.data[self.index[name]]
        else:
            raise RuntimeError(
                f"Map of {name} not in map collection\n"
                "Available names are" + str([k for k in self.index])
            )

    def __setitem__(self, name, value):
        if name in self.names:
            self.data[self.index[name]] = value
        else:
            raise RuntimeError(
                f"Map of {name} not in map collection\n"
                "Available names are" + str([k for k in self.index])
            )

    def plot_map_collage(self, names, nrows=1, titles=None, **kwargs):
        """
        Plot many maps in a compact way using subplots.

        The main goal is providing a quick overview of the fit results,
        without cluttering everything with labels and numbers.

        The following visualization choices are made:
        - No color bars (the fit parameters don't always have
          straightforward units anyway)
        - No pixel number axis labels. If fancy plots with coordinates
          and color bar units are required, those changes should be made
          in the returned figure/axes.

        self.plot_map is called for each panel

        Parameters
        ----------

        names: list of str
            Maps to plot. Might fail if these maps have mostly zeros.

        Returns
        -------
        fig, axes: figure and axes created by calling plt.subplots()
        """
        # find the ideal number of columns
        nplots = len(names)
        q, r = divmod(nplots, nrows)
        ncols = q if r == 0 else q + 1

        fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True, squeeze=False)
        kwargs.setdefault("with_title", False)

        allaxes = axs.flatten()

        biggest_x = 0
        biggest_y = 0
        for i, (n, ax) in enumerate(zip(names, allaxes)):
            plot_info = self.plot_map(n, axes=ax, **kwargs)
            biggest_x = max(plot_info["image"].shape[1], biggest_x)
            biggest_y = max(plot_info["image"].shape[0], biggest_y)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            if titles is not None:
                ax.set_title(titles[i])

        # because of the shared axes, any imshow command will result
        # in the FOV of the last plot. Make sure that we use the
        # BIGGEST plot.
        axs[0, 0].set_xlim(0, biggest_x)
        axs[0, 0].set_ylim(0, biggest_y)

        fig.subplots_adjust(wspace=0, left=0, right=1)
        return fig, axs

    def plot_map(
        self,
        name,
        log_color=False,
        axes=None,
        with_title=True,
        rotate_angle=None,
        autocrop=False,
        manualcrop=None,
        colorbar=False,
    ):
        if np.count_nonzero(self[name]) == 0:
            raise ValueError(
                f"Map with name {name} has all zeros and can therefore not be plotted"
            )

        """Utility function dealing with the most common plotting issues"""
        if axes is None:
            plt.figure()
            ax = plt.gca()
        else:
            ax = axes

        map_data = self[name]
        vmin = np.amin(map_data)
        vmax = np.percentile(map_data, 99)
        kwargs = {"origin": "lower"}
        if log_color:
            vmin = np.amin(map_data[map_data > 0])
            kwargs["norm"] = LogNorm(vmin=vmin, vmax=vmax)
        else:
            kwargs["vmin"] = vmin
            kwargs["vmax"] = vmax

        # imshow assumes (y, x), so transpose
        image_data = self[name].T
        # rotate if requested
        if rotate_angle is not None:
            # rotating causes some artifacts with very small nonzero
            # values. This messes with the autocrop below. We can cut
            # these values out by remembering the minimum nonzero data
            # value, and setting anything below that to zero in the
            # final image.
            nonzero = image_data > 0
            if nonzero.any():
                vmin_nonzero = np.amin(image_data[image_data > 0])
            else:
                # avoid problems when everything is zero
                vmin_nonzero = 0

            image_data = ndimage.rotate(image_data, rotate_angle)
            image_data[image_data < vmin_nonzero] = 0
            # remember the rotation matrix for this rotation
            radians = rotate_angle * np.pi / 180
            c = np.cos(radians)
            s = np.sin(radians)
            rotation_matrix = np.array([[c, -s], [s, c]])
        else:
            rotation_matrix = np.identity(2)

        # Remember the shape after rotation. Important for calculating
        # the reverse transformation.
        center_after_rotation = np.array(image_data.T.shape) / 2

        do_crop = False
        if autocrop:
            keep_i = np.where(np.sum(np.square(image_data), axis=1) != 0)[0]
            keep_j = np.where(np.sum(np.square(image_data), axis=0) != 0)[0]
            # print("keep rows", keep_i)
            # print("keep cols", keep_j)
            if len(keep_i) > 2 and len(keep_j) > 2:
                min_i = keep_i[0]
                max_i = keep_i[-1]
                min_j = keep_j[0]
                max_j = keep_j[-1]
                do_crop = True
                print("Suggested crop range: ", (min_i, max_i, min_j, max_j))
            else:
                print("Something weird with data! Skipping autocrop")
        elif manualcrop:
            min_i, max_i, min_j, max_j = (
                manualcrop[0],
                manualcrop[1],
                manualcrop[2],
                manualcrop[3],
            )
            do_crop = True

        if do_crop:
            image_data = image_data[min_i:max_i, min_j:max_j]
            crop_translate_xy = np.array([-min_j, -min_i])
        else:
            crop_translate_xy = np.zeros(2)

        cax = ax.imshow(image_data, **kwargs)
        if colorbar:
            ax.figure.colorbar(cax, ax=ax)

        if with_title:
            ax.set_title(name)

        # define a function which can be used to convert the displayed
        # XY to the XY of the original cube
        def transform_imagexy_to_mapxy(x, y):
            new_xy = np.array([x, y])
            print("start", new_xy)
            # undo the translation due to cropping
            undo_crop_xy = new_xy - crop_translate_xy
            print("undo crop", undo_crop_xy)
            # shift from center of frame to origin
            origin_xy = undo_crop_xy - center_after_rotation
            print("to origin", origin_xy)
            # apply inverted rotation matrix
            rot_origin_xy = rotation_matrix.dot(origin_xy)
            print("rotated", rot_origin_xy)
            # shift back to original center position BEFORE rotation
            center_original = np.array(map_data.shape) / 2
            final_xy = rot_origin_xy + center_original
            print("back to center", final_xy)
            return final_xy

        plot_info = dict(transform=transform_imagexy_to_mapxy, image=image_data)
        return plot_info

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


def merge(mc1, mc2):
    if mc1.shape != mc2.shape:
        raise ValueError(
            f"MapCollections have shapes {mc1.shape} and {mc2.shape} but should be equal"
        )

    mc_new = MapCollection([n for n in it.chain(mc1.names, mc2.names)], mc1.shape)
    for n in mc1.names:
        mc_new[n] = mc1[n]
    for n in mc2.names:
        mc_new[n] = mc2[n]

    return mc_new
