import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm
import itertools as it
from astropy.table import Table
from scipy import ndimage
from astropy.wcs import WCS


class MapCollection:
    """Object representing a collection of nd named maps with the same wcs."""

    def __init__(self, keys, shape):
        """
        keys: list of str

        shape: nx, ny
        """
        self.shape = shape
        self.data = np.zeros((len(keys), *self.shape))
        self.index = {k: i for i, k in enumerate(keys)}

    def remove_zero_maps(self):
        """Remove all maps that contain only zeros. Destructive.

        Typically done after running a fit, and the data has been filled
        in. After this operation, to get all maps back, you will need to
        make a new object.
        """
        keys_to_remove = [
            k for k, i in self.index.items() if np.count_nonzero(self.data[i]) == 0
        ]
        self.remove_unused_maps(keys_to_remove)

    def remove_unused_maps(self, keys):
        """Remove maps with the given keys.

        Redundant keys to be removed are allowed, meaning that if given
        key is not in map collection, this is safe and nothing will
        happen. This not was written here, because this makes
        implementing some things easier.

        """
        new_index = {}
        keep_map = np.full(self.data.shape[0], True)

        j = 0
        for k, i in self.index.items():
            # remove if key was specified, or if all values are zero
            if k in keys:
                keep_map[i] = False
            else:
                new_index[k] = j
                j += 1

        # keep only nonzero slices
        self.index = new_index
        self.data = self.data[keep_map]

    def keys(self):
        return self.index.keys()

    def __getitem__(self, name):
        """Access one of the maps using [name]"""
        if name in self.index:
            return self.data[self.index[name]]
        else:
            raise RuntimeError(
                f"Map of {name} not in map collection\n"
                "Available keys are" + str([k for k in self.index])
            )

    def __setitem__(self, name, value):
        if name in self.index:
            self.data[self.index[name]] = value
        else:
            raise RuntimeError(
                f"Map of {name} not in map collection\n"
                "Available keys are" + str([k for k in self.index])
            )

    def plot_map_collage(self, keys, nrows=1, titles=None, colorbar=False, **kwargs):
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

        keys: list of str
            Maps to plot. Might fail if these maps have mostly zeros.

        nrows: int
            How many rows to use. Might be useful for many plots

        titles: list of str
            Title for each panel, to make it clear what each panel means.

        colorbar: bool
            Add a colorbar. WARNING: will use values of the last plot.
            To ensure the colorbar numbers are correct, use kwargs to
            set the rescale or (TODO) vmin/vmax options for plot_map().

        Returns
        -------
        fig, axes: figure and axes created by calling plt.subplots()
        """
        # find the ideal number of columns
        nplots = len(keys)
        q, r = divmod(nplots, nrows)
        ncols = q if r == 0 else q + 1

        fig, axs = plt.subplots(nrows, ncols, sharex=True, sharey=True, squeeze=False)
        kwargs.setdefault("with_title", False)

        allaxes = axs.flatten()

        biggest_x = 0
        biggest_y = 0
        ims = []
        for i, (n, ax) in enumerate(zip(keys, allaxes)):
            plot_info = self.plot_map(n, axes=ax, **kwargs)
            ims.append(plot_info["imshow_return"])
            biggest_x = max(plot_info["image"].shape[1], biggest_x)
            biggest_y = max(plot_info["image"].shape[0], biggest_y)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            if titles is not None:
                ax.set_title(titles[i], rotation=90)

        # because of the shared axes, any imshow command will result
        # in the FOV of the last plot. Make sure that we use the
        # BIGGEST plot.
        axs[0, 0].set_xlim(0, biggest_x)
        axs[0, 0].set_ylim(0, biggest_y)

        fig.tight_layout()
        fig.subplots_adjust(wspace=0, left=0, right=1)
        if colorbar:
            cax = axs[-1, -1].inset_axes([1.04, 0.2, 0.05, 0.6])
            fig.colorbar(ims[-1], ax=axs[-1, -1], cax=cax)

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
        rescale=False,
    ):
        """
        Plot map with some commonly useful options.

        Some of this is definitely overengineered.

        Parameters
        ----------

        name: str
            Key of the map to plot

        log_color: bool
            Use log scale for color map

        axes: Axes
            If None, create axes using plt.figure(). Else use the provided axes.

        with_title: bool
            Use map name as title

        rotate_angle: float
            Apply rotation angle in degrees to image array before
            plotting.

        autocrop: bool
            Automatically crop off rows and columns at the edges that
            are zero. Context: After rotating, the array is expanded to
            fit the rotated rectangle. The empty spaces are filled with
            zeros.

        manualcrop: [int, int, int, int]
            Slice the (optinally rotated) array using the given range
            [xmin, xmax, ymin, ymax]. If autocrop is activated, a set of
            numbers will be printed, which can be copy-pasted to be used
            as this argument to reuse the same cropping later.

        rescale: bool
            Rescale and shift the map values so that they span from 0 to
            1, where the latter number represents the linear position
            between the 1st and 99th percentiles. Useful when a shared
            colorbar is used between different plots. TODO: implement
            configurable pmin, pmax.

        """
        if np.count_nonzero(self[name]) == 0:
            raise ValueError(
                f"Map with name {name} has all zeros and can therefore not be plotted"
            )

        if axes is None:
            plt.figure()
            ax = plt.gca()
        else:
            ax = axes

        map_data = self[name].copy()

        pmin = 1
        pmax = 99
        vmin = np.nanpercentile(map_data, pmin)
        vmax = np.nanpercentile(map_data, pmax)
        kwargs = {"origin": "lower"}

        if rescale:
            # map (vmin vmax) to (0 1)
            map_data = (map_data - vmin) / (vmax - vmin)
            # now map (0 1) to (pmin pmax)
            map_data = (1 - map_data) * pmin + map_data * pmax
            vmin = pmin
            vmax = pmax

        if log_color:
            vmin = np.amin(map_data[map_data > 0])
            kwargs["norm"] = LogNorm(vmin=vmin, vmax=vmax)
        else:
            kwargs["vmin"] = vmin
            kwargs["vmax"] = vmax

        # imshow assumes (y, x), so transpose
        image_data = map_data.T
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

        cropped = False
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
                cropped = True
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
            cropped = True

        if cropped:
            image_data = image_data[min_i:max_i, min_j:max_j]
            crop_translate_xy = np.array([-min_j, -min_i])
        else:
            crop_translate_xy = np.zeros(2)

        im = ax.imshow(image_data, **kwargs)
        if colorbar:
            cax = ax.inset_axes([1.04, 0.2, 0.05, 0.6])
            ax.figure.colorbar(im, ax=ax, cax=cax)

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

        plot_info = dict(
            transform=transform_imagexy_to_mapxy, image=image_data, imshow_return=im
        )
        return plot_info

    def save(self, wcs: WCS, fits_fn, transpose=True, meta=None):
        """Save to one (or many?) files.

        Parameters
        ----------

        wcs : WCS
            celestial WCS that matches the last two dimensions of the data shape

        fits_fn : str
            file name ending in ".fits"

        transpose : bool
            Swap the spatial dimensions of the data. Use when the WCS
            has the y direction as axis 0, and the x direction as axis
            1.

        meta : dict {str: value}
            Metadata to add as extra header entries in the fits file.
            Specified as dictionary {HEADER_KEY: HEADER_VALUE, ...}

        """
        header = wcs.to_header()
        # add extra header keywords if requested
        if meta is not None:
            for k, v in meta.items():
                header[k] = v

        new_hdul = fits.HDUList()
        for k, i in self.index.items():
            data_slice = self.data[i]
            image_array = data_slice.T if transpose else data_slice
            hdu = fits.ImageHDU(data=image_array, header=header, name=k)
            new_hdul.append(hdu)

        new_hdul.writeto(fits_fn, overwrite=True)

    @classmethod
    def load(cls, fits_fn, transpose=True):
        """Load MapCollection object from fits file

        Easier way to recover the maps, compared to loading in all the
        fit results again using a CubeModel object. This way, it is
        straightforward to make the collage overview plot based on one
        file.

        Parameters
        ----------

        fits_fn : str
            file name ending in '.fits'. File is in same format as
            output of 'save()'

        transpose : bool
            Whether to transpose the images contained in the file upon
            loading. True by default, as fits prefers yx while we use xy
            (which is the specutils convention)

        """
        # Look at save function, and try to do the reverse
        with fits.open(fits_fn) as hdul:
            #  skip the primary hdu at '0'
            keys = [hdu.name for hdu in hdul[1:]]
            shape = hdul[1].shape[::-1] if transpose else hdul[1].shape
            instance = cls(keys, shape)

            # copy all image hdus
            for i in range(1, len(hdul)):
                hdu = hdul[i]
                instance[hdu.name] = hdu.data.T if transpose else hdu.data

        return instance

    def save_as_table(
        self, fn, wcs: WCS = None, ignore_patterns=None, **table_write_kwargs
    ):
        """Save the map as an astropy table.

        This way, it becomes easy to explore correlations and their
        relationship to the spatial distribution, with tools such as
        Glue or Topcat. The table contains one row per pixel, has the
        following columns:

        x, y, <all map values>

        RA and DEC columns optional by providing spatial wcs (not
        implemented yet)

        ignore_patterns: list of str
            Maps of which the name contains any of the given substrings
            will be ignored.

        """
        if ignore_patterns is None:
            keep_keys = [name for name in self.index]
        else:
            keep_keys = [
                name
                for name in self.index
                if not any(pattern in name for pattern in ignore_patterns)
            ]

        nx, ny = self.data.shape[1:]
        num_rows = nx * ny
        keys = ["x", "y"]
        if wcs is not None:
            keys += ["ra", "dec"]
        start_rest = len(keys)
        keys += keep_keys
        table_data = np.zeros((num_rows, len(keys)))

        # x and y
        for i, (x, y) in enumerate(it.product(range(nx), range(ny))):
            table_data[i, 0] = x
            table_data[i, 1] = y

        # RA and Dec if wcs is given
        if wcs is not None:
            ra, dec = wcs.pixel_to_world_values(table_data[:, 0], table_data[:, 1])
            table_data[:, 2] = ra
            table_data[:, 3] = dec

        # one column for every map
        for j, map_name in enumerate(keys[start_rest:], start=start_rest):
            m = self[map_name]
            for i, (x, y) in enumerate(it.product(range(nx), range(ny))):
                table_data[i, j] = m[x, y]

        # do not write row when all map values are zero
        nonzero = np.any(table_data[:, 2:], axis=1)
        n = np.count_nonzero(nonzero)
        print(f"{n}/{len(nonzero)} rows are nonzero")
        t = Table(data=table_data[nonzero], names=keys)
        t.write(fn, **table_write_kwargs)


def merge(mc1, mc2):
    if mc1.shape != mc2.shape:
        raise ValueError(
            f"MapCollections have shapes {mc1.shape} and {mc2.shape} but should be equal"
        )

    mc_new = MapCollection([n for n in it.chain(mc1.keys, mc2.keys)], mc1.shape)
    for n in mc1.keys:
        mc_new[n] = mc1[n]
    for n in mc2.keys:
        mc_new[n] = mc2[n]

    return mc_new
