from itertools import product
from specutils import Spectrum1D
from pahfit.model import Model
import numpy as np
from os.path import isfile
from astropy.io import fits
from matplotlib import pyplot as plt
import astropy
from multiprocess.pool import Pool
from tqdm import tqdm


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

    def plot_map(self, name):
        plt.figure()
        array = self[name]
        vmin = np.amin(array)
        vmax = np.percentile(array, 99)
        # imshow assumes (y, x), so transpose
        plt.imshow(self[name].T, origin="lower", vmin=vmin, vmax=vmax)
        plt.colorbar()

    def save(self, wcs, fits_fn):
        """Save to one (or many?) files."""
        header = wcs.to_header()
        new_hdul = fits.HDUList()
        for k in self.index:
            i = self.index[k]
            hdu = fits.ImageHDU(data=self.data[i], header=header, name=k)
            new_hdul.append(hdu)
        new_hdul.writeto(fits_fn, overwrite=True)


def unique_feature_dict(model: Model):
    """Flatten features table into dict of {unique_feature_name: value}"""
    d = {}
    for row in model.features:
        base_name = row["name"]
        for col in ["temperature", "tau", "wavelength", "power", "fwhm"]:
            # do not add unused parameters. But make exception for line FWHM.
            if not row[col].mask[0] or row["kind"] == "line" and col == "fwhm":
                d[f"{base_name}_{col}"] = row[col][0]
    return d


class CubeModel:
    """The shape of the map collection is determined by the spectral
    cube (spatial dimensions) and the PAHFIT model (number of
    parameters). It therefore makes sense to combine these in one
    object."""

    def __init__(self, model: Model):
        self.model = model
        self.flat_feature_names = unique_feature_dict(model).keys()
        self.maps = None

    def fit(self, cube: Spectrum1D, checkpoint_prefix=None, maxiter=1000, j=1):
        """Fit the same PAHFIT model to each spaxel of a given cube.

        checkpoint_prefix : str
            e.g. "path/to/dir/prefix_"

            Where progress files are stored, so fitting can be
            interrupted and continued later.

            Files that are already there (and match the pattern) will be
            loaded instead of fitted. This way, fit also acts as a
            loading function.

        """
        # I'm being explicit here as a reminder that spectrum1D works
        # with nx, ny, nw! Note that a FITS HDU has nw, ny, nx!
        nx, ny = cube.shape[:-1]
        self.maps = MapCollection(self.flat_feature_names, (nx, ny))
        self.models = {}

        if j > 1:
            self.model.save("./temp.ecsv", overwrite=True)
            with Pool(j) as p:
                # generates arguments for wrapper
                #
                args_it = (
                    dict(
                        x=x,
                        y=y,
                        # cube=cube,
                        model="temp.ecsv",
                        checkpoint_prefix=checkpoint_prefix,
                        maxiter=maxiter,
                    )
                    for x, y in product(range(nx), range(ny))
                )
                # generates models by calling wrapper on arguments
                results_it = p.imap(wrapper, args_it)

                # do the parallel iteration
                for model_xy, args in zip(results_it, args_it):
                    if model_xy is not None:
                        self._ingest_single_model(model_xy, args["x"], args["y"])

        else:
            for x, y in tqdm(product(range(nx), range(ny))):
                self._fit_spaxel(cube, x, y, maxiter, checkpoint_prefix)

    def _ingest_single_model(self, model, x, y):
        # print("ingesting result")
        # print(model)
        self.models[(x, y)] = model
        flat_result = unique_feature_dict(model)
        for k, v in flat_result.items():
            self.maps[k][x, y] = v

    def _fit_spaxel(self, cube, x, y, maxiter, checkpoint_prefix):
        spec = cube[x, y]
        if self._skip(spec):
            print(f"Skipping ({x, y}), too many bad values")
            return

        model_xy = None
        try:
            # with checkpoints
            if checkpoint_prefix is not None:
                fn = f"{checkpoint_prefix}_xy_{x}_{y}.ecsv"
                if isfile(fn):
                    print(f"Loading {fn}")
                    model_xy = Model.from_saved(fn)
                else:
                    model_xy = self.model.copy()
                    model_xy.fit(spec, maxiter=maxiter)
                    model_xy.save(fn)

            # without checkpoints, use simpler logic
            else:
                model_xy = self.model.copy()
                model_xy.fit(spec, maxiter=maxiter)
            # store all the models for later reference
            self.models[(x, y)] = model_xy
            # also write to mapcollection
            self._ingest_single_model(model_xy, x, y)

        except astropy.modeling.fitting.NonFiniteValueError as e:
            print(f"Fit failed for pixel {x} {y} due to {e}")

    def _skip(self, spec):
        num_wav = len(spec.wavelength)
        too_many_zeros = np.count_nonzero(spec.flux.value <= 0) > 0.2 * num_wav
        has_nan = not all(np.isfinite(spec.flux.value))
        return too_many_zeros or has_nan

    def plot_spaxel(self, cube, x, y):
        try:
            fig = self.models[(x, y)].plot(cube[x, y])
            return fig
        except ValueError as error:
            print(error)
            print(f"Skipping plot x{x}y{y} due to the above error")
            return None


def wrapper(args):
    x = args["x"]
    y = args["y"]
    cube = args["cube"]
    model = args["model"]
    checkpoint_prefix = args["checkpoint_prefix"]
    maxiter = args["maxiter"]

    if all(cube.flux[x, y] == 0) or not all(np.isfinite(cube.flux[x, y])):
        return None

    model_xy = None
    try:
        # with checkpoints
        if checkpoint_prefix is not None:
            fn = f"{checkpoint_prefix}_xy_{x}_{y}.ecsv"
            if isfile(fn):
                model_xy = Model.from_saved(fn)
            else:
                model_xy = Model.from_saved(model)
                model_xy.fit(cube[x, y], maxiter=maxiter)
                model_xy.save(fn)
        # without checkpoints, use simpler logic
        else:
            model_xy = Model.from_saved(model)
            model_xy.fit(cube[x, y], maxiter=maxiter)
        # store all the models for later reference
        return None
    except astropy.modeling.fitting.NonFiniteValueError as e:
        print(f"Fit failed for pixel {x} {y} due to {e}")
        return None
