from itertools import product
from specutils import Spectrum1D
from pahfit.model import Model
import numpy as np
from os.path import isfile
from astropy.io import fits
from matplotlib import pyplot as plt
import astropy


class MapCollection:
    """Object representing a collection of nd named maps with the same wcs."""

    def __init__(self, names, shape):
        self.names = names
        self.data = np.zeros((len(names), *shape))
        self.index = {name: i for i, name in enumerate(self.names)}

    def __getitem__(self, name):
        """Access one of the maps using [name]"""
        if name in self.names:
            return self.data[self.index[name]]
        else:
            raise RuntimeError(f"Map of {name} not in map collection\n"
                               "Available names are" + str([k for k in self.index]))

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
            if not row[col].mask[0] or row["kind"] == 'line' and col == "fwhm":
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

    def fit(self, cube: Spectrum1D, checkpoint_prefix=None, maxiter=1000):
        """The main fitting routine. Fits the same PAHFIT model to each
        spaxel of the given cube

        checkpoint_prefix : str
            Where progress files are stored, so fitting can be
            interrupted and continued later.
        """
        self.maps = MapCollection(self.flat_feature_names, cube.shape[:-1])
        self.models = {}

        # make many copies of model and fit them
        nx, ny = cube.shape[:-1]
        for x, y in product(range(nx), range(ny)):
            if all(cube.flux[x, y] == 0) or not all(np.isfinite(cube.flux[x, y])):
                continue
            model_xy = None

            try:
                # with checkpoints
                if checkpoint_prefix is not None:
                    fn = f"{checkpoint_prefix}_xy_{x}_{y}.ecsv"
                    if isfile(fn):
                        model_xy = Model.from_saved(fn)
                    else:
                        model_xy = self.model.copy()
                        model_xy.fit(cube[x, y], maxiter=maxiter)
                        model_xy.save(fn)

                # without checkpoints, use simpler logic
                else:
                    model_xy = self.model.copy()
                    model_xy.fit(cube[x, y], maxiter=maxiter)
                # store all the models for later reference
                self.models[(x, y)] = model_xy
                # also write to mapcollection
                self._ingest_single_model(model_xy, x, y)

            except astropy.modeling.fitting.NonFiniteValueError as e:
                print(f"Fit failed for pixel {x} {y} due to {e}")

    def _ingest_single_model(self, model, x, y):
        # print("ingesting result")
        # print(model)
        flat_result = unique_feature_dict(model)
        for k, v in flat_result.items():
            self.maps[k][x, y] = v

    def plot_spaxel(self, cube, x, y, save_fn=None):
        plt.figure()
        try:
            fig, axs = self.models[(x, y)].plot(cube[x, y])
        except ValueError as error:
            print(error)
            print(f"Skipping plot x{x}y{y} due to the above error")
            # raise error
            # continue

        if save_fn is not None:
            fig.savefig(save_fn)

        plt.close(fig)
