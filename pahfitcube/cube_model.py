from itertools import product
from specutils import Spectrum1D
from pahfit.model import Model
import numpy as np
from os.path import isfile
import astropy
from multiprocess.pool import Pool
from tqdm import tqdm
from pathlib import Path
import re

# from astropy import units as u

from pahfitcube.map_collection import MapCollection


def unique_feature_dict(model: Model):
    """Flatten features table into dict of {unique_feature_name: value}"""
    d = {}
    for row in model.features:
        base_name = row["name"]
        for col in ["temperature", "tau", "wavelength", "power", "fwhm"]:
            # do not add unused parameters. But make exception for line FWHM.
            if (
                not (row[col] is np.ma.masked)
                or row["kind"] == "line"
                and col == "fwhm"
            ):
                d[f"{base_name}_{col}"] = row[col][0]
    return d


class CubeModel:
    """The shape of the map collection is determined by the spectral
    cube (spatial dimensions) and the PAHFIT model (number of
    parameters). It therefore makes sense to combine these in one
    object."""

    def __init__(self, model: Model):
        self.init_model = model
        self.flat_feature_names = unique_feature_dict(model).keys()
        self.maps = None
        self.models = {}
        self.prefix = None

    @classmethod
    def load(cls, nx, ny, prefix):
        """Load the files already there, ignoring pixels that still need fitting.

        Shape equal to the original fit must be provided for
        consistency (hint: if you have the spectrum, use
        *spec.shape[:2]

        """
        # list all files matching the prefix
        files = sorted(str(p) for p in Path(".").glob(prefix + "_xy_*.ecsv"))

        # set up based on first file
        model0 = Model.from_saved(files[0])
        instance = cls(model0)
        instance.maps = MapCollection(instance.flat_feature_names, (nx, ny))

        # ingest all models
        for fn in tqdm(files):
            model = Model.from_saved(fn)
            # determine x and y based on file name (TODO: put this in
            # the meta of the saved table)
            m = re.match(f"{prefix}_xy_([0-9]+)_([0-9]+).*?", fn)
            x = int(m[1])
            y = int(m[2])
            instance._ingest_single_model(model, x, y)

        # prefix is remembered so it can be included in metadata when
        # results are exported to file
        instance.prefix = prefix
        return instance

    def fit(
        self,
        cube: Spectrum1D,
        prefix=None,
        maxiter=1000,
        j=1,
        pahfit_guess=False,
    ):
        """Fit the same PAHFIT model to each spaxel of a given cube.

        prefix : str
            e.g. "path/to/dir/prefix_"

            Where progress files are stored, so fitting can be
            interrupted and continued later.

            Files that are already there (and match the pattern) will be
            loaded instead of fitted. This way, fit also acts as a
            loading function.

        pahfit_guess : bool
            If True, use the guess() method of the PAHFIT Model (which
            means the initial guess given to this CubeModel at setup
            will be ignored).

        """
        # I'm being explicit here as a reminder that spectrum1D works
        # with nx, ny, nw! Note that a FITS HDU has nw, ny, nx!
        nx, ny = cube.shape[:-1]
        self.maps = MapCollection(self.flat_feature_names, (nx, ny))

        # Generator for arguments that will be passed to parallel
        # wrapper function. Need to unwrap a lot of things to avoid
        # unpicklable things.
        args_it = (
            dict(
                x=x,
                y=y,
                model=self.init_model,
                prefix=prefix,
                maxiter=maxiter,
                pahfit_guess=pahfit_guess,
                spectral_axis=cube.spectral_axis,
                flux=cube.flux[x, y],
                uncertainty=cube.uncertainty[x, y],
                mask=cube.mask[x, y],
                # meta can contain unpicklable things, so unpack the necessary parts here
                instrument=cube.meta["instrument"],
                # header=cube.meta["header"] if "header" in cube.meta else None,
            )
            for x, y in product(range(nx), range(ny))
        )

        # Same loop for serial and parallel, just use different iterator
        def fit_loop(results_it):
            for x, y, model_xy in tqdm(results_it, total=nx * ny):
                # print(x, y)
                self._ingest_single_model(model_xy, x, y)

        if j > 1:
            with Pool(j) as p:
                # version with parallel Pool.imap
                fit_loop(p.imap(wrapper, args_it, chunksize=1))

        else:
            # version with regular loop
            fit_loop((wrapper(args) for args in args_it))

    def _ingest_single_model(self, model, x, y):
        # print("ingesting result")
        # print(model)
        if model is not None:
            self.models[(x, y)] = model
            flat_result = unique_feature_dict(model)
            for k, v in flat_result.items():
                self.maps[k][x, y] = v

    def get_spaxel_model(self, x, y):
        """Get the individual PAHFIT model for a spaxel

        Returns
        -------
        pahfit.Model fit to spaxel x, y, or None if no model was fit to
        this spaxel.

        """
        k = (x, y)
        if k in self.models:
            return self.models[k]
        else:
            return None

    def plot_spaxel(self, cube, x, y, **kwargs):
        """Default PAHFIT plot for one of the spaxels

        cube : the spectrum to show

        x, y : the xy coordinates in the cube / cube model. Should be
            consistent with what was passed to fit()

        kwargs : options passed to pahfit.model.Model.plot()

        """
        # try:
        fig = self.models[(x, y)].plot(cube[x, y], **kwargs)
        return fig
        # except ValueError as error:
        #     print(error)
        #     print(f"Skipping plot x{x}y{y} due to the above error")
        #     return None

    def make_derived_maps(self, inst, z, dust_cont_wavelength):
        """Generate a number of specialty maps.

        Some examples
        - AIB: sum of all the dust features
        - dust_continuum: sum of all dust continua

        Parameters
        ----------

        spec : Spectrum1D
            the spectral data is needed because we need the
            instrument and the redshift to evaluate the model
            function...

        dust_cont_wavelength : Quantity
            wavelength at which the dust continuum should be evaluated

        Returns
        -------

        MapCollection

        """
        mc = MapCollection(["dust_continuum"], shape=self.maps.shape)

        for (x, y), model in self.models.items():
            tab_spec = model.tabulate(
                inst,
                z,
                [dust_cont_wavelength],
                feature_mask=model.features["kind"] == "dust_continuum",
            )
            dust_continuum_value = tab_spec.flux[0].value
            mc["dust_continuum"][x, y] = dust_continuum_value

        return mc

    def remove_unused_components(self, zero_fraction):
        """Removes maps of all parameters for component's whose power is zero in too many pixels.

        This utility makes the results more readable, by removing
        clutter in the MapCollection object stored at self.maps.

        zero_fraction: if the fraction of zero pixels in a power map,
        above which a map will be flagged for removal

        """
        keys_to_remove = []
        # go over all the component names, and check if that component
        # is mostly zero everywhere. If so, remove the appropriate keys
        # relating to that component (name_wavelength, name_fwhm, etc)
        for name, kind in zip(
            self.init_model.features["name"], self.init_model.features["kind"]
        ):
            if (
                kind == "line" or kind == "dust_feature"
            ) and name + "_power" in self.maps.keys():
                power_map = self.maps[name + "_power"]
                if np.count_nonzero(power_map) < zero_fraction * np.size(power_map):
                    keys_to_remove.append(name + "_wavelength")
                    keys_to_remove.append(name + "_power")
                    keys_to_remove.append(name + "_fwhm")

        self.maps.remove_unused_maps(keys_to_remove)

    def export_maps(self, wcs, fits_fn):
        """Save fit results as maps

        Wrapper around MapCollection.save(). Will include the location
        of the PAHFIT models in the fits header as a reminder

        Parameters
        ----------

        wcs : WCS
            celestial WCS that matches the last two dimensions of the data shape

        fits_fn : str
            file name ending in ".fits"

        """
        self.maps.save(wcs, fits_fn, meta={"PAHFITSD": self.prefix})


def _skip(spec):
    num_wav = len(spec.wavelength)
    too_many_zeros = np.count_nonzero(spec.flux.value <= 0) > 0.2 * num_wav
    too_many_nan = np.count_nonzero(~np.isfinite(spec.flux.value)) > 0.2 * num_wav
    return too_many_zeros or too_many_nan


def _result_filename(x, y, prefix):
    if prefix is None:
        return None

    return f"{prefix}_xy_{x}_{y}.ecsv"


def _load_fit_save(x, y, spec, model: Model, maxiter, prefix, pahfit_guess):
    model_xy = None

    # shortcut: if model file already exists at given directory, just
    # load the model and return it
    fn = _result_filename(x, y, prefix)
    if fn is not None and isfile(fn):
        return Model.from_saved(fn)

    # analyze spectrum and decide if we will skip or not. Can take
    # about 1s for big spectra, so do this only if we cannot load.
    if _skip(spec):
        print(f"Skipping ({x, y}), too many bad values")
        return None

    try:
        # copy so we don't overwrite the initial conditions
        model_xy = Model(model.features.copy())

        if pahfit_guess:
            model_xy.guess(spec)

        try:
            model_xy.fit(spec, maxiter=maxiter, verbose=False)
        except TypeError as e:
            print(f"caught {e} in cube_model fit")
            print("x, y = ", (x, y))
            print("spec = ", spec)
            raise e
        # save result if prefix was provided
        if fn is not None:
            model_xy.save(fn)
    except astropy.modeling.fitting.NonFiniteValueError as e:
        print(f"Fit failed for pixel {x} {y} due to the following Exception: {e}")

    return model_xy


def wrapper(args):
    """Call to the fit function, compatible with pickling.

    Should pack the following in to args:

    "x": int

    "y": int

    "model": the pahfit Model

    "prefix": str

    "maxiter": int

    "pahfit_guess": bool

    'Unpacked' version of the spectral cube
    ---------------------------------------
    "flux": array

    "uncertainty": StdDevUncertainty

    "mask": mask, if any, to make the fit ignore points

    "instrument": str
        The contents of spec.meta["instrument"]

    "header": ?
        The contents of spec.meta["header"]
    """
    x = args["x"]
    y = args["y"]
    model = args["model"]
    prefix = args["prefix"]
    maxiter = args["maxiter"]
    pahfit_guess = args["pahfit_guess"]

    # Spectrum1D not picklable at the moment, so recreate it here (i.e.,
    # on the parallel process)
    spectral_axis = args["spectral_axis"]
    flux = args["flux"]
    uncertainty = args["uncertainty"]
    mask = args["mask"]
    spec = Spectrum1D(flux, spectral_axis, uncertainty=uncertainty)
    if mask is not None:
        spec.mask = mask
    spec.meta = {k: args[k] for k in ("header", "instrument") if k in args}

    return (
        x,
        y,
        _load_fit_save(x, y, spec, model, maxiter, prefix, pahfit_guess=pahfit_guess),
    )
