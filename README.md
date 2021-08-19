# PAHFIT-cube

Collection of scripts to experiment with running PAHFIT (citation) on IFU data of the ISM. This
is in the context of the ERS PDR program for JWST (http://jwst-ism.org/).

## Test data

We do our experiments and testing with Spitzer IFU cubes from the SAGE-Spec (citation) program.
From this dataset, one HII region was taken, for which the wavelength range from 35 ~ 5 micron
is covered by 4 IFU cubes: LL1, LL2, SL1, and SL2 (in the order from long to short wavelengths).

If you are on the ERS-PDRs Slack space, you can look for the data there. Or ask me for it. The
list of files should look like this

    ls -1 PAHFIT-cube/data/sage-spec_hii1_hii8_4dec08
    hii1_hii8_ll_LL1_cube.fits
    hii1_hii8_ll_LL1_cube_unc.fits
    hii1_hii8_ll_LL2_cube.fits
    hii1_hii8_ll_LL2_cube_unc.fits
    hii1_hii8_ll_LL3_cube.fits
    hii1_hii8_ll_LL3_cube_unc.fits
    hii1_hii8_sl_SL1_cube.fits
    hii1_hii8_sl_SL1_cube_unc.fits
    hii1_hii8_sl_SL2_cube.fits
    hii1_hii8_sl_SL2_cube_unc.fits
    hii1_hii8_sl_SL3_cube.fits
    hii1_hii8_sl_SL3_cube_unc.fits

## Create input cube

Because PAHFIT needs the whole wavelength range to work with, the 4 Spitzer cubes need to be
merged first. The script `resample_and_merge_spitzer_cubes.py` takes care this. Since this is an
early test, it does not have command line arguments yet, so it only works for the data above. Simply run

    python resample_and_merge_spitzer_cubes.py

### Reprojection
Because the field of view is different for each IFU cube, the Python package `reproject` is used
to resample all the slices of LL1, LL2, SL1, and SL2 onto the same spatial grid.

(picture here)

### Spectral Order Stitching
A basic spectral order stitching step has been implemented. For every pair of cubes adjacent in
wavelength space, the spectra are rescaled on a pixel-per-pixel basis, according to the average
values in the wavelength overlap region. This is done in order of decreasing wavelength:
- scale LL2 to match LL1
- scale SL1 to match LL2
- scale SL2 to match SL1

After this step, the slices of the 4 cubes are merged and sorted by wavelength.

## Run PAHFIT-cube
The script `run_pahfit_cube.py` has the same command line arguments as the normal PAHFIT run
script. The reprojected and merged data cube can be used as input, and the Spitzer extragalactic
pack should be used until new calibration settings are developed. Running the script can be done
as follows

    ipython run_pahfit_cube.py reprojected.fits scipack_ExGal_SpitzerIRSSLLL.ipac --fit_maxiter 50

For debugging the script, it is best to set `--fit_maxiter` to a low number, so that everything
goes faster.

### Output map
The main output is a multi-extension fits file `reprojected_parameter_maps.fits`, with the same
WCS as the input spectral cube. Each extension contains a map of one of the fit parameters. It
can be viewed in DS9 when opened as a multi-extension cube.

### Output per pixel and resume option
For every pixel, the normal PAHFIT output is written out (fit parameters and a plot). This is
useful for diagnosing problems with individual pixels.

Currently there's a basic `--resume` option, which skips pixels when their output file is
already present. It loads the fit results for those pixels from the output files instead.

### Multiprocessing
The `-j` option controls the number of processes. Looking at the implementation, one can see
that parallelizing the fitting is relatively simple.
