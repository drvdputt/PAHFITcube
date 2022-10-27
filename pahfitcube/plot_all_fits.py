from argparse import ArgumentParser
from pahfit.model import Model
from specutils import Spectrum1D
import re

ap = ArgumentParser()
ap.add_argument("cube")
ap.add_argument("instrumentname")
ap.add_argument("fit_result_files", nargs="+")
args = ap.parse_args()

spec = Spectrum1D.read(args.cube)
spec.meta["instrument"] = args.instrumentname

# open every fit result, and plot it
for f in args.fit_result_files:
    # determine x and y coords from filename
    m = re.match(".*?xy_([0-9]+)_([0-9]+).*?", f)
    x = int(m[1])
    y = int(m[2])

    # get spectrum for this x and y
    s = spec[y, x]

    # load model and plot
    m = Model.from_saved(f)

    fig = m.plot(s)
    fig.savefig(f + ".pdf")
