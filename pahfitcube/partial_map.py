"""Tool to load all models, and inspect output which is not done yet as a function of x and y
(e.g. because the code stopped somewhere)

"""

from argparse import ArgumentParser
from pahfit.model import Model
import re
import numpy as np
from matplotlib import pyplot as plt

ap = ArgumentParser()
ap.add_argument("fit_result_files", nargs="+")
args = ap.parse_args()

# open every fit result
xs = []
ys = []
models = []

for f in args.fit_result_files:
    # determine x and y coords from filename
    m = re.match(".*?xy_([0-9]+)_([0-9]+).*?", f)
    x = int(m[1])
    y = int(m[2])
    model = Model.from_saved(f)

    xs.append(x)
    ys.append(y)
    models.append(model)

# plot something here
nx = max(xs) + 1
ny = max(ys) + 1
map_array = np.zeros((ny, nx))

for x, y, model in zip(xs, ys, models):
    map_array[y, x] = model.features.loc["PAH_3.4"]["power"][0]

plt.imshow(map_array)
plt.colorbar()
plt.show()
