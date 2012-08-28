import os
import exceptions
import gc
import matplotlib.pyplot as plt
import json
from glob import glob

import healpy as hp

root_folder = "../dx9"
out_folder = "../dx9null/images"

try:
    os.mkdir(out_folder)
except:
    pass

print "HALFRINGS"

def plot_figure(metadata):
    try:
        allmap = hp.ma(hp.read_map(os.path.join(root_folder, metadata["file_name"]), (0,1,2)))
    except exceptions.IndexError:
        allmap = [hp.ma(hp.read_map(os.path.join(root_folder, metadata["file_name"])))]
    for comp, m in zip("IQU", allmap):
        if comp in "QU":
            plot_range = 30
        else:
            plot_range = 30
        fig = plt.figure(figsize=(9, 6), dpi=100)
        hp.mollview(m * 1e6, min=-plot_range, max=plot_range, unit="muK", title=metadata["title"] + " %s" % comp, xsize=900, hold=True)
        plt.savefig(os.path.join(out_folder, metadata["file_name"].replace(".fits", "_%s.jpg" % comp)), dpi=100)
        plt.close()
        fig = plt.figure(figsize=(9, 6), dpi=20)
        fig.add_axes([0.01, 0.01, 0.98, 0.98])
        hp.mollview(m * 1e6, min=-plot_range, max=plot_range, cbar=False, title="", xsize=180, hold=True)
        plt.savefig(os.path.join(out_folder, metadata["file_name"].replace(".fits", "_%s_thumb.jpg" % comp)), dpi=20)
        plt.close()

for fold in ["halfrings", "surveydiff", "chdiff"]:
    try:
        os.mkdir(os.path.join(out_folder, fold))
    except:
        pass

for f in glob(os.path.join(root_folder, "*", "*map.json")):
    print f
    plot_figure(json.load(open(f)))
