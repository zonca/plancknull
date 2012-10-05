import os
import exceptions
import matplotlib.pyplot as plt
import json
from glob import glob
import sys
sys.path.append("..")

import healpy as hp

root_folder = "../dx8_10deg"
out_folder = "../dx9null/images"

try:
    os.mkdir(out_folder)
except:
    pass

def plot_figure(metadata):
    try:
        allmap = hp.ma(hp.read_map(os.path.join(root_folder, metadata["file_name"]), (0,1,2)))
    except exceptions.IndexError:
        allmap = [hp.ma(hp.read_map(os.path.join(root_folder, metadata["file_name"])))]
    for comp, m in zip("IQU", allmap):
        if comp in "QU":
            plot_range = 20
        else:
            plot_range = 20
        if len(allmap) == 1: #only T, single ch
            plot_range = 20
        is_single_channel = isinstance(metadata["channel"], basestring) and len(metadata["channel"])==6
        if is_single_channel:
            if int(metadata["channel"][3:5]) < 24: # 70GHz
                plot_range = 100
        test_type = metadata["base_file_name"].split("/")[0]
        if isinstance(metadata["channel"], list):
            metadata["channel"] = "_".join(metadata["channel"])
        if metadata["channel"].find("_") < 0:
            try:
                if int(metadata["channel"]) > 70:
                    plot_range = 5
                    if int(metadata["channel"]) >= 353:
                        if comp in "QU":
                            plot_range = 100
                        else:
                            plot_range = 30
                    if int(metadata["channel"]) >= 545 and comp in "QU":
                        plot_range = 1e6
                    if int(metadata["channel"]) >= 857:
                        plot_range = 1e4
                        if comp in "QU":
                            plot_range = 1e6
                    if test_type == "surveydiff" and int(metadata["channel"]) >= 545 and comp == "Q":
                        plot_range *= 1e2
            except exceptions.ValueError:
                pass

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

for f in glob(os.path.join(root_folder, "*", "*SS5-SS4*map.json")):
    print f
    plot_figure(json.load(open(f)))
