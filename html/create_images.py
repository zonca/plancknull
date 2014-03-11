import os
import exceptions
import matplotlib.pyplot as plt
import matplotlib
import json
from glob import glob
import sys
sys.path.append("..")

import healpy as hp

root_folder = "../osgtv_10deg_dstfull/"
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
                plot_range = 20
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
                    if int(metadata["channel"]) >= 545:
                        if comp in "QU":
                            plot_range = 1e6
                        else:
                            plot_range = 500
                    if int(metadata["channel"]) >= 857:
                        plot_range = 1e5
                        if comp in "QU":
                            plot_range = 1e6
                    if test_type == "surveydiff" and int(metadata["channel"]) >= 545 and comp == "Q":
                        plot_range *= 1e2
            except exceptions.ValueError:
                pass

        fig = plt.figure(figsize=(9, 6), dpi=100)
        matplotlib.rcParams.update({'font.size': 14})
        hp.mollview(m * 1e6, min=-plot_range, max=plot_range, unit="uK", title=metadata["title"] + " %s" % comp, xsize=900, hold=True)
        plt.savefig(os.path.join(out_folder, metadata["file_name"].replace(".fits", "_%s.jpg" % comp)), dpi=100)
        plt.close()
        fig = plt.figure(figsize=(9, 6), dpi=20)
        fig.add_axes([0.01, 0.01, 0.98, 0.98])
        matplotlib.rcParams.update({'font.size': 30})
        hp.mollview(m * 1e6, min=-plot_range, max=plot_range, cbar=True, title="", xsize=180, hold=True)
        plt.savefig(os.path.join(out_folder, metadata["file_name"].replace(".fits", "_%s_thumb.jpg" % comp)), dpi=20)
        plt.close()

for fold in ["halfrings", "surveydiff", "chdiff"]:
    try:
        os.mkdir(os.path.join(out_folder, fold))
    except:
        pass

for f in sorted(glob(os.path.join(root_folder, "*", "*"  + "*map.json"))):
    print f
    plot_figure(json.load(open(f)))
