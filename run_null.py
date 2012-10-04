import numpy as np
import healpy as hp
from glob import glob
import logging as log
import os
from differences import halfrings, surveydiff, chdiff 
import reader
from ConfigParser import SafeConfigParser

import utils
import sys

if len(sys.argv) < 2:
    print "Launch script as: python run_null.py run_*.conf"
    sys.exit(1)

#read configuration
config = SafeConfigParser()
config.read(sys.argv[1])

paral = config.getboolean("paral")
root_folder = config.get("output_folder")
try:
    os.mkdir(root_folder)
except:
    pass
if paral:
    from IPython.parallel import Client


log.root.level = log.DEBUG

# create map reader
mapreader = reader.DXReader(config.get("run", "reader_conf"), nside=config.getint("run", "nside"))
smooth_combine_config = dict(fwhm=np.radians(config.getfloat("smooth_combine", "smoothing")), degraded_nside=config.getint("smooth_combine", "degraded_nside"), spectra=True)

if paral:
    tasks = []
    tc = Client()
    lview = tc.load_balanced_view() # default load-balanced view

if config.getboolean("run", "run_halfrings"):
    print "HALFRINGS"
    survs = ["nominal", "full"]
    freqs = [30,44,70]
    for freq in freqs:
        chtags = [""]
        if freq == 70:
            chtags += ["18_23", "19_22", "20_21"]
        for chtag in chtags:
            for surv in survs:
                if paral:
                    tasks.append(lview.apply_async(halfrings,freq, chtag,
                                                   surv, pol='IQU',
                                                   smooth_combine_config=smooth_combine_config,
                                                   root_folder=root_folder,log_to_file=True,
                                                   mapreader=mapreader))
                else:
                    halfrings(freq, chtag, surv, pol='IQU',
                              smooth_combine_config=smooth_combine_config,
                              root_folder=root_folder,log_to_file=False,
                             mapreader=mapreader)

if config.getboolean("run", "run_surveydiff"):
    print "SURVDIFF"
    survs = [1,2,3,4,5]
    freqs = [30, 44, 70]
    for bp_corr in [False, True]:
        for freq in freqs:
            chtags = [""]
            #chtags = []; log.warning("Disabled full freq")
            if freq == 70:
                chtags += ["18_23", "19_22", "20_21"]
            #log.warning("Disabled quadruplets")
            chtags += utils.chlist(freq)
            for chtag in chtags:
                if bp_corr and chtag: # no corr for single ch
                    continue
                if chtag and chtag.find("_") < 0:
                    pol='I'
                else:
                    pol="IQU"
                if paral:
                    tasks.append(lview.apply_async(surveydiff,freq, chtag,
                                                   survs, pol=pol,
                                                   smooth_combine_config=smooth_combine_config,
                                                   root_folder=root_folder,log_to_file=True,
                                                   bp_corr=bp_corr,
                                                   mapreader=mapreader))
                else:
                    surveydiff(freq, chtag, survs, pol=pol,
                               smooth_combine_config=smooth_combine_config,
                               root_folder=root_folder,log_to_file=False,
                               bp_corr=bp_corr, mapreader=mapreader)

if config.getboolean("run", "run_chdiff"):
    print "CHDIFF"
    survs = [1,2,3,4,5]
    freqs = [30, 44, 70]
        
    for freq in freqs:
        for surv in survs:
            if paral:
               tasks.append(lview.apply_async(chdiff,freq, ["LFI%d" % h for
                                                            h in
                                                            utils.HORNS[freq]],
                                              surv, pol='I',
                                              smooth_combine_config=smooth_combine_config,
                                              root_folder=root_folder,
                                              log_to_file=True,
                                              mapreader=mapreader))
            else:
               chdiff(freq, ["LFI%d" % h for h in utils.HORNS[freq]], surv,
                      pol='I', smooth_combine_config=smooth_combine_config,
                      root_folder=root_folder,
                      log_to_file=False,
                      mapreader=mapreader)

if config.getboolean("run", "compute_union_mask"):
    assert ~paral
    print "UNIONMASK"
    survs = [1,2,3,4,5]
    freqs = [30, 44, 70]
    for freq in freqs:
        mask = utils.read_mask(
            glob(os.path.join(os.environ["DDX9_LFI"], "MASKs",
                          'destripingmask_%d.fits' % freq))[-1]
                         , nside= )
        mask |= utils.read_mask("/global/project/projectdirs/planck/user/zonca/masks/wmap_polarization_analysis_mask_r9_7yr_v4.fits", nside=NSIDE)
        chtags = [""]
        if freq == 70:
            chtags += ["18_23", "19_22", "20_21"]
        for chtag in chtags:
            for surv in survs:
                mask |= mapreader(freq, surv, chtag, halfring=0, pol='Q').mask
        hp.write_map("union_mask_%d.fits" % freq, mask)

if paral:
    print("Wait for %d tasks to complete" % len(tasks))
    tc.wait(tasks)
