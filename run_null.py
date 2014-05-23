import numpy as np
import logging as log
import os
from differences import halfrings, surveydiff, chdiff 
import reader
from ConfigParser import SafeConfigParser, NoOptionError
import exceptions
import json

import utils
import sys

if len(sys.argv) < 2:
    print "Launch script as: python run_null.py ,6,7run_*.conf"
    sys.exit(1)

#read configuration
config = SafeConfigParser()
config.read(sys.argv[1])

paral = config.getboolean("run", "paral")
root_folder = config.get("run", "output_folder")
try:
    os.mkdir(root_folder)
except:
    pass
if paral:
    from IPython.parallel import Client

log.root.level = log.DEBUG

# create map reader
mapreader = reader.DXReader(config.get("run", "reader_conf"), nside=config.getint("smooth_combine", "nside"), debug=config.getboolean("run", "debug"))
smooth_combine_config = dict(fwhm=np.radians(config.getfloat("smooth_combine", "smoothing")), degraded_nside=config.getint("smooth_combine", "degraded_nside"), spectra=config.getboolean("smooth_combine", "spectra"), chi2=config.getboolean("smooth_combine", "chi2"))

survs = [1,2,3,4,5,6,7,8,9]
survs = [1,2,3,4,5,6,7,8]

if paral:
    tasks = []
    tc = Client()
    lview = tc.load_balanced_view() # default load-balanced view

# get list of frequencies
freqs = json.loads(config.get("run", "frequency"))

if config.getboolean("run", "run_halfrings"):
    print "HALFRINGS"
    halfrings_survs = ["full"]
    for freq in freqs:
        chtags = [""]
        pol = "IQU"
        if freq == 70:
            chtags += ["18_23", "19_22", "20_21"]
        if freq >= 545:
            pol = "I"
        for chtag in chtags:
            for surv in halfrings_survs:
                if paral:
                    tasks.append(lview.apply_async(halfrings,freq, chtag,
                                                   surv, pol=pol,
                                                   smooth_combine_config=smooth_combine_config,
                                                   root_folder=root_folder,log_to_file=True,
                                                   mapreader=mapreader))
                else:
                    try:
                        halfrings(freq, chtag, surv, pol=pol,
                                  smooth_combine_config=smooth_combine_config,
                                  root_folder=root_folder,log_to_file=False,
                                 mapreader=mapreader)
                    except (NoOptionError, exceptions.IOError) as e:
                        log.error("SKIP TEST: " + e.message)

if config.getboolean("run", "run_surveydiff"):
    print "SURVDIFF"
    for bp_corr in [False]:
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
                if (freq >= 545) or (chtag and chtag.find("_") < 0):
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
                    try:
                        surveydiff(freq, chtag, survs, pol=pol,
                                   smooth_combine_config=smooth_combine_config,
                                   root_folder=root_folder,log_to_file=False,
                                   bp_corr=bp_corr, mapreader=mapreader)
                    except (NoOptionError, exceptions.IOError) as e:
                        log.error("SKIP TEST: " + e.message)

if config.getboolean("run", "run_chdiff"):
    print "CHDIFF"
        
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
                try:
                    chdiff(freq, ["LFI%d" % h for h in utils.HORNS[freq]], surv,
                          pol='I', smooth_combine_config=smooth_combine_config,
                          root_folder=root_folder,
                          log_to_file=False,
                          mapreader=mapreader)
                except (NoOptionError, exceptions.IOError) as e:
                    log.error("SKIP TEST: " + e.message)

if paral:
    print("Wait for %d tasks to complete" % len(tasks))
    tc.wait(tasks)
