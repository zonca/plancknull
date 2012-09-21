import numpy as np
import logging as log
import os
from differences import halfrings, surveydiff, chdiff 

paral = True
if paral:
    from IPython.parallel import Client

from reader import SingleFolderDXReader

NSIDE = 1024

HORNS = {30:[27,28], 44:[24,25,26], 70:list(range(18,23+1))}

def chlist(freq):
    horns = HORNS[freq]
    chs = []
    for horn in horns:
        chs += ["LFI%dM" % horn, "LFI%dS" % horn]
    return chs

log.root.level = log.DEBUG
mapreader = SingleFolderDXReader(os.environ["DDX9_LFI"])

if __name__ == '__main__':

    if paral:
        tasks = []
        tc = Client()
        lview = tc.load_balanced_view() # default load-balanced view

    root_folder = "ddx9_1deg"
    run_halfrings = True
    run_surveydiff = True
    run_chdiff = False

    if run_halfrings:
        print "HALFRINGS"
        survs = ["nominal", "full"]
        freqs = [30,44,70]
        for freq in freqs:
            smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
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

    if run_surveydiff:
        print "SURVDIFF"
        survs = [1,2,3,4,5]
        freqs = [30, 44, 70]
        for bp_corr in [False, True]:
            for freq in freqs:
                smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
                chtags = [""]
                if freq == 70:
                    chtags += ["18_23", "19_22", "20_21"]
                #chtags += chlist(freq)
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

    if run_chdiff:
        print "CHDIFF"
        survs = [1,2,3,4,5]
        freqs = [30, 44, 70]
            
        for freq in freqs:
            smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
            for surv in survs:
                if paral:
                   tasks.append(lview.apply_async(chdiff,freq, ["LFI%d" % h for
                                                                h in
                                                                HORNS[freq]],
                                                  surv, pol='I',
                                                  smooth_combine_config=smooth_combine_config,
                                                  root_folder=root_folder,
                                                  log_to_file=True,
                                                  mapreader=mapreader))
                else:
                   chdiff(freq, ["LFI%d" % h for h in HORNS[freq]], surv,
                          pol='I', smooth_combine_config=smooth_combine_config,
                          root_folder=root_folder,
                          log_to_file=True,
                          mapreader=mapreader)

    if paral:
        print("Wait for %d tasks to complete" % len(tasks))
        tc.wait(tasks)
