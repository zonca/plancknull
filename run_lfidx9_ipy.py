import numpy as np
import logging as log
import healpy as hp
from glob import glob
import os
from differences import halfrings, surveydiff, chdiff 
from IPython.parallel import Client

def read_dpc_masks(freq):
    import numpy as np
    import os
    import healpy as hp
    from glob import glob
    ps_mask = np.logical_not(np.floor(hp.ud_grade( 
    hp.read_map(
        glob(os.path.join(os.environ["DX9_LFI"], "MASKs",'mask_ps_%dGHz_*.fits' % freq))[0]), 1024))
    ).astype(np.bool)
    gal_filename = glob(os.path.join(
        os.environ["DX9_LFI"], "MASKs",
        'destripingmask_%d.fits' % freq))[0]
    gal_mask = np.logical_not(hp.read_map(gal_filename)).astype(np.bool)
    return ps_mask, gal_mask

HORNS = {30:[27,28], 44:[24,25,26], 70:list(range(18,23+1))}

def chlist(freq):
    horns = HORNS[freq]
    chs = []
    for horn in horns:
        chs += ["LFI%dM" % horn, "LFI%dS" % horn]
    return chs

NSIDE = 1024
log.root.level = log.DEBUG

paral = True
if __name__ == '__main__':

    if paral:
        tasks = []
        tc = Client()
        lview = tc.load_balanced_view() # default load-balanced view

    root_folder = "dx9"
    run_halfrings = True
    run_surveydiff = True
    run_chdiff = True

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
                        tasks.append(lview.apply_sync(halfrings,freq, chtag, surv, pol='IQU', smooth_combine_config=smooth_combine_config, root_folder=root_folder,read_masks=read_dpc_masks,log_to_file=True))
                    else:
                        halfrings(freq, chtag, surv, pol='IQU', smooth_combine_config=smooth_combine_config, root_folder=root_folder,read_masks=read_dpc_masks,log_to_file=False)

    if run_surveydiff:
        print "SURVDIFF"
        survs = [1,2,3,4,5]
        freqs = [30, 44, 70]
        for freq in freqs:
            smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
            chtags = [""]
            if freq == 70:
                chtags += ["18_23", "19_22", "20_21"]
            chtags += chlist(freq)
            for chtag in chtags:
                if chtag and chtag.find("_") < 0:
                    pol='I'
                else:
                    pol="IQU"
                if paral:
                    tasks.append(lview.apply_sync(surveydiff,freq, chtag, survs, pol=pol, smooth_combine_config=smooth_combine_config, root_folder=root_folder,read_masks=read_dpc_masks,log_to_file=True))
                else:
                    surveydiff(freq, chtag, survs, pol=pol, smooth_combine_config=smooth_combine_config, root_folder=root_folder,read_masks=read_dpc_masks,log_to_file=True)

    if run_chdiff:
        print "CHDIFF"
        survs = [1,2,3,4,5]
        freqs = [30, 44, 70]
            
        for freq in freqs:
            smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
            for surv in survs:
                if paral:
                   tasks.append(lview.apply_sync(chdiff,freq, ["LFI%d" % h for h in HORNS[freq]], surv, pol='I', smooth_combine_config=smooth_combine_config, root_folder=root_folder, read_masks=read_dpc_masks,log_to_file=True))
                else:
                   chdiff(freq, ["LFI%d" % h for h in HORNS[freq]], surv, pol='I', smooth_combine_config=smooth_combine_config, root_folder=root_folder, read_masks=read_dpc_masks,log_to_file=True)

    if paral:
        print("Wait for %d tasks to complete" % len(tasks))
        tc.wait(tasks)
