import numpy as np
import logging as log
import healpy as hp
from glob import glob
import os
from differences import halfrings, surveydiff, chdiff 
from reader import SingleFolderDXReader
from IPython.parallel import Client

os.environ["DX9_LFI"] = "/global/project/projectdirs/planck/data/mission/DPC_maps/dx9/lfi/"

def read_dpc_masks(freq):
    import numpy as np
    import os
    import healpy as hp
    from glob import glob
    if freq > 70:
        nside = 2048
    else:
        nside = 1024
    ps_mask = np.logical_not(np.floor(hp.ud_grade( 
    hp.read_map(
        glob(os.path.join(os.environ["DX9_LFI"], "MASKs",'mask_ps_%dGHz_*.fits' % freq))[0]), nside))
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
mapreader = None
log.root.level = log.DEBUG

if __name__ == '__main__':

    tasks = []
    tc = Client()
    lview = tc.load_balanced_view() # default load-balanced view

    #print "HALFRINGS"
    #survs = ["nominal", "full"]
    #freqs = [30,44,70]
    #for freq in freqs:
    #    smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
    #    chtags = [""]
    #    if freq == 70:
    #        chtags += ["18_23", "19_22", "20_21"]
    #    for chtag in chtags:
    #        for surv in survs:
    #            tasks.append(lview.apply_async(halfrings,freq, chtag, surv, pol='IQU', smooth_combine_config=smooth_combine_config, mapreader=mapreader, output_folder="dx9/halfringspar/",read_masks=read_dpc_masks))

    #print "SURVDIFF"
    #survs = [1,2,3,4,5]
    #freqs = [30, 44, 70]
    #for freq in freqs:
    #    smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
    #    chtags = [""]
    #    if freq == 70:
    #        chtags += ["18_23", "19_22", "20_21"]
    #    chtags += chlist(freq)
    #    for chtag in chtags:
    #        if chtag and chtag.find("_") < 0:
    #            pol='I'
    #        else:
    #            pol="IQU"
    #        tasks.append(lview.apply_async(surveydiff,freq, chtag, survs, pol=pol, smooth_combine_config=smooth_combine_config, mapreader=mapreader, output_folder="dx9/surveydiffpar/",read_masks=read_dpc_masks))

    print "CHDIFF"
    survs = [1,2,3,4,5]
    freqs = [30, 44, 70]
        
    for freq in freqs:
        smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128)
        for surv in survs:
           tasks.append(lview.apply_async(chdiff,freq, ["LFI%d" % h for h in HORNS[freq]], surv, pol='I', smooth_combine_config=smooth_combine_config, mapreader=mapreader, output_folder="dx9/chdiff/", read_masks=read_dpc_masks,log_to_file=True))

    print("Wait for %d tasks to complete" % len(tasks))
    tc.wait(tasks)
