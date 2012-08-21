import numpy as np
import logging as log
import healpy as hp
from glob import glob
import os
from differences import halfrings, surveydiff, chdiff 
from reader import SingleFolderDXReader

def read_dpc_masks(freq, nside):
    ps_mask = np.logical_not(np.floor(hp.ud_grade( 
    hp.read_map(
        glob(os.path.join(os.environ["DX9_LFI"], "MASKs",'mask_ps_%dGHz_*.fits' % freq))[0]), nside))
    )
    gal_mask = np.logical_not(hp.read_map(
        glob(os.path.join(os.environ["DX9_LFI"], "MASKs",'destriping_mask_%d.fits' % freq))[0])
        )
    return ps_mask, gal_mask

NSIDE = 1024

log.root.level = log.DEBUG
mapreader = SingleFolderDXReader(os.environ["DX9_LFI"])

survs = ["nominal", "full"]
for freq in [30, 44, 70]:
    ps_mask, gal_mask = read_dpc_masks(freq, NSIDE)
    smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128,smooth_mask=ps_mask, spectra_mask=gal_mask)
    chtags = [""]
    if freq == 70:
        chtags += ["18_23", "19_22", "20_21"]
    for chtag in chtags:
        for surv in survs:
            halfrings(freq, chtag, surv, pol='IQU', smooth_combine_config=smooth_combine_config, mapreader=mapreader, output_folder="dx9/halfrings/")
