from glob import glob
import utils
import os
import reader
import sys

import healpy as hp
mapreader = reader.DXReader(sys.argv[1])

print "UNIONMASK"
nside = 1024
survs = [1,2,3,4,5]
freqs = [30, 44, 70]
for freq in freqs:
    mask = utils.read_mask(
        glob(os.path.join(os.environ["DDX9_LFI"], "MASKs",
                      'destripingmask_%d.fits' % freq))[-1]
                     , nside=nside )
    mask |= utils.read_mask("/global/project/projectdirs/planck/user/zonca/masks/wmap_polarization_analysis_mask_r9_7yr_v4.fits", nside=nside)
    chtags = [""]
    if freq == 70:
        chtags += ["18_23", "19_22", "20_21"]
    for chtag in chtags:
        for surv in survs:
            mask |= mapreader(freq, surv, chtag, halfring=0, pol='Q').mask
    hp.write_map("union_mask_%d.fits" % freq, mask)
