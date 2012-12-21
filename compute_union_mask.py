import utils
import reader
import sys

import healpy as hp
mapreader = reader.DXReader(sys.argv[1])
from planck.Planck import Planck
pl = Planck()

print "UNIONMASK"
nside = 2048
survs = [1,2,3,4,5]
freqs = pl.inst["HFI"].f.keys()

for freq in freqs:
    mask = utils.read_mask("/global/homes/z/zonca/p/masks/mask_4.fits", nside)
    #mask |= utils.read_mask("/global/project/projectdirs/planck/user/zonca/masks/wmap_polarization_analysis_mask_r9_7yr_v4.fits", nside=nside)
    chtags = [""]
    if freq == 70:
        chtags += ["18_23", "19_22", "20_21"]
    for chtag in chtags:
        for surv in survs:
            mask |= mapreader(freq, surv, chtag, halfring=0, pol='Q').mask
    hp.write_map("union_mask_%d.fits" % freq, mask)
