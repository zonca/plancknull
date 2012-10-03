import os
import json
import numpy as np
import healpy as hp

import sys
sys.path.append("../../")
from plancknull.differences import smooth_combine
from plancknull import reader

def test_smoothcombine():

    nside = 256
    freq = 30
    m = hp.ma(np.random.standard_normal(hp.nside2npix(nside)))
    m.mask = np.zeros(len(m), dtype=np.bool)

    Reader = reader.Readers[os.environ["NULLTESTS_ENV"]]

    mapreader = Reader(os.environ["DDX9_LFI"], nside=nside,
                       baseline_length=None)

    ps_mask, gal_mask = mapreader.read_masks(freq)
    #ps_mask = False
    #gal_mask = False
    smooth_combine([(m, 1)], [(np.ones_like(m), 1)], spectra=True, smooth_mask=ps_mask, spectra_mask=gal_mask, metadata={"file_type":""})


    # check chi-square
    metadata = json.load(open("out_map.json", "r"))
    assert np.abs(metadata["map_unsm_chi2"] - 1) < .01
    assert np.abs(metadata["map_chi2"] - 1) < .1
    # check spectrum
    cl = hp.read_cl("out_cl.fits")
    realization_wn = cl[200:].mean()
    assert np.abs(realization_wn - metadata["whitenoise_cl"]) < 1e-5

