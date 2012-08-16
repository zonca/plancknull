from glob import glob
import exceptions
import logging as log
import os.path
import numpy as np
import healpy as hp

from planck.Planck import Planck

PLANCK = Planck()
stokes = "IQU"

class BaseMapReader:
    """Abstract class, all readers should provide this
    interface"""

    def __call__(self, freq, surv, chtag='', halfring=0, pol='I'):
        """Read a map and return the array of pixels.
        
        Parameters
        ----------
        freq : int
            frequency
        chtag : string 
            can be "" for frequency, radiometer("LFI18S"), horn("LFI18"), quadruplet("LFI18-LFI23"), detset("detset_1")
        surv : int or string
            "nominal", "full", or survey number

        Returns
        -------
        maps : array or tuple of arrays
            maps as returned by healpy.read_map

        """
        return np.zeros(hp.nside2npix(1024))

class SingleFolderDXReader:
    """All maps in a single folder, DX9 naming convention"""

    def __init__(self, folder):
        self.folder = folder

    def __call__(self, freq, surv, chtag='', halfring=0, pol="I"):
        # stokes component
        components = [stokes.index(p) for p in pol]
        if len(components) == 1:
            components = components[0]

        # subset tag
        subset_halfring_tag = chtag
        if subset_halfring_tag:
            subset_halfring_tag = "_" + subset_halfring_tag

        # halfrings
        if halfring != 0:
            subset_halfring_tag += "_ringhalf_%d" % halfring

        # survey
        if isinstance(surv, int):
            surv = "survey_%d" % surv

        # read_map
        filename = glob(os.path.join(self.folder, "???_%d_????_????????%s_%s.fits*" % (freq, subset_halfring_tag, surv)))
        log.info("Reading %s" % filename)
        return hp.read_map(filename, components)
