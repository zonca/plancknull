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

    def __call__(self, freq, surv, nside=None, chtag='', halfring=0, pol='I'):
        """Read a map and return the array of pixels.
        
        Parameters
        ----------
        freq : int
            frequency
        surv : int or string
            "nominal", "full", or survey number
        nside : None or int
            if None matches any nside, otherwise integer nside
        chtag : string 
            can be "" for frequency, radiometer("LFI18S"), horn("LFI18"), quadruplet("18_23"), detset("detset_1")

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

    def __call__(self, freq, surv, nside=None, chtag='', halfring=0, pol="I"):
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

        # nside
        if nside:
            nside_tag = "%d" % nside
        else:
            nside_tag = "????"

        # read_map
        filename_pattern = os.path.join(self.folder, "???_%d_%s_????????%s_%s.fits*" % (freq, nside_tag, subset_halfring_tag, surv))
        filename = glob(filename_pattern)
        if len(filename) == 1:
            filename = filename[0]
        else: 
            if len(filename) == 0:
                error_log = "No match for pattern " + filename_pattern
            else:
                error_log = "Multiple matches for pattern " + filename_pattern
            log.fatal(error_log)
            raise exceptions.IOError(error_log)

        log.info("Reading %s" % os.path.basename(filename))
        #return hp.read_map(filename, components)
        return None
