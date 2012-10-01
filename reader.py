from glob import glob
import exceptions
import logging as log
import os.path
import numpy as np
import healpy as hp

stokes_IQU = "IQUHABCDEF"
stokes_I = "IHA" 
# H for hits,
# ABCDEF 6 components of VARIANCE matrix
# II, IQ, IU, QQ, QU, UU

def arcmin2rad(arcmin):
    return np.radians(arcmin/60.)

FWHM = { 30:arcmin2rad(33.15873215), 44:arcmin2rad(28.08523439), 70:arcmin2rad(13.08124258)}


def type_of_channel_set(ch):
    """Returns a string that identifies the set of channels"""
    if ch == "":
        return "frequency"
    elif ch.find('_') >= 0:
        return "detset"
    else:
        return "single_ch"


class BaseMapReader:
    """Abstract class, all readers should provide this
    interface"""

    def __call__(self, freq, surv, chtag='', nside=None, halfring=0, pol="I"):
        """See docstrings of the child classes"""
        return np.zeros(hp.nside2npix(1024))


class DXReader(BaseMapReader):
    """All maps in a single folder, DX9 naming convention"""


    def __init__(self, folder, nside=None, debug=False):
        """
        nside : None or int
            if None matches any nside, otherwise integer nside
        """
        self.folder = folder
        self.nside = nside
        self.subfolder = { "ps_masks":os.path.join(self.folder, "MASKs") }
        self.subfolder["spectra_masks"] = self.subfolder["ps_masks"]
        self.subfolder["bandpass_corrections"] = os.path.join(self.folder, "bandpass_correction")
        self.subfolder["single_channels"] = os.path.join(self.folder, "channels")
        self.debug = debug

    def read_masks(self, freq):
        result = []
        file_names = [
            glob(os.path.join(self.subfolder["ps_masks"],
                              'mask_ps_%dGHz_*.fits' % freq))[-1],
            glob(os.path.join(self.subfolder["spectra_masks"],
                              'union_mask_%d.fits' % freq))[-1]]

        for file_name in file_names:
            result.append(np.logical_not(np.floor(hp.ud_grade(hp.read_map(file_name), self.nside)).astype(np.bool)))
        return tuple(result)

    def __call__(self, freq, surv, chtag='', halfring=0, pol="I", bp_corr=False, baseline_length="1s"):
        """Read a map and return the array of pixels.

        Parameters
        ----------
        freq : int
            frequency
        surv : int or string
            "nominal", "full", or survey number
        chtag : string
            can be "" for frequency, radiometer("LFI18S"), horn("LFI18"), quadruplet("18_23"), detset("detset_1")
        halfring : int
            0 for full, 1 and 2 for first and second halfrings
        pol : string
            required polarization components, e.g. 'I', 'Q', 'IQU'
        baseline_length : string
            baseline length string, typically "1m" or "1s"
            if None, it is not used to match the filename

        Returns
        -------
        maps : array or tuple of arrays
            single map or tuple of maps as returned by healpy.read_map
        """

        # single channel
        is_single_channel = chtag and chtag.find('_') < 0 # single channels do not have underscores
        is_horn = chtag.find('_')<0 and len(chtag) == 5
        is_subset = chtag and chtag.find('_')>0 #quadruplets for LFI, subsets for HFI

        # stokes component
        stokes = stokes_I if is_single_channel else stokes_IQU
        if isinstance(pol, str):
            components = [stokes.index(p) for p in pol]
            if len(components) == 1:
                components = components[0]
        else:
            components = pol

        # folder
        folder = self.folder


        if is_single_channel or is_horn:
            if freq > 70:
                freq = chtag
                chtag = ""
            else:
                chtag = chtag.translate(None, "LFI") # remove LFI from channel name

        # subset tag
        subset_halfring_tag = chtag
        if subset_halfring_tag:
            subset_halfring_tag = "_" + subset_halfring_tag


        if is_single_channel:
            folder = self.subfolder.get("single_channels", folder)
        elif is_subset and isinstance(surv, int):
            folder = self.subfolder.get("subsets_surveys", folder)
        elif is_subset:
            folder = self.subfolder.get("subsets", folder)
        elif isinstance(surv, int):
            folder = self.subfolder.get("surveys", folder)

        # halfrings
        if halfring != 0:
            folder = self.subfolder.get("halfrings", folder) #overwrites is_subset, surv folder
            subset_halfring_tag += "_ringhalf_%d" % halfring

        # survey
        if isinstance(surv, int):
            surv = "survey_%d" % surv

        # nside
        nside_tag = "????"

        # baseline length
        baseline_length_tag = ''
        if baseline_length:
            baseline_length_tag = '_' + baseline_length

        # read_map
        output_map = []
        if is_horn:
            if freq > 70:
                tags = [subset_halfring_tag.replace(chtag, chtag+'a'), subset_halfring_tag.replace(chtag, chtag+'b')]
            else:
                tags = [subset_halfring_tag.replace(chtag, chtag+'M'), subset_halfring_tag.replace(chtag, chtag+'S')]
        else:
            tags = [subset_halfring_tag]
        for tag in tags:
            filename_pattern = os.path.join(folder, "???_%s_%s_????????%s_%s%s.fits*" % (str(freq), nside_tag, tag, surv, baseline_length_tag))
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

            log.info("Reading %s, components %s" % (os.path.basename(filename), str(components)))
            if not self.debug:
                output_map.append(hp.ma(hp.read_map(filename, components)))
            else:
                output_map.append(
                    np.zeros((np.size(components), hp.nside2npix(1024)))
                                 )

        if bp_corr:
            bp_corr_filename = "iqu_bandpass_correction_%d_" % freq
            if surv in ["nominal", "full"]:
                bp_corr_filename += surv + "survey"
            else:
                bp_corr_filename += surv.replace("survey_", "ss")
            bp_corr_filename += ".fits"
            log.info("Applying bandpass correction: " + bp_corr_filename)
            if not self.debug:
                corr_map = hp.ma(hp.read_map(os.path.join(self.subfolder["bandpass_corrections"], bp_corr_filename), (0,1,2)))
            else:
                corr_map = [np.zeros(hp.nside2npix(1024)) for c in [0,1,2]]
            for comp, corr in zip(output_map[0], corr_map):
                comp += corr

        if is_horn:
            log.info("Combining maps in horn map")
            out = .5 * (output_map[0] + output_map[1])
        else:
            out = output_map[0]

        if self.nside:
            out = hp.ud_grade(out, self.nside)
        return out

class DPCDXReader(DXReader):

    def __init__(self, folder, nside=None, debug=False):
        """
        nside : None or int
            if None matches any nside, otherwise integer nside
        """
        self.folder = folder
        self.nside = nside
        self.subfolder = { "ps_masks":os.path.join(self.folder, "MASKs") }
        self.subfolder["spectra_masks"] = "/planck/sci_ops1/null_test_area/"
        self.subfolder["halfrings"] = os.path.join(self.folder, "JackKnife")
        self.subfolder["surveys"] = os.path.join(self.folder, "Surveys")
        self.subfolder["subsets"] = os.path.join(self.folder, "Couple_horn")
        self.subfolder["subsets_surveys"] = os.path.join(self.folder, "Couple_horn_Surveys")
        self.subfolder["bandpass_corrections"] = os.path.join(self.folder, "IQU_Corrections_Maps")
        self.subfolder["single_channels"] = os.path.join(self.folder, "Single_Radiometer")
        self.debug = debug

Readers = {"LFIDPC":DPCDXReader, "NERSC":DXReader}

import unittest

class TestReader(unittest.TestCase):
    def setUp(self):
        self.reader = DPCDX9Reader(os.environ["DX9_LFI"])

    def test_reader(self):
        for frequency in (30, 44):
            for survey in ("nominal", "full", 1):
                cur_map = self.reader(frequency, survey, "", pol = "IQU")
                self.assertEqual(np.ndim(cur_map), 2)

        for survey in (1, 2):
            cur_map = self.reader(70, survey, "LFI18S")
            self.assertEqual(np.ndim(cur_map), 1)

            cur_map = self.reader(70, survey, "LFI18")
            self.assertEqual(np.ndim(cur_map), 1)

            cur_map = self.reader(70, survey, "18_23", pol = "IQU")
            self.assertEqual(np.ndim(cur_map), 2)

        cur_map = self.reader(30, "nominal", "", halfring=1)
        self.assertEqual(np.ndim(cur_map), 2)

if __name__ == "__main__":
    unittest.main()
