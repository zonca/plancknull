from glob import glob
from ConfigParser import SafeConfigParser
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


def get_filename(filename_pattern):
    for pattern in [filename_pattern, filename_pattern.replace("_full","")]:
        filename = glob(pattern)
        if len(filename) == 1:
            filename = filename[0]
            log.debug("File: " + filename)
            return filename
    if len(filename) == 0:
        error_log = "No match for pattern " + filename_pattern
    else:
        error_log = "Multiple matches for pattern " + filename_pattern
    log.fatal(error_log)
    raise exceptions.IOError(error_log)

def arcmin2rad(arcmin):
    return np.radians(arcmin/60.)

FWHM = { 30:arcmin2rad(33.15873215), 44:arcmin2rad(28.08523439), 70:arcmin2rad(13.08124258)}


def type_of_channel_set(ch):
    """Returns a string that identifies the set of channels"""
    if ch == "":
        return "frequency"
    elif ch.find('_') >= 0:
        return "detset"
    elif len(ch) == 5: # this works for lfi only
        return "horn"
    else:
        return "channel"


class BaseMapReader:
    """Abstract class, all readers should provide this
    interface"""

    def __call__(self, freq, surv, chtag='', nside=None, halfring=0, pol="I"):
        """See docstrings of the child classes"""
        return np.zeros(hp.nside2npix(1024))


class DXReader(BaseMapReader):
    """All maps in a single folder, DX9 naming convention"""


    def __init__(self, config_filename, nside=None, debug=False):
        """
        nside : None or int
            if None matches any nside, otherwise integer nside
        """
        self.config = SafeConfigParser(); self.config.read(config_filename)
        self.nside = nside
        self.debug = debug

    def read_masks(self, freq):
        result = []
        filenames = [get_filename(self.config.get("Templates", mask_type).format(frequency=freq)) for mask_type in ["ps_mask", "spectra_mask"]]

        for file_name in filenames:
            result.append(np.logical_not(np.floor(hp.ud_grade(hp.read_map(file_name), self.nside)).astype(np.bool)))
        return tuple(result)

    def __call__(self, freq, surv, chtag='', halfring=0, pol="I", bp_corr=False):
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

        Returns
        -------
        maps : array or tuple of arrays
            single map or tuple of maps as returned by healpy.read_map
        """

        # type of channel
        channel_type = type_of_channel_set(chtag)

        # type of map
        is_survey = isinstance(surv, int)
        is_halfring = halfring != 0

        # stokes component
        stokes = stokes_I if channel_type in ["channel", "horn"] else stokes_IQU
        if isinstance(pol, str):
            components = [stokes.index(p) for p in pol]
            if len(components) == 1:
                components = components[0]
        else:
            components = pol

        if channel_type in ["channel", "horn"]:
            if freq > 70:
                freq = chtag
                chtag = ""
            else:
                chtag = chtag.translate(None, "LFI") # remove LFI from channel name

        # read_map
        output_map = []
        file_template_list = ["map", channel_type]
        if is_survey:
            file_template_list.append("survey")
        if is_halfring:
            file_template_list.append("halfring")
        file_template = "_".join(file_template_list)

        if channel_type == "horn":
            # horn maps created summing channel maps
            file_template = file_template.replace("horn", "channel")
            if freq > 70:
                tags = [chtag+'a', chtag+'b']
            else:
                tags = [chtag+'M', chtag+'S']
        else:
            tags = [chtag]

        file_parameters = {"frequency":freq, "survey":surv}
        if is_halfring:
            file_parameters["halfring"] = halfring

        for tag in tags:
            filename_pattern = self.config.get("Templates", file_template).format(channel=tag, **file_parameters)
            filename = get_filename(filename_pattern)
            log.info("components %s" % (str(components)))
            if not self.debug:
                output_map.append(hp.ma(hp.read_map(filename, components)))
            else:
                output_map.append(
                    np.zeros((np.size(components), hp.nside2npix(1024)))
                                 )
        if bp_corr:
            bp_corr_file_template = "map_iqucorrection"
            if is_survey:
                bp_corr_file_template += "_survey"
            bp_corr_filename_pattern = self.config.get("Templates", bp_corr_file_template).format(frequency=freq, survey=surv)
            bp_corr_filename = get_filename(bp_corr_filename_pattern)
            if not self.debug:
                corr_map = hp.ma(hp.read_map(bp_corr_filename, (0,1,2)))
            else:
                corr_map = [np.zeros(hp.nside2npix(1024)) for c in [0,1,2]]
            for comp, corr in zip(output_map[0], corr_map):
                comp += corr

        if channel_type == "horn":
            log.info("Combining maps in horn map")
            out = .5 * (output_map[0] + output_map[1])
        else:
            out = output_map[0]

        if self.nside:
            log.info("Downgrading to nside %d" % self.nside)
            power = None
            if pol in "ADF":
                log.info("Downgrading a covariance matrix")
                power = 2
            out = hp.ud_grade(out, self.nside,power=power)
        return out
