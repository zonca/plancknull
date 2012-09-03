from glob import glob
import exceptions
import logging as log
import os.path
import numpy as np
import healpy as hp
import re

stokes = "IQU"

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

class SingleFolderDXReader(BaseMapReader):
    """All maps in a single folder, DX9 naming convention"""

    def __init__(self, folder):
        self.folder = folder

    def __call__(self, freq, surv, chtag='', nside=None, halfring=0, pol="I", bp_corr=False):
        """Read a map and return the array of pixels.
        
        Parameters
        ----------
        freq : int
            frequency
        surv : int or string
            "nominal", "full", or survey number
        chtag : string 
            can be "" for frequency, radiometer("LFI18S"), horn("LFI18"), quadruplet("18_23"), detset("detset_1")
        nside : None or int
            if None matches any nside, otherwise integer nside
        halfring : int
            0 for full, 1 and 2 for first and second halfrings
        pol : string
            required polarization components, e.g. 'I', 'Q', 'IQU'

        Returns
        -------
        maps : array or tuple of arrays
            single map or tuple of maps as returned by healpy.read_map
        """
        # stokes component
        components = [stokes.index(p) for p in pol]
        if len(components) == 1:
            components = components[0]

        # folder
        folder = self.folder

        # single channel
        if chtag and chtag.find('_') < 0:
            # single channels do not have underscores
            if freq > 70:
                freq = chtag
                chtag = ""
            else:
                folder = os.path.join(self.folder, "channels")
                chtag = chtag.translate(None, "LFI") # remove LFI from channel name

        if chtag.find('_')<0 and len(chtag) == 2: #horn
            is_horn = True
        else:
            is_horn = False

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
        output_map = []
        if is_horn:
            if freq > 70:
                tags = [subset_halfring_tag.replace(chtag, chtag+'a'), subset_halfring_tag.replace(chtag, chtag+'b')]
            else:
                tags = [subset_halfring_tag.replace(chtag, chtag+'M'), subset_halfring_tag.replace(chtag, chtag+'S')]
        else:
            tags = [subset_halfring_tag]
        for tag in tags:
            filename_pattern = os.path.join(folder, "???_%s_%s_????????%s_%s.fits*" % (str(freq), nside_tag, tag, surv))
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
            output_map.append(hp.ma(hp.read_map(filename, components)))

        if bp_corr:
            bp_corr_filename = "iqu_bandpass_correction_%d_" % freq
            if surv in ["nominal", "full"]:
                bp_corr_filename += surv + "survey"
            else:
                bp_corr_filename += surv.replace("survey_", "ss")
            bp_corr_filename += ".fits"
            log.info("Applying bandpass correction: " + bp_corr_filename)
            corr_map = hp.ma(hp.read_map(os.path.join(folder, "bandpass_correction", bp_corr_filename), (0,1,2)))
            for comp, corr in zip(output_map[0], corr_map):
                comp -= corr

        if is_horn:
            log.info("Combining maps in horn map")
            return .5 * (output_map[0] + output_map[1])
        else:
            return output_map[0]

class DPCDX9Reader(BaseMapReader):
    """All maps in a single folder, DX9 naming convention"""

    def __init__(self, folder):
        self.folder = folder

    def __call__(self, freq, surv, chtag='', nside=None, halfring=0, pol="I"):
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
        halfring : int
            0 for full, 1 and 2 for first and second halfrings
        pol : string
            required polarization components, e.g. 'I', 'Q', 'IQU'

        Returns
        -------
        maps : array or tuple of arrays
            single map or tuple of maps as returned by healpy.read_map
        """

        radiometer_regex = re.compile("^LFI([0-9][0-9][MS])$")
        horn_regex = re.compile("^LFI([0-9][0-9])$")
        quadruplet_regex = re.compile("^([0-9][0-9]_[0-9][0-9])$")

        base_path = self.folder

        format_dict = {'freq': '',
                       'nside': '',
                       'survey': '',
                       'halfring': ''}

        if freq in (30, 44, 70):
            format_dict['freq'] = freq
        else:
            format_dict['freq'] = '*'

        if nside is None:
            format_dict['nside'] = '*'  # This will be filled by 'glob'
        else:
            format_dict['nside'] = nside

        if surv == "nominal":
            format_dict['survey'] = "nominal"
        elif surv == "full":
            format_dict['survey'] = "full"
        elif type(surv) is int:
            format_dict['survey'] = "survey_%d" % surv
        else:
            raise ValueError("Unknown tag '{0}' for 'surv'".format(surv))

        if halfring == 1:
            format_dict['halfring'] = "ringhalf_1_"
        elif halfring == 2:
            format_dict['halfring'] = "ringhalf_2_"

        radiometer_match = radiometer_regex.match(chtag)
        horn_match = horn_regex.match(chtag)
        quadruplet_match = quadruplet_regex.match(chtag)

        filenames = []
        if chtag == "":
            # We look for a frequency map
            if type(surv) is int:
                base_path = os.path.join(base_path, "Surveys_DX9")

            filenames = glob(os.path.join(base_path,
                                          "LFI_{freq}_{nside}_????????_{halfring}{survey}.fits"
                                          .format(**format_dict)))

            if len(filenames) > 1:
                # Take the last one, as it is likely to be the most recent
                filenames = [filenames.sorted()[-1]]

        elif radiometer_match:
            # We look for a radiometer map
            filenames = glob(os.path.join(base_path, "SINGLE_horn_Survey",
                                          "LFI_{freq}_{nside}_????????_{rad}_{halfring}{survey}.fits"
                                          .format(rad=radiometer_match.group(1),
                                                  **format_dict)))

            if len(filenames) > 1:
                # Take the last one, as it is likely to be the most recent
                filenames = [filenames.sorted()[-1]]

        elif horn_match:
            # We look for a horn map
            output_map = None
            for rad in [horn_match.group(1) + arm for arm in ('M', 'S')]:
                mask = ("LFI_{freq}_{nside}_????????_{rad}_{halfring}{survey}.fits"
                        .format(rad=rad, **format_dict))
                match = glob(os.path.join(base_path, "SINGLE_horn_Survey", mask))
                filenames.append(sorted(match)[-1])

        elif quadruplet_match:
            # We look for a quadruplet map
            if surv in ("nominal", "full"):
                base_path = os.path.join(base_path, "Couple_horn_DX9")
            else:
                base_path = os.path.join(base_path, "Couple_horn_Surveys_DX9")

            filenames = glob(os.path.join(base_path,
                                          "LFI_{freq}_{nside}_????????_{quadruplet}_{halfring}{survey}.fits"
                                          .format(quadruplet=quadruplet_match.group(1),
                                                  **format_dict)))

            if len(filenames) > 1:
                # Take the last one, as it is likely to be the most recent
                filenames = [filenames.sorted()[-1]]

        if not filenames:
            raise RuntimeError(("Unable to find a match (freq: '{freq}', "
                                "nside: '{nside}', survey: '{survey}', "
                                "halfring: '{halfring}')")
                               .format(**format_dict))

        components = [stokes.index(p) for p in pol]
        if len(components) == 1:
            components = components[0]

        output_map = None
        for cur_filename in filenames:
            log.info("Reading file {}".format(os.path.basename(cur_filename)))
            
            cur_map = hp.ma(hp.read_map(cur_filename, components))
            if output_map is None:
                output_map = cur_map
            else:
                if len(cur_map) == 1:
                    output_map = output_map + cur_map
                else:
                    for idx in range(len(cur_map)):
                        output_map[idx] = output_map[idx] + cur_map[idx]

        if len(output_map) == 1:
            return output_map * (1.0 / len(filename))
        else:
            return tuple([output_map[idx] * (1.0 / len(filename))
                          for idx in range(len(output_map))])
