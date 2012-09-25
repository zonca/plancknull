from glob import glob
import exceptions
import logging as log
import os.path
import numpy as np
import healpy as hp
import re

stokes = "IQUHABCDEF" # H for hits,
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


class SingleFolderDXReader(BaseMapReader):
    """All maps in a single folder, DX9 naming convention"""

    def __init__(self, folder, nside=None):
        self.folder = folder
        self.nside = nside

    def read_masks(self, freq):
        result = []
        file_names = [
            glob(os.path.join(self.folder, "MASKs",
                              'mask_ps_%dGHz_*.fits' % freq))[-1],
            glob(os.path.join(self.folder, "MASKs",
                              'union_mask_%d.fits' % freq))[-1]]

        for file_name in file_names:
            result.append(np.logical_not(np.floor(hp.ud_grade(hp.read_map(file_name), 1024)).astype(np.bool)))

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
        nside : None or int
            if None matches any nside, otherwise integer nside
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
        if self.nside:
            nside_tag = "%d" % self.nside
        else:
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
                comp += corr

        if is_horn:
            log.info("Combining maps in horn map")
            out = .5 * (output_map[0] + output_map[1])
        else:
            out = output_map[0]

        if self.nside:
            out = hp.ud_grade(out, self.nside)
        return out

################################################################################


class DPCDX9Reader(BaseMapReader):
    """All maps in a single folder, DX9 naming convention"""

    def __init__(self, folder, default_nside = 1024, debug_mode = False):
        self.folder = folder
        self.default_nside = default_nside
        self.debug_mode = debug_mode

    def read_map(self, path, components):
        if not self.debug_mode:
            return hp.ma(hp.read_map(path, components))
        else:
            log.info("Reading file '{0}'...".format(path))
            z = hp.ma(np.zeros(self.default_nside * self.default_nside * 12))
            if type(components) is int:
                return z
            else:
                return [z] * len(components)

    def read_masks(self, freq):
        result = []
        file_names = [
            glob(os.path.join(self.folder, "MASKs",
                              'mask_ps_{0}GHz_*.fits'.format(freq)))[-1],
            '/planck/sci_ops1/null_test_area/destriping_mask_{0}.fits'.format(freq)]

        for file_name in file_names:
            log.info("Reading file '%s'", file_name)
            result.append(np.logical_not(np.floor(hp.ud_grade(hp.read_map(file_name),
                                                              self.default_nside))
                                         .astype(np.bool)))

        return tuple(result)

    def __call__(self, freq, surv, chtag='', nside=None, halfring=0, pol="I",
                 bp_corr = False):
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
                       'halfring': '',
                       'map_type': ''}

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

        list_of_filenames = []
        if chtag == "":
            format_dict['map_type'] = 'frequency map'
            # We look for a frequency map
            if halfring > 0:
                base_path = os.path.join(base_path, "JackKnife_DX9")
            elif type(surv) is int:
                base_path = os.path.join(base_path, "Surveys_DX9")

            list_of_filenames = glob(os.path.join(base_path,
                                          "LFI_{freq}_{nside}_????????_{halfring}{survey}.fits"
                                          .format(**format_dict)))

            if len(list_of_filenames) > 1:
                # Take the last one, as it is likely to be the most recent
                list_of_filenames = [list_of_filenames.sorted()[-1]]

        elif radiometer_match:
            format_dict['map_type'] = 'single radiometer map'
            # We look for a radiometer map
            list_of_filenames = glob(os.path.join(base_path, "SINGLE_horn_Survey",
                                          "LFI_{freq}_{nside}_????????_{rad}_{halfring}{survey}.fits"
                                          .format(rad=radiometer_match.group(1),
                                                  **format_dict)))

            if len(list_of_filenames) > 1:
                # Take the last one, as it is likely to be the most recent
                list_of_filenames = [list_of_filenames.sorted()[-1]]

        elif horn_match:
            format_dict['map_type'] = 'single horn map'
            # We look for a horn map
            output_map = None
            for rad in [horn_match.group(1) + arm for arm in ('M', 'S')]:
                mask = ("LFI_{freq}_{nside}_????????_{rad}_{halfring}{survey}.fits"
                        .format(rad=rad, **format_dict))
                match = glob(os.path.join(base_path, "SINGLE_horn_Survey", mask))
                list_of_filenames.append(sorted(match)[-1])

        elif quadruplet_match:
            format_dict['map_type'] = 'horn pair map'
            # We look for a quadruplet map
            if halfring > 0:
                base_path = os.path.join(base_path, "JackKnife_DX9")
            elif surv in ("nominal", "full"):
                base_path = os.path.join(base_path, "Couple_horn_DX9")
            elif type(surv) is int:
                base_path = os.path.join(base_path, "Couple_horn_Surveys_DX9")

            list_of_filenames = glob(os.path.join(base_path,
                                          "LFI_{freq}_{nside}_????????_{quadruplet}_{halfring}{survey}.fits"
                                          .format(quadruplet=quadruplet_match.group(1),
                                                  **format_dict)))

            if len(list_of_filenames) > 1:
                # Take the last one, as it is likely to be the most recent
                list_of_filenames = [list_of_filenames.sorted()[-1]]

        if not list_of_filenames:
            raise RuntimeError(("Unable to find a match for {map_type} "
                                "(freq: '{freq}', "
                                "nside: '{nside}', survey: '{survey}', "
                                "halfring: '{halfring}')")
                               .format(**format_dict))

        components = [stokes.index(p) for p in pol]
        if len(components) == 1:
            components = components[0]

        output_map = None
        for cur_filename in list_of_filenames:
            log.info("Reading file {0}".format(os.path.basename(cur_filename)))

            cur_map = self.read_map(cur_filename, components)
            if output_map is None:
                output_map = cur_map
            else:
                if np.ndim(cur_map) == 1:
                    # This is just a temperature map
                    output_map = output_map + cur_map
                else:
                    # This is a IQU map, so we have to iterate over
                    # the three components
                    for idx in range(len(cur_map)):
                        output_map[idx] = output_map[idx] + cur_map[idx]

        if len(list_of_filenames) > 1:
            # If we loaded more than one file, average all the maps
            # into one. We assume here that in this case all the maps
            # were temperature-only. (In fact, we load more than one
            # file only when building a horn map from two radiometer
            # maps.)
            output_map = (1.0 / len(list_of_filenames)) \
                * np.sum(np.array(output_map), axis=0)

        # Apply the bandpass correction
        if bp_corr:
            bp_corr_filename = "iqu_bandpass_correction_%d_" % freq
            if surv in ["nominal", "full"]:
                bp_corr_filename += format_dict["survey"] + "survey"
            else:
                bp_corr_filename += format_dict["survey"].replace("survey_", "ss")
            bp_corr_filename += ".fits"
            log.info("Applying bandpass correction: " + bp_corr_filename)
            corr_map = self.read_map(os.path.join(self.folder,
                                                  "IQU_Corrections_Maps",
                                                  bp_corr_filename),
                                     (0, 1, 2))
            for comp, corr in zip(output_map, corr_map):
                comp += corr

        return output_map

################################################################################


class SingleFolderToastReader(BaseMapReader):
    """All maps in a single folder, Toast naming convention"""

    def __init__(self, folder, nside=None):
        self.folder = folder
        self.nside = nside

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

        if chtag.find('_')<0 and len(chtag) == 5: #horn
            is_horn = True
        else:
            is_horn = False

        # subset tag
        subset_halfring_tag = chtag
        if chtag:
            if chtag.find('_') > 0: #is quadruplet
                subset_halfring_tag = "_".join( ["LFI%sM_LFI%sS" % (horn,horn) for horn in chtag.split("_")])
        else:
            subset_halfring_tag = "%03d" % freq

        # halfrings
        halfring_tag = ""
        if halfring != 0:
            halfring_tag = "_subchunk_%d_of_2" % halfring

        # survey
        if isinstance(surv, int):
            surv = "survey%d" % surv

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
            filename_pattern = os.path.join(folder, "map_ddx9_%s_%s%s.fits" % (tag, surv, halfring_tag))
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
            #output_map.append(None)

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
                comp += corr

        if is_horn:
            log.info("Combining maps in horn map")
            out = .5 * (output_map[0] + output_map[1])
        else:
            out = output_map[0]

        if self.nside:
            out = hp.ud_grade(out, self.nside)
        return out

################################################################################

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
