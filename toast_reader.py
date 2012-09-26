from reader import *

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

