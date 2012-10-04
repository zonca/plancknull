import os
import json
import exceptions
import itertools
import numpy as np
import logging as log
import healpy as hp
import reader

import utils

def configure_file_logger(base_filename):
    rl = log.root
    handler = log.FileHandler(base_filename + ".log", mode='w')
    handler.setLevel(log.DEBUG)
    handler.setFormatter(log.Formatter('%(asctime)s %(levelname)-8s %(message)s'))
    for h in rl.handlers:
        h.stream.close()
        rl.removeHandler(h)
    rl.level = log.DEBUG
    rl.addHandler( handler )
    log.info("Start logging")


def combine_maps(maps_and_weights):
    """Combine maps with given weights

    """
    is_IQU = len(maps_and_weights[0][0]) == 3
    # combine
    if is_IQU:
        combined_map = [maps_and_weights[0][0][comp] * maps_and_weights[0][1] for comp in range(3)]  
        for comp in range(3):
            for m,w in maps_and_weights[1:]:
                combined_map[comp] += m[comp] * w
    else: # single component only
        # combined_map is a list of 1 element
        combined_map = [maps_and_weights[0][0] * maps_and_weights[0][1]]
        for m,w in maps_and_weights[1:]:
            combined_map[0] += m * w

    return combined_map

def smooth_combine(maps_and_weights, variance_maps_and_weights, fwhm=np.radians(2.0), degraded_nside=32, spectra=False, smooth_mask=False, spectra_mask=False, base_filename="out", root_folder=".", metadata={}):
    """Combine, smooth, take-spectra, write metadata

    The maps (I or IQU) are first combined with their own weights, then smoothed and degraded.
    This function writes a combined smoothed and degraded map, a spectra 1 or 6 components (not degraded) and a json file with metadata
    
    Parameters
    ----------
    maps_and_weights : list of tuples
        [(map1_array, map1_weight), (map2_array, map2_weight), ...]
        each tuple contains a I or IQU map to be combined with its own weight to give the final map
    variance_maps_and_weights : list of tuples
        same as maps_and_weights but containing variances
    fwhm : double
        smoothing gaussian beam width in radians
    degraded_nside : integer
        nside of the output map
    spectra : bool
        whether to compute and write angular power spectra of the combined map
    smooth_mask, spectra_mask : bool array
        masks for smoothing and spectra, same nside of input maps, masks shoud be true *inside* the masked region. spectra are masked with both masks. Typically smooth_mask should be a point source mask, while spectra_mask a galaxy plane mask.
    base_filename : string
        base filename of the output files
    root_folder : string
        root path of the output files
    metadata : dict
        initial state of the metadata to be written to the json files

    Returns
    -------
    None : all outputs are written to fits files
    """

    log.debug("smooth_combine")
    # check if I or IQU
    is_IQU = len(maps_and_weights[0][0]) == 3
    if not is_IQU:
        assert hp.isnpixok(len(maps_and_weights[0][0])), "Input maps must have either 1 or 3 components"

    combined_map = combine_maps(maps_and_weights)
    combined_variance_map = combine_maps(variance_maps_and_weights)

    # apply mask
    for m in combined_map + combined_variance_map:
        m.mask |= smooth_mask

    monopole_I, dipole_I = hp.fit_dipole(combined_map[0], gal_cut=30)
    # remove monopole, only I
    combined_map[0] -= monopole_I

    if spectra:
        # save original masks
        orig_mask = [m.mask.copy() for m in combined_map] 

        # spectra
        log.debug("Anafast")
        for m in combined_map:
            m.mask |= spectra_mask
        # dividing by two in order to recover the same noise as the average map (M1 - M2)/2
        cl = hp.anafast([m/2. for m in combined_map])
        # sky fraction
        sky_frac = (~combined_map[0].mask).sum()/float(len(combined_map[0]))

        if is_IQU:
            for cl_comp in cl:
                cl_comp /= sky_frac
        else:
            cl /= sky_frac

        # write spectra
        log.debug("Write cl: " + base_filename + "_cl.fits")
        try:
            hp.write_cl(os.path.join(root_folder, base_filename + "_cl.fits"), cl)
        except exceptions.NotImplementedError:
            log.error("Write IQU Cls to fits requires more recent version of healpy")
        del cl

        # expected cl from white noise
        # /4. to have same normalization of cl
        metadata["whitenoise_cl"] = utils.get_whitenoise_cl(combined_variance_map[0]/4., mask=combined_map[0].mask) / sky_frac
        if is_IQU:
            # /2. is the mean, /4. is the half difference in power
            metadata["whitenoise_cl_P"] = utils.get_whitenoise_cl((combined_variance_map[1] + combined_variance_map[2])/2./4., mask=combined_map[1].mask | combined_map[2].mask) / sky_frac 

        # restore masks
        for m, mask in zip(combined_map, orig_mask):
            m.mask = mask

    # smooth
    log.debug("Smooth")

    smoothed_map = hp.smoothing(combined_map, fwhm=fwhm)

    if is_IQU:
        smoothed_variance_map = [utils.smooth_variance_map(var, fwhm=fwhm) for var in combined_variance_map]
        for comp,m,var in zip("IQU", smoothed_map, smoothed_variance_map):
             metadata["map_chi2_%s" % comp] = np.mean(m**2 / var) 
        for comp,m,var in zip("IQU", combined_map, combined_variance_map):
             metadata["map_unsm_chi2_%s" % comp] = np.mean(m**2 / var) 
    else:
        smoothed_variance_map = utils.smooth_variance_map(combined_variance_map[0], fwhm=fwhm)
        metadata["map_chi2"] = np.mean(smoothed_map**2 / smoothed_variance_map) 
        metadata["map_unsm_chi2"] = np.mean(combined_map[0]**2 / combined_variance_map[0]) 

    del smoothed_variance_map
    # removed downgrade of variance
    # smoothed_variance_map = hp.ud_grade(smoothed_variance_map, degraded_nside, power=2)

    # fits
    log.info("Write fits map: " + base_filename + "_map.fits")
    smoothed_map = hp.ud_grade(smoothed_map, degraded_nside)
    hp.write_map(os.path.join(root_folder, base_filename + "_map.fits"), smoothed_map)

    # metadata
    metadata["base_file_name"] = base_filename
    metadata["file_name"] = base_filename + "_cl.fits"
    metadata["file_type"] += "_cl"
    metadata["removed_monopole_I"] = monopole_I
    metadata["dipole_I"] = tuple(dipole_I)

    if spectra:
        metadata["sky_fraction"] = sky_frac
        with open(os.path.join(root_folder, base_filename + "_cl.json"), 'w') as f:
            json.dump(metadata, f)

    metadata["file_name"] = base_filename + "_map.fits"
    metadata["file_type"] = metadata["file_type"].replace("_cl","_map")

    metadata["smooth_fwhm_deg"] = "%.2f" % np.degrees(fwhm)
    metadata["out_nside"] = degraded_nside
    if is_IQU:
        for comp,m in zip("IQU", smoothed_map):
             metadata["map_p2p_%s" % comp] = m.ptp()
             metadata["map_std_%s" % comp] = m.std()
    else:
        metadata["map_p2p_I"] = smoothed_map.ptp()
        metadata["map_std_I"] = smoothed_map.std()

    with open(os.path.join(root_folder, base_filename + "_map.json"), 'w') as f:
        json.dump(metadata, f)


def halfrings(freq, ch, surv, pol='I', smooth_combine_config=None, root_folder="out/",log_to_file=False, mapreader=None):
    """Half ring differences
    
    Parameters
    ----------
    freq : integer
        channels frequency
    chtag : string 
        can be "" for frequency, radiometer("LFI18S"), horn("LFI18"), quadruplet("18_23"), detset("detset_1")
    surv : int or string
        "nominal", "full", or survey number
    pol : string
        required polarization components, e.g. 'I', 'Q', 'IQU'
    smooth_combine_config : dict
        configuration for smooth_combine, see its docstring
    root_folder : string
        root path of the output files
    log_to_file : bool
        whether log to file
    """

    try:
        os.makedirs(os.path.join(root_folder, "halfrings"))
    except:
        pass

    if ch:
        chtag = ch
    else:
        chtag = str(freq)

    base_filename = os.path.join("halfrings", "%s_SS%s" % (chtag, str(surv)))
    if log_to_file:
        configure_file_logger(os.path.join(root_folder, base_filename))

    ps_mask, gal_mask = mapreader.read_masks(freq)

    metadata = dict( 
        file_type="halfring_%s" % (reader.type_of_channel_set(ch),),
        channel=chtag,
        survey=surv,
        title="Halfring difference survey %s ch %s" % (str(surv), chtag),
        )
    log.info("Call smooth_combine")
    var_pol = 'A' if len(pol) == 1 else 'ADF' # for I only read sigma_II, else read sigma_II, sigma_QQ, sigma_UU
    smooth_combine(
            [(mapreader(freq, surv, ch, halfring=1, pol=pol), 1), 
             (mapreader(freq, surv, ch, halfring=2, pol=pol), -1)],
            [(mapreader(freq, surv, ch, halfring=1, pol=var_pol), 1.), 
             (mapreader(freq, surv, ch, halfring=2, pol=var_pol), 1.)],
              base_filename=base_filename,
              metadata=metadata,
              root_folder=root_folder,
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
            **smooth_combine_config)
    log.info("Completed")

def surveydiff(freq, ch, survlist=[1,2,3,4,5], pol='I', root_folder="out/", smooth_combine_config=None, log_to_file=False, bp_corr=False, mapreader=None):
    """Survey differences

    for a specific channel or channel set, produces all the possible combinations of the surveys in survlist
    
    Parameters
    ----------
    survlist : list of survey id (1..5, "nominal", "full")

    see the halfrings function for other parameters
    """
    try:
        os.makedirs(os.path.join(root_folder, "surveydiff"))
    except:
        pass

    if ch:
        chtag = ch
    else:
        chtag = str(freq)

    if log_to_file:
        logfilename = "%s_SSdiff" % chtag
        if bp_corr:
            logfilename += "_bpcorr"
        configure_file_logger(os.path.join(root_folder, "surveydiff", logfilename))

    # read all maps
    maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=pol, bp_corr=bp_corr)) for surv in survlist])

    log.debug("Read variance")
    var_pol = 'A' if len(pol) == 1 else 'ADF' # for I only read sigma_II, else read sigma_II, sigma_QQ, sigma_UU
    variance_maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=var_pol, bp_corr=False)) for surv in survlist])
    for var_m in variance_maps.values():
        assert np.all(var_m >= 0)

    log.debug("All maps read")

    ps_mask, gal_mask = mapreader.read_masks(freq)

    log.debug("Metadata")

    metadata = dict( 
        file_type="surveydiff_%s" % (reader.type_of_channel_set(ch),),
        channel=chtag,
        )

    combs = list(itertools.combinations(survlist, 2))
    for comb in combs:
        # in case of even-odd, swap to odd-even. Do the same for
        # combinations like e.g. SS3-SS1 (-> SS1-SS3)
        if (comb[1] % 2 != 0 and comb[0] % 2 == 0) or (comb[1] < comb[0]):
            comb = (comb[1], comb[0])

        metadata["title"]="Survey difference SS%s-SS%s ch %s" % (str(comb[0])[:4], str(comb[1])[:4], chtag)
        base_filename = os.path.join("surveydiff", "%s_SS%d-SS%d" % (chtag, comb[0], comb[1]))
        if bp_corr:
            metadata["title"] += " BPCORR"
            base_filename += "_bpcorr"

        log.debug("Launching smooth_combine")
        smooth_combine(
                [ (maps[comb[0]],  1),
                  (maps[comb[1]], -1) ],
                [ (variance_maps[comb[0]], 1),
                  (variance_maps[comb[1]], 1) ],
                base_filename=base_filename,
                root_folder=root_folder,
                metadata=dict(metadata.items() + {"surveys": comb}.items()),
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
                **smooth_combine_config )
    log.info("Completed")

def chdiff(freq, chlist, surv, pol='I', smooth_combine_config=None, root_folder="out/", log_to_file=False, mapreader=None):
    """Channel difference

    for a specific survey, produces all the possible combinations of the channels in chlist
    
    Parameters
    ----------
    chlist : list of channel tags (see reader or halfrings documentation)

    see the halfrings function for other parameters
    """

    try:
        os.makedirs(os.path.join(root_folder, "chdiff"))
    except:
        pass

    base_filename=os.path.join("chdiff", "%d_SS%s" % (freq, surv))
    if log_to_file:
        configure_file_logger(os.path.join(root_folder, base_filename))

    # read all maps
    maps = dict([(ch, mapreader(freq, surv, ch, halfring=0, pol=pol)) for ch in chlist])

    log.debug("Read variance")
    var_pol = 'A' if len(pol) == 1 else 'ADF' # for I only read sigma_II, else read sigma_II, sigma_QQ, sigma_UU
    variance_maps = dict([(ch, mapreader(freq, surv, ch, halfring=0, pol=var_pol, bp_corr=False)) for ch in chlist])
    for var_m in variance_maps.values():
        assert np.all(var_m >= 0)

    ps_mask, gal_mask = mapreader.read_masks(freq)

    metadata = dict( 
        file_type="chdiff",
        survey=surv
        )

    combs = list(itertools.combinations(chlist, 2))
    for comb in combs:
        metadata["title"]="Channel difference %s-%s SS%s" % (comb[0], comb[1], surv)
        metadata["channel"] = comb
        smooth_combine(
                [ (maps[comb[0]],  1),
                  (maps[comb[1]], -1) ],
                [ (variance_maps[comb[0]], 1),
                  (variance_maps[comb[1]], 1) ],
                base_filename=os.path.join("chdiff", "%s-%s_SS%s" % (comb[0],   comb[1], surv)),
                root_folder=root_folder,
                metadata=metadata,
                smooth_mask=ps_mask,
                spectra_mask=gal_mask,
                **smooth_combine_config )
    log.info("Completed")
