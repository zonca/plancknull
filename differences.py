import os
import json
import exceptions
import itertools
import numpy as np
import logging as log
import healpy as hp
import reader

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

def smooth_combine(maps_and_weights, fwhm=np.radians(2.0), degraded_nside=32, spectra=False, smooth_mask=False, spectra_mask=False, base_filename="out", root_folder=".", metadata={}):
    """Combine, smooth, take-spectra, write metadata

    The maps (I or IQU) are first combined with their own weights, then smoothed and degraded.
    This function writes a combined smoothed and degraded map, a spectra 1 or 6 components (not degraded) and a json file with metadata
    
    Parameters
    ----------
    maps_and_weights : list of tuples
        [(map1_array, map1_weight), (map2_array, map2_weight), ...]
        each tuple contains a I or IQU map to be combined with its own weight to give the final map
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

    # apply mask
    for m in combined_map:
        m.mask |= smooth_mask

    monopole_I, dipole_I = hp.fit_dipole(combined_map[0], gal_cut=30)
    # remove monopole, only I
    combined_map[0] -= monopole_I

    # smooth
    log.debug("Smooth")

    smoothed_map = hp.smoothing(combined_map, fwhm=fwhm)
    smoothed_map = hp.ud_grade(smoothed_map, degraded_nside)

    # fits
    log.info("Write fits map: " + base_filename + "_map.fits")
    hp.write_map(os.path.join(root_folder, base_filename + "_map.fits"), smoothed_map)

    # metadata
    metadata["base_file_name"] = base_filename
    metadata["file_name"] = base_filename + "_cl.fits"
    metadata["file_type"] += "_cl"
    metadata["removed_monopole_I"] = monopole_I
    metadata["dipole_I"] = tuple(dipole_I)

    if spectra:
        # spectra
        log.debug("Anafast")
        for m in combined_map:
            m.mask |= spectra_mask
        cl = hp.anafast(combined_map)
        # write spectra
        log.debug("Write cl: " + base_filename + "_cl.fits")
        try:
            hp.write_cl(os.path.join(root_folder, base_filename + "_cl.fits"), cl)
        except exceptions.NotImplementedError:
            log.error("Write IQU Cls to fits requires more recent version of healpy")
        del cl

        with open(os.path.join(root_folder, base_filename + "_cl.json"), 'w') as f:
            json.dump(metadata, f)

    metadata["file_name"] = base_filename + "_map.fits"
    metadata["file_type"] = metadata["file_type"].replace("_cl","_map")

    metadata["smooth_fwhm_deg"] = "%.2f" % np.degrees(fwhm)
    metadata["out_nside"] = degraded_nside
    if is_IQU:
        for comp,m in zip("IQU", smoothed_map):
             metadata["map_p2p_%s" % comp] = "%.2e" % m.ptp()
             metadata["map_std_%s" % comp] = "%.2e" % m.std()
    else:
        metadata["map_p2p"] = "%.2e" % smoothed_map.ptp()
        metadata["map_std"] = "%.2e" % smoothed_map.std()

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
    smooth_combine(
            [(mapreader(freq, surv, ch, halfring=1, pol=pol), .5), 
             (mapreader(freq, surv, ch, halfring=2, pol=pol), -.5)],
              base_filename=base_filename,
              metadata=metadata,
              root_folder=root_folder,
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
            **smooth_combine_config)

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
        configure_file_logger(os.path.join(root_folder, "surveydiff", "%s_SSdiff" % chtag))

    # read all maps
    maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=pol, bp_corr=bp_corr)) for surv in survlist])

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
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=base_filename,
                root_folder=root_folder,
                metadata=dict(metadata.items() + {"surveys": comb}.items()),
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
                **smooth_combine_config )

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
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join("chdiff", "%s-%s_SS%s" % (comb[0],   comb[1], surv)),
                root_folder=root_folder,
                metadata=metadata,
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
                **smooth_combine_config )
