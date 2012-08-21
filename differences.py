import os
import json
import exceptions
import itertools
import numpy as np
import logging as log
from glob import glob

import healpy as hp

def smooth_combine(maps_and_weights, fwhm=np.radians(2.0), degraded_nside=32, spectra=False, smooth_mask=False, spectra_mask=False, plot=True, base_filename="smooth_output/out", metadata={}):
    """Combine, smooth, take-spectra, write metadata"""

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

    # remove monopole, only I
    combined_map[0] -= hp.fit_monopole(combined_map[0], gal_cut=30)

    # smooth
    log.debug("Smooth")

    smoothed_map = hp.smoothing(combined_map, fwhm=fwhm)
    smoothed_map = hp.ud_grade(smoothed_map, degraded_nside)

    # fits
    log.debug("Write fits map")
    hp.write_map(base_filename + "_map.fits", smoothed_map)

    # spectra
    log.debug("Anafast")
    for m in combined_map:
        m.mask |= spectra_mask
    cl = hp.anafast(combined_map)
    # write spectra
    try:
        hp.write_cl(base_filename + "_cl.fits", cl)
    except exceptions.NotImplementedError:
        log.error("Write IQU Cls to fits requires more recent version of healpy")
    del cl

    # metadata
    metadata["smooth_fwhm_deg"] = "%.2f" % np.degrees(fwhm)
    metadata["out_nside"] = degraded_nside
    metadata["base_file_name"] = base_filename

    if is_IQU:
        for comp,m in zip("IQU", smoothed_map):
             metadata["map_p2p_%s" % comp] = "%.2e" % m.ptp()
             metadata["map_std_%s" % comp] = "%.2e" % m.std()
    else:
        metadata["map_p2p"] = "%.2e" % smoothed_map.ptp()
        metadata["map_std"] = "%.2e" % smoothed_map.std()

    with open(base_filename + ".json", 'w') as f:
        json.dump(metadata, f)

def type_of_channel_set(ch):
    """Returns a string that identifies the set of channels"""
    if ch == "":
        return "frequency"
    elif ch.find('_') >= 0:
        return "detset"
    else:
        return "single_ch"

def halfrings(freq, ch, surv, pol='I', smooth_combine_config=None, degraded_nside=None, mapreader=None, output_folder="halfrings/", read_masks=None):
    """Half ring differences"""
    try:
        os.mkdir(output_folder)
    except:
        pass

    if ch:
        chtag = ch
    else:
        chtag = str(freq)

    if not read_masks is None:
        ps_mask, gal_mask = read_masks(freq)
    else:
        ps_mask = None; gal_mask = None

    metadata = dict( 
        file_type="halfring_%s" % (type_of_channel_set(ch),),
        channel=chtag,
        title="Halfring difference survey %s ch %s" % (str(surv), chtag),
        )
    smooth_combine(
            [(mapreader(freq, surv, ch, halfring=1, pol=pol), .5), 
             (mapreader(freq, surv, ch, halfring=2, pol=pol), -.5)],
              base_filename=os.path.join(output_folder, "%s_SS%s" % (chtag, str(surv))),
              metadata=metadata,
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
            **smooth_combine_config)
    return 1

def surveydiff(freq, ch, survlist=[1,2,3,4,5], pol='I', output_folder="survdiff/", mapreader=None, smooth_combine_config=None, read_masks=None):
    """Survey differences"""
    try:
        os.mkdir(output_folder)
    except:
        pass

    # read all maps
    maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=pol)) for surv in survlist])

    if ch:
        chtag = ch
    else:
        chtag = str(freq)

    if not read_masks is None:
        ps_mask, gal_mask = read_masks(freq)
    else:
        ps_mask = None; gal_mask = None

    metadata = dict( 
        file_type="surveydiff_%s" % (type_of_channel_set(ch),),
        channel=chtag,
        )

    combs = list(itertools.combinations(survlist, 2))
    for comb in combs:
        # in case of even-odd, swap to odd-even
        if comb[1] % 2 != 0 and comb[0] % 2 == 0:
            comb = (comb[1], comb[0])

        metadata["title"]="Survey difference SS%s-SS%s ch %s" % (str(comb[0])[:4], str(comb[1])[:4], chtag),
        smooth_combine(
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join(output_folder, "%s_SS%d-SS%d" % (chtag, comb[0], comb[1])) ,
                metadata=metadata,
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
                **smooth_combine_config )
    return 1

def chdiff(freq, chlist, surv, pol='I', mapreader=None, smooth_combine_config=None, output_folder="chdiff/", read_masks=None):
    try:
        os.mkdir(output_folder)
    except:
        pass
    # read all maps
    maps = dict([(ch, mapreader(freq, surv, ch, halfring=0, pol=pol)) for ch in chlist])

    if not read_masks is None:
        ps_mask, gal_mask = read_masks(freq)
    else:
        ps_mask = None; gal_mask = None

    metadata = dict( 
        file_type="chdiff",
        )

    combs = list(itertools.combinations(chlist, 2))
    for comb in combs:
        metadata["title"]="Channel difference %s-%s SS%s" % (comb[0], comb[1], surv)
        metadata["channel"] = comb
        smooth_combine(
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join(output_folder, "%s-%s_SS%s" % (comb[0],   comb[1], surv)),
                metadata=metadata,
              smooth_mask=ps_mask,
              spectra_mask=gal_mask,
                **smooth_combine_config )
    return 1

if __name__ == "__main__":
    log.root.level = log.DEBUG
    from reader import SingleFolderDXReader
    freq = 30
    NSIDE = 1024

    ps_mask = np.logical_not(np.floor(hp.ud_grade( 
    hp.read_map(
        glob(os.path.join(os.environ["DX9_LFI"], "MASKs",'mask_ps_%dGHz_*.fits' % freq))[0]), NSIDE))
    )
    gal_mask = np.logical_not(hp.read_map(
        glob(os.path.join(os.environ["DX9_LFI"], "MASKs",'destriping_mask_%d.fits' % freq))[0])
        )

    mapreader = SingleFolderDXReader(os.environ["DX9_LFI"])
    smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128,smooth_mask=ps_mask, spectra_mask=gal_mask)
    halfrings(30, "", "nominal", pol='I', smooth_combine_config=smooth_combine_config, mapreader=mapreader, output_folder="halfrings/")
    #surveydiff(30, "", survlist=[1,2], pol='I', output_folder="survdiff/", mapreader=mapreader, smooth_combine_config=smooth_combine_config)
