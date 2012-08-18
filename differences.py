import os
import json
import exceptions
import itertools
import numpy as np
import logging as log
from glob import glob

import healpy as hp

def smooth_combine(maps_and_weights, fwhm=np.radians(2.0), degraded_nside=32, spectra=False, ps_mask=False, gal_mask=False, plot=True, base_filename="smooth_output/out"):
    """Combine, smooth, take-spectra, write metadata"""

    # check if I or IQU
    is_IQU = len(maps_and_weights[0][0]) == 3
    if not is_IQU:
        assert hp.isnpixok(len(maps_and_weights[0][0])), "Input maps must have either 1 or 3 components"

    # combine
    if is_IQU:
        combined_map = [maps_and_weights[0][0][comp] * maps_and_weights[0][1][comp] for comp in range(3)]  
        for comp in range(3):
            for m,w in maps_and_weights[1:]:
                combined_map[comp] += m[comp] * w
    else: # single component only
        # combined_map is a list of 1 element
        combined_map = [maps_and_weights[0][0] * maps_and_weights[0][1]]
        for m,w in maps_and_weights[1:]:
            combined_map[0] += m * w

    hp.mollview(combined_map[0].filled(), min=-1e-3, max=1e-3)
    # apply mask
    for m in combined_map:
        m.mask |= ps_mask

    # remove monopole, only I
    combined_map[0] -= hp.fit_monopole(combined_map[0].filled(), gal_cut=30)

    # smooth
    log.debug("Smooth")
    alms = hp.map2alm([m.filled() for m in combined_map])
    hp.smoothalm(alms, fwhm=fwhm, inplace=True) # inplace!
    smoothed_map = hp.alm2map(alms, degraded_nside, pixwin = False)
    if is_IQU:
        for sm,m in zip(smoothed_map, combined_map):
            sm[m.mask] = hp.UNSEEN
    else:
        smoothed_map[combined_map[0].mask]=hp.UNSEEN
        
    del alms

    # fits
    log.debug("Write fits map")
    hp.write_map(base_filename + "_map.fits", smoothed_map)

    # spectra
    log.debug("Anafast")
    for m in combined_map:
        m.mask |= gal_mask
    cl = hp.anafast([m.filled() for m in combined_map])
    # write spectra
    try:
        hp.write_cl(base_filename + "_cl.fits", cl)
    except exceptions.NotImplementedError:
        log.error("Write IQU Cls to fits requires more recent version of healpy")
    del cl

    # metadata
    metadata = dict([
             ("smooth_fwhm_deg" , "%.2f" % np.degrees(fwhm)),
             ("out_nside" , degraded_nside),
    ])

    if is_IQU:
        smoothed_map = [hp.ma(m) for m in smoothed_map]
        for comp,m in zip("IQU", smoothed_map):
             metadata["map_p2p_%s" % comp] = "%.2e" % m.ptp()
             metadata["map_std_%s" % comp] = "%.2e" % m.std()
    else:
        metadata["map_p2p"] = "%.2e" % smoothed_map.ptp()
        metadata["map_std"] = "%.2e" % smoothed_map.std()

    with open(base_filename + ".json", 'w') as f:
        json.dump(metadata, f)

def halfrings(freq, ch, surv, pol='I', smooth_combine_config=None, degraded_nside=None, mapreader=None, output_folder="halfrings/"):
    """Half ring differences"""
    try:
        os.mkdir(output_folder)
    except:
        pass

    if ch:
        chtag = ch
    else:
        chtag = str(freq)

    smooth_combine(
            [(mapreader(freq, surv, ch, halfring=1, pol=pol), .5), 
             (mapreader(freq, surv, ch, halfring=2, pol=pol), -.5)],
              base_filename=os.path.join(output_folder, "%s_SS%s" % (chtag, str(surv))),
            **smooth_combine_config)

def surveydiff(freq, ch, survlist=[1,2,3,4,5], pol='I', output_folder="survdiff/", mapreader=None, smooth_combine_config=None):
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

    combs = list(itertools.combinations(survlist, 2))
    for comb in combs:
        # in case of even-odd, swap to odd-even
        if comb[1] % 2 != 0 and comb[0] % 2 == 0:
            comb = (comb[1], comb[0])

        smooth_combine(
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join(output_folder, "%s_SS%d-SS%d" % (chtag, comb[0], comb[1])) ,
                **smooth_combine_config )

def chdiff(freq, chlist, surv, pol='I', mapreader=None, smooth_combine_config=None, output_folder="chdiff/"):
    try:
        os.mkdir(output_folder)
    except:
        pass
    # read all maps
    maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=pol)) for ch in chlist])

    if ch:
        chtag = ch
    else:
        chtag = str(freq)

    combs = list(itertools.combinations(chlist, 2))
    for comb in combs:
        smooth_combine(
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join(output_folder, "%s_SS%d-SS%d" % (chtag, comb[0], comb[1])) ,
                **smooth_combine_config )

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
    smooth_combine_config = dict(fwhm=np.radians(1.), degraded_nside=128,ps_mask=ps_mask, gal_mask=gal_mask)
    #halfrings(30, "", "nominal", pol='I', smooth_combine_config=smooth_combine_config, mapreader=mapreader, output_folder="halfrings/")
    surveydiff(30, "", survlist=[1,2], pol='I', output_folder="survdiff/", mapreader=mapreader, smooth_combine_config=smooth_combine_config)
