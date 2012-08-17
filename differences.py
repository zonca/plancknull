import os
import json
import exceptions
import itertools
import numpy as np
from collections import SortedDict
import logging as log

import healpy as hp

def smooth_combine(maps_and_weights, fwhm=np.radians(2.0), degraded_nside=32, spectra=False, mask=False, plot=True, base_filename="smooth_output/out"):
    """Combine, smooth, take-spectra, write metadata"""

    # check if I or IQU
    is_IQU = len(maps_and_weights[0][0]) == 3
    if not is_IQU:
        assert hp.isnpixok(len(maps_and_weights[0][0])), "Input maps must have either 1 or 3 components"

    # combine
    if is_IQU:
        combined_map = [
                np.ma.sum(
                    [m[comp]*w for m,w in maps_and_weights]
                    , axis=1)
                for comp in [0,1,2]
                ]
    else: # single component only
        combined_map = [np.ma.sum(
            [m*w for m,w in maps_and_weights]
            , axis=1)]

    # apply mask
    for m in combined_map:
        m.mask |= mask

    # remove monopole, only I
    combined_map[0] -= hp.fit_monopole(combined_map[0].filled(), gal_cut=30)

    # spectra
    cl, alms = hp.anafast([m.filled() for m in combined_map], alm=True)
    # write spectra
    try:
        hp.write_cl(base_filename + "_cl.fits", cl)
    except exceptions.NotImplementedError:
        log.error("Write IQU Cls to fits requires more recent version of healpy")
    del cl

    # smooth
    hp.smoothalm(alms, fwhm=fwhm, inplace=True) # inplace!
    smoothed_map = hp.alm2map(alms, degraded_nside, pixwin = False)
    del alms

    # fits
    hp.write_map(base_filename + "_map.fits", smoothed_map)

    # metadata
    smoothed_map = hp.ma(smoothed_map)
    metadata = SortedDict([
             ("smooth_fwhm" , "%.2f" % np.degrees(fwhm)),
             ("out_nside" , degraded_nside),
             ("map_p2p" , smoothed_map.ptp()),
             ("map_std" , smoothed_map.std()),
    ])

    with open(base_filename + ".json") as f:
        json.dump(metadata, f)

def halfrings(freq, ch, surv, pol='I', smooth_combine_config=None, degraded_nside=None, mapreader=None, output_folder="halfrings/"):
    """Half ring differences"""
    try:
        os.mkdir(output_folder)
    except:
        pass

    smooth_combine(
            [(mapreader(freq, surv, ch, halfring=1, pol=pol), .5), 
             (mapreader(freq, surv, ch, halfring=2, pol=pol), -.5)],
              base_filename=os.path.join(output_folder, "%s_SS%s" % (ch, str(surv))),
            **smooth_combine_config)

def surveydiff(freq, ch, survlist=[1,2,3,4,5], pol='I', output_folder="survdiff/", mapreader=None, smooth_combine_config=None):
    """Survey differences"""
    try:
        os.mkdir(output_folder)
    except:
        pass

    # read all maps
    maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=pol)) for surv in survlist])

    combs = list(itertools.combinations(survlist, 2))
    for comb in combs:
        # in case of even-odd, swap to odd-even
        if comb[1] % 2 != 0 and comb[0] % 2 == 0:
            comb = (comb[1], comb[0])

        smooth_combine(
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join(output_folder, "%s_SS%d-SS%d" % (ch, comb[0], comb[1])) ,
                **smooth_combine_config )

def chdiff(freq, chlist, surv, pol='I', mapreader=None, smooth_combine_config=None, output_folder="chdiff/"):
    try:
        os.mkdir(output_folder)
    except:
        pass
    # read all maps
    maps = dict([(surv, mapreader(freq, surv, ch, halfring=0, pol=pol)) for ch in chlist])

    combs = list(itertools.combinations(chlist, 2))
    for comb in combs:
        smooth_combine(
                [ (maps[comb[0]],  .5),
                  (maps[comb[1]], -.5) ],
                base_filename=os.path.join(output_folder, "%s_SS%d-SS%d" % (ch, comb[0], comb[1])) ,
                **smooth_combine_config )

if __name__ == "__main__":
    pass
