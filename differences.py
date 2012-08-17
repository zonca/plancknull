import itertools
import numpy as np

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
                    [m[comp]*w for m[comp],w in maps_and_weights]
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
    cl, alm = hp.anafast([m.filled() for m in combined_map], alm=True)
    # plot spectra

    # smooth
    hp.smoothalm(alm, fwhm=fwhm, inplace=True) # inplace!
    smoothed_map = hp.alm2map(alms, degraded_nside, pixwin = False)

    # fits
    hp.write_map(base_filename + "_map.fits", smoothed_map)
    # plot map

    # metadata


def halfrings(freq, ch, surv, pol='I', smooth_combine_config=None, degraded_nside=None, mapreader=None):
    """Half ring differences"""

    smooth_combine(
            [(mapreader(freq, surv, ch, halfring=1, pol=pol), .5), 
             (mapreader(freq, surv, ch, halfring=2, pol=pol), -.5)],
            **smooth_combine_config)

def surveydiff(freq, ch, survlist=[1,2,3,4,5], pol='I', output_folder="survdiff/", mapreader=None, smooth_combine_config=None):
    """Survey differences"""

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
                **smooth_combine_config )

def chdiff(chlist, surv, pol='I'):

if __name__ == "__main__":
    pass
