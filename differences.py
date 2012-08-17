def smooth_combine(maps_and_weights, fwhm=np.radians(2.0), degraded_nside=32, spectra=False):
    """Combine, smooth, take-spectra, write metadata"""
    pass

def halfrings(ch, surv, pol='I', smooth_combine_config=None, degraded_nside=None, mapreader=None):
    """Half ring differences"""

    smooth_combine(
            [(mapreader(freq, surv, chtag, halfring=1, pol=pol), .5), 
             (mapreader(freq, surv, chtag, halfring=2, pol=pol), -.5)],
            **smooth_combine_config)
        

def surveydiff(ch, survlist=[1,2,3,4,5], pol='I'):
    pass

def chdiff(chlist, surv, pol='I'):
    pass


if __name__ == "__main__":
    pass
