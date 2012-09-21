import numpy as np
import healpy as hp

def smooth_variance_map(var_m, orig_fwhm, fwhm):
    """Smooth a variance map

    Algorithm from 'Pixel errors in convolved maps'
    J.P. Leahy, version 0.2

    Parameters
    ----------
    var_m : array
        input variance map
    orig_fwhm : float (radians)
        original fwhm
    fwhm : float (radians)
        target fwhm

    Returns
    -------
    smoothed_var_m : array
        smoothed variance map
    """

    # smooth map
    fwhm_variance = fwhm / np.sqrt(2)
    smoothed_var_m = hp.smoothing(var_m, fwhm=fwhm_variance)

    # normalization factor
    pix_area = hp.nside2pixarea(hp.npix2nside(len(var_m)))
    orig_beam_width = orig_fwhm/np.sqrt(8*np.log(2))
    A_vb = pix_area / (4. * np.pi * orig_beam_width**2)
    smoothed_var_m *= A_vb

    return smoothed_var_m
