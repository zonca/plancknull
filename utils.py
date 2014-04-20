import numpy as np
import healpy as hp
import logging as log
import exceptions

HORNS = {30:[27,28], 44:[24,25,26], 70:list(range(18,23+1))}

def chlist(freq):
    try:
        from planck.Planck import Planck
        pl = Planck()
        return [ch.tag for ch in pl.f[freq].ch]
    except exceptions.ImportError:
        horns = HORNS[freq]
        chs = []
        for horn in horns:
            chs += ["LFI%dM" % horn, "LFI%dS" % horn]
        return chs

def get_chisq(m, var):
    return np.mean(m**2/var)

def get_whitenoise_cl(var, mask):
    """White noise C_ell

    Computes the C_ell's of variance map as
    mean variance multiplied by the pixel area
    
    Parameters
    ----------
    var : array
        variance map, 1 component only
    mask : array
        mask to be applied, True or 1 if pixel IS masked
    """
    log.info("Masked pixels: %d" % mask.sum())
    return (var.filled() * ~mask).mean() * 4 * np.pi / len(var)

def smooth_variance_map(var_m, fwhm):
    """Smooth a variance map

    Algorithm from 'Pixel errors in convolved maps'
    J.P. Leahy, version 0.2

    Parameters
    ----------
    var_m : array
        input variance map
    fwhm : float (radians)
        target fwhm

    Returns
    -------
    smoothed_var_m : array
        smoothed variance map
    """

    # smooth map
    fwhm_variance = fwhm / np.sqrt(2)
    smoothed_var_m = hp.smoothing(var_m, fwhm=fwhm_variance, regression=False)

    # normalization factor
    pix_area = hp.nside2pixarea(hp.npix2nside(len(var_m)))
    orig_beam_width = fwhm/np.sqrt(8*np.log(2))
    A_vb = pix_area / (4. * np.pi * orig_beam_width**2)
    smoothed_var_m *= A_vb

    return smoothed_var_m

def read_mask(filename, nside):
    return np.logical_not(
               np.floor(
                    hp.ud_grade(
                        hp.read_map(filename), nside_out=nside
                    )
               ).astype(np.bool)
           )
