import numpy as np
from astropy.convolution import Box1DKernel, convolve
from astropy.io import ascii as asc

from .. import scalefact, wavelength_dict
from ..logging import log_stdout

con1 = wavelength_dict['OII_3729'] / wavelength_dict['OII_3726']  # Ratio of OII doublet line


def gauss(x, xbar, s, a, c):
    """
    Purpose:
      Function providing Gaussian profile for curve_fit

    :param x: wavelength array
    :param xbar: central wavelength of Gaussian fit
    :param s: sigma (width) of Gaussian fit
    :param a: amplitude of Gaussian fit
    :param c: continuum constant for Gaussian fit

    :return:
        Gaussian fit
    """
    return a * np.exp(-(x - xbar) ** 2 / (2 * s ** 2)) + c


def double_gauss(x, xbar, s1, a1, c, s2, a2):
    """
    Purpose:
      Function providing double Gaussian profile (emission and absorption)
      for curve_fit

    :param x: wavelength array
    :param xbar: central wavelength of Gaussian fit
    :param s1: sigma (width) of first Gaussian fit
    :param a1: amplitude of first Gaussian fit
    :param c: continuum constant for Gaussian fit

    :param s2: sigma (width) of second Gaussian fit
    :param a2: amplitude of second Gaussian fit

    :return:
        Double Gaussian fit
    """

    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - xbar) ** 2 / (2 * s2 ** 2))


def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    """
    Purpose:
      Function providing [OII] doublet Gaussian profile for curve_fit

    :param x:  wavelength array
    :param xbar: central wavelength of [OII]3726 Gaussian fit
    :param s1: sigma (width) of [OII]3726 Gaussian fit
    :param a1: amplitude of [OII]3726 Gaussian fit
    :param c: continuum constant for Gaussian fit
    :param s2: sigma (width) of [OII]3728 Gaussian fit
    :param a2: amplitude of [OII]3728 Gaussian fit

    :return:
            [OII] doublet Gaussian fit

    """

    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - (xbar * con1)) ** 2 / (2 * s2 ** 2))


def rms_func(wave, dispersion, lambda_in, y0, sigma_array, mask_flag):
    """
    Purpose:
      Compute rms in the spectra

    :param wave: numpy array containing rest wavelengths
    :param dispersion: spectral dispersion in AA/pix
    :param lambda_in: central wavelength of fit
    :param y0: numpy array containing spectra in units of erg/s/cm2/AA
    :param sigma_array: Gaussian sigma
    :param mask_flag: numpy array indicating spectra that are masked for OH
           skyline contamintion

    :return:
    """

    x_idx = np.where((np.abs(wave - lambda_in) <= 100) & (mask_flag == 0))[0]

    sigma = np.std(y0[x_idx])

    if sigma_array == 0:
        RMS = sigma/scalefact
        return RMS
    else:
        pix = 5 * sigma_array / dispersion
        s_pix = np.sqrt(pix)
        ini_sig = s_pix * sigma * dispersion
        RMS_pix = sigma * dispersion / scalefact

        return ini_sig / scalefact, RMS_pix


def OIII4363_flux_limit(combine_flux_ascii, log=None):
    """
    Purpose:
      Determine 3-sigma limit on [OIII]4363 based on H-gamma measurements

    :param combine_flux_ascii: filename of ASCII file containing emission-line
                               flux measurements
    :param log: LogClass object

    :return: numpy array containing 3-sigma flux limit
    """

    if log is None:
        log = log_stdout()

    try:
        combine_fits = asc.read(combine_flux_ascii)
    except FileNotFoundError:
        log.info("File not found! "+combine_flux_ascii)
        return

    Hgamma    = combine_fits['HGAMMA_Flux_Observed'].data
    Hgamma_SN = combine_fits['HGAMMA_S/N'].data

    OIII4363_flux_limit = (Hgamma / Hgamma_SN) * 3

    return OIII4363_flux_limit


def movingaverage_box1D(values, width, boundary='fill', fill_value=0.0):
    """
    Purpose:
      Applies as boxcar kernel to smooth spectra

    :param values: numpy array containing the spectrum
    :param width: float containing width for smoothing
    :param boundary: handling of boundary values.
                     Options are: 'None', 'fill', 'wrap', and 'extend'
                     See astropy.convolution.convolve for more information
    :param fill_value: float to indicate fill value for default boundary='fill'

    :return smooth: numpy array contained the smoothed/convolved spectrum
    """

    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary, fill_value=fill_value)
    return smooth
