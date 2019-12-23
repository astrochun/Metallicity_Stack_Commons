import numpy as np

from . import scalefact

con1 = 3728.91 / 3726.16 # Ratio of OII doublet line


def gauss(x, xbar, s, a, c):
    '''

    Function providing Gaussian profile for curve_fit

    :param x: wavelength array
    :param xbar: central wavelength of Gaussian fit
    :param s: sigma (width) of Gaussian fit
    :param a: amplitude of Gaussian fit
    :param c: continuum constant for Gaussian fit

    :return:
        Gaussian fit

    '''
    return a * np.exp(-(x - xbar) ** 2 / (2 * s ** 2)) + c


def double_gauss(x, xbar, s1, a1, c, s2, a2):
    '''

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
    '''

    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - xbar) ** 2 / (2 * s2 ** 2))


def oxy2_gauss(x, xbar, s1, a1, c, s2, a2):
    '''

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

    '''

    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - (xbar * con1)) ** 2 / (2 * s2 ** 2))


def rms_func(wave, dispersion, lambda_in, y0, sigma_array, mask_flag):
    '''

    Purpose:
        Compute rms in the spectra

    :param wave: numpy array containing rest wavelengths
    :param dispersion: spectral dispersion in AA/pix
    :param lambda_in: central wavelength of fit
    :param y0: numpy array containing spectra in units of erg/s/cm2/AA
    :param sigma_array: Gaussian sigma
    :param scalefact: Spectra normalization factor
    :param mask_flag: numpy array indicating spectra that are masked for OH
           skyline contamintion

    :return:
    '''

    x_idx = np.where((np.abs(wave - lambda_in) <= 100) & (mask_flag == 0))[0]

    sigma = np.std(y0[x_idx])
    pix = 5 * sigma_array / dispersion
    s_pix = np.sqrt(pix)

    ini_sig = s_pix * sigma * dispersion
    RMS_pix = sigma * dispersion / scalefact

    return ini_sig / scalefact, RMS_pix
