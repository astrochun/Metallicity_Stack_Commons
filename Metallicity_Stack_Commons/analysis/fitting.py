from logging import Logger
from typing import Union
import numpy as np
from astropy.convolution import Box1DKernel, convolve
from astropy.io import ascii as asc

from .. import scalefact, wavelength_dict
from ..logging import log_stdout, log_verbose

# Ratio of OII doublet line
con1 = wavelength_dict['OII_3729'] / wavelength_dict['OII_3726']


def gauss(x: np.ndarray, xbar: float, s: float, a: float, c: float) \
        -> np.ndarray:
    """
    Function providing Gaussian profile for curve_fit

    :param x: Wavelength array
    :param xbar: Central wavelength of Gaussian fit
    :param s: Sigma (width) of Gaussian fit
    :param a: Amplitude of Gaussian fit
    :param c: Continuum constant for Gaussian fit

    :return: Gaussian fit
    """
    return a * np.exp(-(x - xbar) ** 2 / (2 * s ** 2)) + c


def double_gauss(x: np.ndarray, xbar: float, s1: float, a1: float, c: float,
                 s2: float, a2: float) -> np.array:
    """
    Function providing double Gaussian profile (emission and absorption)
    for curve_fit

    :param x: Wavelength array
    :param xbar: Central wavelength of Gaussian fit
    :param s1: Sigma (width) of first Gaussian fit (positive)
    :param a1: Amplitude of first Gaussian fit
    :param c: Continuum constant for Gaussian fit
    :param s2: Sigma (width) of second Gaussian fit (absorption)
    :param a2: Amplitude of second Gaussian fit

    :return: Double Gaussian fit
    """

    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - xbar) ** 2 / (2 * s2 ** 2))


def oxy2_gauss(x: np.ndarray, xbar: float, s1: float, a1: float, c: float,
               s2: float, a2: float) -> np.ndarray:
    """
    Function providing [OII] doublet Gaussian profile for curve_fit

    :param x: Wavelength array
    :param xbar: Central wavelength of [OII]3726 Gaussian fit
    :param s1: Sigma (width) of [OII]3726 Gaussian fit
    :param a1: Amplitude of [OII]3726 Gaussian fit
    :param c: Continuum constant for Gaussian fit
    :param s2: Sigma (width) of [OII]3728 Gaussian fit
    :param a2: Amplitude of [OII]3728 Gaussian fit

    :return: [OII] doublet Gaussian fit
    """

    return a1 * np.exp(-(x - xbar) ** 2 / (2 * s1 ** 2)) + c + \
           a2 * np.exp(-(x - (xbar * con1)) ** 2 / (2 * s2 ** 2))


def rms_func(wave: np.ndarray, dispersion: float, lambda_in: float,
             y0: np.ndarray, sigma_array: float, mask_flag: np.ndarray,
             verbose: bool = False, log: Logger = log_stdout()):
    """
    Compute rms in the spectra

    :param wave: Array of rest wavelengths
    :param dispersion: Spectral dispersion in AA/pix
    :param lambda_in: Central wavelength of fit
    :param y0: Array of fluxes in units of erg/s/cm2/AA
    :param sigma_array: Gaussian sigma (AA)
    :param mask_flag: Indicates spectra are masked for OH skyline contamination
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: logging.Logger object

    :return:
    """

    log_verbose(log, "starting ...", verbose=verbose)

    x_idx = np.where((np.abs(wave - lambda_in) <= 100) & (mask_flag == 0))[0]

    sigma = np.std(y0[x_idx])

    if sigma_array == 0:
        RMS = sigma/scalefact

        log_verbose(log, "finished.", verbose=verbose)
        return RMS
    else:
        pix = 5 * sigma_array / dispersion
        s_pix = np.sqrt(pix)
        ini_sig = s_pix * sigma * dispersion
        RMS_pix = sigma * dispersion / scalefact

        log_verbose(log, "finished.", verbose=verbose)
        return ini_sig / scalefact, RMS_pix


def OIII4363_flux_limit(combine_flux_file: str, verbose: bool = False,
                        log: Logger = log_stdout()) -> \
        Union[None, np.ndarray]:
    """
    Determine 3-sigma limit on [OIII]4363 based on H-gamma measurements

    :param combine_flux_file: Filename of ASCII file containing emission-line
                              flux measurements
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: logging.Logger object

    :return: Array containing 3-sigma flux limit
    """

    log_verbose(log, "starting ...", verbose=verbose)

    try:
        combine_fits = asc.read(combine_flux_file)
    except FileNotFoundError:
        log.warning(f"File not found! {combine_flux_file}")
        return

    Hgamma    = combine_fits['HGAMMA_Flux_Observed'].data
    Hgamma_SN = combine_fits['HGAMMA_S/N'].data

    flux_limit = (Hgamma / Hgamma_SN) * 3

    log_verbose(log, "finished.", verbose=verbose)
    return flux_limit


def movingaverage_box1D(values: np.ndarray, width: float,
                        boundary: str = 'fill', fill_value: float = 0.0)\
        -> np.ndarray:
    """
    Applies as boxcar kernel to smooth spectra

    :param values: Array containing the spectrum
    :param width: Width for smoothing
    :param boundary: Handling of boundary values.
                     Options are: 'None', 'fill', 'wrap', and 'extend'
                     See astropy.convolution.convolve for more information
    :param fill_value: Indicate fill value for default boundary='fill'

    :return: Array contained the smoothed/convolved spectrum
    """

    box_kernel = Box1DKernel(width)
    smooth = convolve(values, box_kernel, boundary=boundary,
                      fill_value=fill_value)
    return smooth
