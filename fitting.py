import numpy as np

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