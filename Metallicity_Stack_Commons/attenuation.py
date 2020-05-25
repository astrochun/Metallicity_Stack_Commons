from astropy.io import ascii as asc
from astropy.table import Table
import numpy as np

from . import k_dict


HgHb_CaseB = 0.468  # Hg/Hb ratio for zero reddening

k_HBETA  = k_dict['HBETA']
k_HGAMMA = k_dict['HGAMMA']


def compute_EBV(fitspath, combine_asc):
    """
    Purpose:
      Determines E(B-V) from Hg/Hb flux ratio

    :param fitspath: str containing root path
    :param combine_asc: Astropy table containing emission-line flux
    """

    ID = combine_asc['ID']
    HBETA  = combine_asc['HBETA_Flux_Observed'].data
    HGAMMA = combine_asc['HGAMMA_Flux_Observed'].data

    EBV = -2.5 * np.log10((HGAMMA / HBETA) / HgHb_CaseB) / (k_HGAMMA - k_HBETA)

    out_ascii = fitspath + '/dust_attenuation_values.tbl'

    tab1 = Table([ID, EBV], names=('ID', 'E(B-V)'))
    asc.write(tab1, out_ascii, format='fixed_width_two_line')


def compute_A(EBV):

    k_arr  = np.array(list(k_dict.values()))

    A_arr  = k_arr * EBV
    A_dict = dict(zip(list(k_dict.keys()), A_arr))

    return A_dict


def line_ratio_atten(ratio, EBV, wave_top, wave_bottom):

    k_top    = k_dict[wave_top]
    k_bottom = k_dict[wave_bottom]

    ratio_atten = ratio * 10**(0.4*EBV*(k_top - k_bottom))

    return ratio_atten


def Hb_SFR(log_LHb, EBV):
    """
    Purpose:
      Determine dust-corrected SFR using the H-beta luminosity and a
      measurement for nebular attenuation

    Equation below is based on Eq. 2 in Ly et al. (2015), ApJ, 805, 45
      DOI: https://doi.org/10.1088/0004-637X/805/1/45

    :param log_LHb: numpy array or float containing logarithm of H-beta
           luminosity in units of erg/s
    :param EBV: numpy array or float providing E(B-V)

    :return logSFR: numpy array or float containing the SFR in
            logarithmic units of M_sun/yr
    """

    logSFR = np.log10(4.4e-42 * 2.86) + 0.4*EBV*k_dict['HBETA'] + log_LHb

    return logSFR
