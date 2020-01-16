from astropy.io import ascii as asc
from astropy.table import Table
import numpy as np

from . import k_dict

HgHb_CaseB = 0.468 # Hg/Hb ratio for zero reddening

k_HBETA  = k_dict['HBETA']
k_HGAMMA = k_dict['HGAMMA']

def compute_EBV(fitspath, combine_asc):

    ID = combine_asc['ID']
    HBETA  = combine_asc['HBETA_Flux_Observed'].data
    HGAMMA = combine_asc['HGAMMA_Flux_Observed'].data

    EBV = -2.5 * np.log10((HGAMMA / HBETA) / (HgHb_CaseB)) / (k_HGAMMA - k_HBETA)

    out_ascii = fitspath + '/dust_attenuation_values.tbl'

    tab1 = Table([ID, EBV], names=('ID', 'E(B-V)'))
    asc.write(tab1, out_ascii, format='fixed_width_two_line')

def compute_A(EBV):

    k_arr  = np.array(list(k_dict.values()))

    A_arr  = k_arr * EBV
    A_dict = dict(zip(list(k_dict.keys()),A_arr))

    return A_dict

def line_ratio_atten(ratio, EBV, wave_top, wave_bottom):

    k_top    = k_dict[wave_top]
    k_bottom = k_dict[wave_bottom]

    ratio_atten = ratio * 10**(0.4*EBV*(k_top - k_bottom))

    return ratio_atten
