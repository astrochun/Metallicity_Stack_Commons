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

        EBV = np.log10((HBETA / HGAMMA) * (HgHb_CaseB)) * 2.5 * (1 / (k_HGAMMA - k_HBETA))

        out_ascii = fitspath + '/dust_attenuation_values.tbl'
        # if not exists(out_ascii_single):
        n2 = ('ID', 'E(B-V)')
        tab1 = Table([ID, EBV], names=n2)
        asc.write(tab1, out_ascii, format='fixed_width_two_line')