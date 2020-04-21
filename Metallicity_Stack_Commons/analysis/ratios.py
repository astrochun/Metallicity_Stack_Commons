import numpy as np

from .. import line_name_short, OIII_r
from .temp_metallicity_calc import R_calculation


def flux_ratios(flux_dict, EBV=None):
    """
    Purpose:
      Primary code to determine a variety of line ratios based on a dictionary
      containing emission-line fluxes

    :param flux_dict: dictionary containing line ratios
    :param EBV: array of E(B-V). Same dimensions as contents of flux_dict

    :return flux_ratios_dict: dictionary containing flux ratios
    """

    # Retrieve emission line fluxes
    OII  = flux_dict[line_name_short['OII']]
    OIII = flux_dict[line_name_short['OIII']]
    Hb   = flux_dict[line_name_short['HB']]
    OIII4363 = flux_dict[line_name_short['4363']]

    # Define flux ratios
    two_beta = OII/Hb
    three_beta = (1+1/OIII_r) * OIII/Hb
    logR23 = np.log10(two_beta + three_beta)
    logO32 = np.log10((1+1/OIII_r) * OIII/OII)

    # Define dictionary of flux ratios
    flux_ratios_dict = dict()
    flux_ratios_dict['two_beta'] = two_beta
    flux_ratios_dict['three_beta'] = three_beta
    flux_ratios_dict['logR23'] = logR23
    flux_ratios_dict['logO32'] = logO32

    if EBV is None:
        print("Not applying dust attenuation correction")
        EBV = np.zeros(OIII4363.shape)

    flux_ratios_dict['R'] = R_calculation(OIII4363, OIII, EBV)

    return flux_ratios_dict
