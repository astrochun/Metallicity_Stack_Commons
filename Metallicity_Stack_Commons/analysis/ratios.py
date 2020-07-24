import numpy as np

from .. import line_name_short, OIII_r
from .temp_metallicity_calc import R_calculation
from ..column_names import bin_ratios0, indv_names0


def flux_ratios(flux_dict, flux_type='individual', EBV=None, get_R=True):
    """
    Purpose:
      Primary code to determine a variety of line ratios based on a dictionary
      containing emission-line fluxes

    :param flux_dict: dictionary containing line ratios
    :param flux_type: str describing input type. Either 'composite' or 'individual'
    :param EBV: array of E(B-V). Same dimensions as contents of flux_dict
    :param get_R: Boolean to indicate whether to get OIII4363/OIII5007 flux ratio for temp calculation

    :return flux_ratios_dict: dictionary containing flux ratios
    """

    if flux_type not in ['composite', 'individual']:
        print("Incorrect [flux_type]")
        raise ValueError

    two_beta_key = indv_names0[5]
    three_beta_key = indv_names0[6]
    logR23_key = indv_names0[1]
    logO32_key = indv_names0[2]

    if flux_type == 'composite':
        two_beta_key = bin_ratios0[2]
        three_beta_key = bin_ratios0[3]
        logR23_key = bin_ratios0[0]
        logO32_key = bin_ratios0[1]

    # Retrieve emission line fluxes
    OII  = flux_dict[line_name_short['OII']]
    OIII = flux_dict[line_name_short['OIII']]
    Hb   = flux_dict[line_name_short['HB']]
    if get_R:
        OIII4363 = flux_dict[line_name_short['4363']]

    # Define flux ratios
    two_beta = OII/Hb
    three_beta = (1+1/OIII_r) * OIII/Hb
    logR23 = np.log10(two_beta + three_beta)
    logO32 = np.log10((1+1/OIII_r) * OIII/OII)

    # Define dictionary of flux ratios
    flux_ratios_dict = dict()

    flux_ratios_dict[two_beta_key] = two_beta
    flux_ratios_dict[three_beta_key] = three_beta
    flux_ratios_dict[logR23_key] = logR23
    flux_ratios_dict[logO32_key] = logO32

    if EBV is None:
        print("Not applying dust attenuation correction")
        EBV = np.zeros(OIII.shape)

    if get_R:
        flux_ratios_dict['R'] = R_calculation(OIII4363, OIII, EBV)

    return flux_ratios_dict
