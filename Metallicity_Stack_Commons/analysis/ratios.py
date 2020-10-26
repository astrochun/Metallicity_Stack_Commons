import numpy as np

from .. import line_name_short, OIII_r
from .temp_metallicity_calc import R_calculation
from ..column_names import bin_ratios0


def flux_ratios(flux_dict, binned_data=False, get_R=True):
    """
    Purpose:
      Primary code to determine a variety of line ratios based on a dictionary
      containing emission-line fluxes

    :param flux_dict: dictionary containing line ratios
    :param get_R: Boolean to indicate whether to get OIII4363/OIII5007 flux ratio for temp calculation
    :param binned_data: bool for whether to analysis binned data. Default: False
    :return flux_ratios_dict: dictionary containing flux ratios
    """

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

    if not binned_data:
        bin_ratios = [ratios0.replace('_composite', '') for ratios0 in bin_ratios0]
    else:
        bin_ratios = bin_ratios0

    flux_ratios_dict[bin_ratios[2]] = two_beta
    flux_ratios_dict[bin_ratios[3]] = three_beta
    flux_ratios_dict[bin_ratios[0]] = logR23
    flux_ratios_dict[bin_ratios[1]] = logO32

    # Add Balmer decrement
    if line_name_short['HG'] in flux_dict:
        flux_ratios_dict['HgHb'] = flux_dict[line_name_short['HG']] / Hb
    if line_name_short['HD'] in flux_dict:
        flux_ratios_dict['HdHb'] = flux_dict[line_name_short['HD']] / Hb

    if get_R:
        flux_ratios_dict[bin_ratios[4]] = R_calculation(OIII4363, OIII)

    return flux_ratios_dict
