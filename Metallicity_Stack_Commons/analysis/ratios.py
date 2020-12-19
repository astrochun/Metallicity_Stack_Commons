import numpy as np

from .. import line_name_short, OIII_r
from .temp_metallicity_calc import R_calculation
from ..column_names import bin_ratios0, dust0
from ..logging import log_stdout


def flux_ratios(flux_dict, binned_data=False, get_R=True, log=None):
    """
    Purpose:
      Primary code to determine a variety of line ratios based on a dictionary
      containing emission-line fluxes

    :param flux_dict: dict containing line ratios
    :param get_R: bool to indicate to populate with OIII4363/OIII5007 flux ratio
    :param binned_data: bool for whether to analysis binned data. Default: False
    :param log: LogClass or logging object

    :return flux_ratios_dict: dict containing flux ratios
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    if not binned_data:
        bin_ratios = [ratios0.replace('_composite', '')
                      for ratios0 in bin_ratios0]
        dust = [t_dust.replace('_composite', '') for t_dust in dust0]
    else:
        bin_ratios = bin_ratios0
        dust = dust0

    two_beta_key = bin_ratios[2]
    three_beta_key = bin_ratios[3]
    logR23_key = bin_ratios[0]
    logO32_key = bin_ratios[1]
    R_key = bin_ratios[4]

    # Retrieve emission-line fluxes
    OII  = flux_dict[line_name_short['OII']]
    OIII = flux_dict[line_name_short['OIII']]
    Hb   = flux_dict[line_name_short['HB']]

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

    # Add Balmer decrement
    if line_name_short['HG'] in flux_dict:
        flux_ratios_dict[dust[0]] = flux_dict[line_name_short['HG']] / Hb
    if line_name_short['HD'] in flux_dict:
        flux_ratios_dict[dust[1]] = flux_dict[line_name_short['HD']] / Hb

    if get_R:
        OIII4363 = flux_dict[line_name_short['4363']]  # Retrieve OIII4363
        flux_ratios_dict[R_key] = R_calculation(OIII4363, OIII)

    log.debug("finished ...")

    return flux_ratios_dict
