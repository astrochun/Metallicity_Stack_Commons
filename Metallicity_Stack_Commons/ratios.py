from . import line_name

def error_prop_flux_ratios(flux_dict):
    """
    Purpose:
      Primary code to determine a variety of line ratios based on a dictionary
      containing emission-line fluxes

    :param flux_dict: dictionary containing line ratios
    :return:
    """

    two_beta = flux_dict['OII_3727']