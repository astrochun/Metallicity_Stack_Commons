import numpy as np

from . import k_dict

k_4363 = k_dict['OIII_4363']
k_5007 = k_dict['OIII_5007']

# Constants

a = 13205
b = 0.92506
c = 0.98062


def R_calculation(OIII4363, OIII5007, EBV):
    """
    Computes the excitation flux ratio of [OIII]4363 to [OIII]5007.
    Adopts a 3.1-to-1 ratio for 5007/4959

    :param OIII4363: numpy array of OIII4363 fluxes
    :param OIII5007: numpy array of OIII5007 fluxes
    :param EBV: numpy array of E(B-V).  Set to zero if not applying attenuation

    :return R_value: O++ excitation flux ratio
    """

    R_value = OIII4363 / (OIII5007 * (1 + 1 / 3.1)) * 10 ** (0.4 * EBV * (k_4363 - k_5007))

    return R_value


def temp_calculation(R):
    """
    Computes electron temperature (T_e) from O++ excitation flux ratio

    Formula is:
        T_e = a(-log(R)-b)^(-c)
    where a = 13025, b=0.92506, and c=0.98062 (Nicholls et al. 2014)

    :param R: numpy array of O++ excitation flux ratio (see R_calculation)

    :return T_e: numpy array of T_e (Kelvins)
    """

    T_e = a * (-np.log10(R) - b) ** (-1 * c)

    return T_e


def metallicity_calculation(T_e, TWO_BETA, THREE_BETA):
    """
    Determines 12+log(O/H) from electron temperature and [OII]/Hb and [OIII]/Hb flux ratio

    :param T_e: numpy array of temperature (see temp_calculation)
    :param TWO_BETA: numpy array of [OII]/Hb flux ratio
    :param THREE_BETA: numpy array of [OIII]/Hb flux ratio

    :return com_O_log: numpy array of 12+log(O/H)
    :return metal_dict: dictionary containing O+/H, O++/H, log(O+/H), log(O++/H)
    """

    t_3 = T_e * 1e-4
    t_2 = 0.7 * t_3 + 0.17
    x2 = 1e-4 * 1e3 * t_2 ** (-0.5)

    # Equations from Izotov et al. (2006)
    O_s_ion_log = np.log10(TWO_BETA) + 5.961 + 1.676 / t_2 - 0.4 * np.log10(t_2) \
                  - 0.034 * t_2 + np.log10(1 + 1.35 * x2) - 12
    O_d_ion_log = np.log10(THREE_BETA) + 6.200 + 1.251 / t_3 \
                  - 0.55 * np.log10(t_3) - 0.014 * (t_3) - 12

    O_s_ion = 10 ** (O_s_ion_log)
    O_d_ion = 10 ** (O_d_ion_log)
    com_O = O_s_ion + O_d_ion
    com_O_log = np.log10(com_O) + 12

    metal_dict = dict(O_s_ion=O_s_ion, O_d_ion=O_d_ion,
                      O_s_ion_log=O_s_ion_log, O_d_ion_log=O_d_ion_log)

    return com_O_log, metal_dict
