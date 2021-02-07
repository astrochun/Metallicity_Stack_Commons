from typing import Union
import numpy as np

from .. import k_dict, OIII_r
from ..logging import log_stdout, log_verbose

from ..column_names import temp_metal_names0, remove_from_list

k_3727 = k_dict['OII_3727']
k_4861 = k_dict['HBETA']
k_4363 = k_dict['OIII_4363']
k_5007 = k_dict['OIII_5007']

# Constants
a = 13205
b = 0.92506
c = 0.98062


def R_calculation(OIII4363: Union[float, np.ndarray],
                  OIII5007: Union[float, np.ndarray],
                  verbose: bool = False,
                  log: type(log_stdout) = log_stdout()) \
        -> Union[float, np.ndarray]:
    """
    Computes the excitation flux ratio of [OIII]4363 to [OIII]5007.
    Adopts a 3.1-to-1 ratio for 5007/4959

    :param OIII4363: OIII4363 fluxes
    :param OIII5007: OIII5007 fluxes
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return R_value: O++ excitation flux ratio
    """

    log_verbose(log, "starting ...", verbose=verbose)

    R_value = OIII4363 / (OIII5007 * (1 + 1/OIII_r))

    log_verbose(log, "finished.", verbose=verbose)
    return R_value


def temp_calculation(R: np.ndarray, EBV: Union[None, np.ndarray] = None,
                     verbose: bool = False,
                     log: type(log_stdout) = log_stdout()) -> np.ndarray:
    """
    Computes electron temperature (T_e) from O++ excitation flux ratio

    Formula is:
        T_e = a(-log(R)-b)^(-c)
    where a = 13025, b=0.92506, and c=0.98062 (Nicholls et al. 2014)

    :param R: Array of O++ excitation flux ratio (see R_calculation)
    :param EBV: Array of E(B-V). Set to zero if not applying attenuation
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return T_e: Array of T_e (Kelvins)
    """

    log_verbose(log, "starting ...", verbose=verbose)

    arr_shape = R.shape

    if EBV is None:
        log.info("Not applying dust attenuation correction")
        EBV = np.zeros(arr_shape)

    R_corr = R * 10 ** (0.4 * EBV * (k_4363 - k_5007))

    T_e = a * (-np.log10(R_corr) - b) ** (-1 * c)

    log_verbose(log, "finished.", verbose=verbose)
    return T_e


def metallicity_calculation(T_e: np.ndarray, TWO_BETA: np.ndarray,
                            THREE_BETA: np.ndarray,
                            EBV: Union[None, np.ndarray] = None,
                            det3: Union[None, np.ndarray] = None,
                            verbose: bool = False,
                            log: type(log_stdout) = log_stdout()) -> dict:
    """
    Determines 12+log(O/H) from electron temperature and [OII]/Hb and
    [OIII]/Hb flux ratio

    :param T_e: Array of electron temperatures (see ``temp_calculation``)
    :param TWO_BETA: Array of [OII]/Hb flux ratios
    :param THREE_BETA: Array of [OIII]/Hb flux ratios
    :param EBV: Optional array containing EBV distribution
    :param det3: Optional array to pass in to identify those satisfying det3
                 requirements. Default: None means full array is considered

                 Note: for Monte Carlo inputs, a 1-D np.array index satisfying
                       det3 requirements will suffice
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return metal_dict: Contains 12+log(O/H), O+/H, O++/H, log(O+/H), log(O++/H)
    """

    if log is None:
        log = log_stdout()

    log_verbose(log, "starting ...", verbose=verbose)

    arr_shape = T_e.shape
    t_3 = np.zeros(arr_shape)
    t_2 = np.zeros(arr_shape)
    x2 = np.zeros(arr_shape)

    if EBV is None:
        log.info("Not applying dust attenuation correction")
        EBV = np.zeros(arr_shape)

    if det3 is None:
        det3 = np.arange(arr_shape[0])

    t_3[det3] = T_e[det3] * 1e-4
    t_2[det3] = 0.7 * t_3[det3] + 0.17
    x2[det3]  = 1e-4 * 1e3 * t_2[det3] ** (-0.5)

    O_s_ion_log = np.zeros(arr_shape)
    O_d_ion_log = np.zeros(arr_shape)

    # Equations from Izotov et al. (2006)
    O_s_ion_log[det3] = np.log10(TWO_BETA[det3]) + 5.961 + 1.676 / t_2[det3] \
                        - 0.4 * np.log10(t_2[det3]) - 0.034 * t_2[det3] + \
                        np.log10(1 + 1.35 * x2[det3]) - 12
    O_d_ion_log[det3] = np.log10(THREE_BETA[det3]) + 6.200 + 1.251 / t_3[det3] \
                        - 0.55 * np.log10(t_3[det3]) - 0.014 * (t_3[det3]) - 12

    # Apply dust attenuation to derived properties
    O_s_ion_log[det3] += 0.4 * EBV[det3] * (k_3727 - k_4861)
    O_d_ion_log[det3] += 0.4 * EBV[det3] * (k_5007 - k_4861)

    O_s_ion = 10 ** O_s_ion_log
    O_d_ion = 10 ** O_d_ion_log
    com_O = O_s_ion + O_d_ion
    com_O_log = np.log10(com_O) + 12

    key_dict = remove_from_list(temp_metal_names0, ['T_e'])
    key_values = [com_O_log, O_s_ion_log, O_d_ion_log, O_s_ion, O_d_ion]  # Order matters here
    metal_dict = dict(zip(key_dict, key_values))

    log_verbose(log, "finished.", verbose=verbose)
    return metal_dict
