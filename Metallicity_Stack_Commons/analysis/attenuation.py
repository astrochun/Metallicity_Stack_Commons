from typing import Union

import numpy as np

from .. import k_dict, line_name_short
from ..logging import log_stdout, log_verbose

from chun_codes import compute_onesig_pdf

# Balmer decrement Case B, zero reddening
HdHb_CaseB = 0.259  # Hd/Hb ratio
HgHb_CaseB = 0.468  # Hg/Hb ratio
HaHb_CaseB = 2.86   # Ha/Hb ratio

HB = line_name_short['HB']
HG = line_name_short['HG']
HD = line_name_short['HD']

k_HBETA  = k_dict[HB]
k_HGAMMA = k_dict[HG]
k_HDELTA = k_dict[HD]


def compute_EBV(ratio: Union[float, np.ndarray], source: str = 'HgHb',
                zero_neg: bool = True, verbose: bool = False,
                log: type(log_stdout) = log_stdout()):
    """
    Determines E(B-V) from Hg/Hb or Hd/Hb flux ratios using Case B assumptions

    :param ratio: Float or array containing Hg/Hb or Hd/Hb values
    :param source: Indicate ratio type. Either 'HgHb' or 'HdHb'.
                   Default: 'HgHb'
    :param zero_neg: Indicate whether to zero out negative reddening.
                     Default: True
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return EBV: E(B-V) values
            Note: Not correcting for negative reddening
    :rtype: float|np.ndarray
    :return EBV_peak: E(B-V) peak values
    :rtype: float|np.ndarray
    """

    log_verbose(log, "starting ...", verbose=verbose)

    if isinstance(ratio, list):
        log.warning("!!! Incorrect type for input [ratio].  Cannot be list !!!")
        raise TypeError("!!! Incorrect type for input [ratio].  Cannot be list !!!")

    if source not in ['HgHb', 'HdHb']:
        log.warning("!!! Incorrect [source]")
        raise KeyError("!!! Incorrect [source]")

    # For HgHb (default)
    ratio0 = HgHb_CaseB
    k1 = k_HGAMMA

    if 'HdHb' in source:
        ratio0 = HdHb_CaseB
        k1 = k_HDELTA

    EBV = -2.5 * np.log10(ratio/ratio0)/(k1 - k_HBETA)

    if zero_neg:
        if isinstance(EBV, float):
            if EBV < 0.0:
                EBV = 0.0
                log.info('zero substituted for negative reddening')
            return EBV
        else:
            if len(EBV.shape) == 1:
                neg_idx = np.where(EBV < 0.0)[0]
                if len(neg_idx) > 0:
                    EBV[neg_idx] = 0.0
                    log.info('zero substituted for negative reddening')
                log_verbose(log, "finished.", verbose=verbose)
                return EBV
            if len(EBV.shape) == 2:
                EBV_avg = np.average(EBV, axis=0)  # initial guess
                EBV_err, EBV_peak = compute_onesig_pdf(EBV, EBV_avg, usepeak=True)
                neg_idx = np.where(EBV_peak < 0.0)[0]
                if len(neg_idx) > 0:
                    log.info('EBV distribution shifted for peak')
                    EBV[neg_idx, :] -= EBV_peak[neg_idx].reshape((len(neg_idx), 1))
                    EBV_peak[neg_idx] = 0.0
                log_verbose(log, "finished.", verbose=verbose)
                return EBV, EBV_peak
    else:
        log_verbose(log, "finished.", verbose=verbose)
        return EBV


def compute_A(EBV: float, verbose: bool = False,
              log: type(log_stdout) = log_stdout()) -> dict:
    """
    Compute A(Lambda) for all possible emission lines

    :param EBV: E(B-V) value
      Has not been configured to handle a large array.
      Some array handling would be needed

    :param verbose: Write verbose message to stdout.
           Default: file only
    :param log: LogClass or logging object

    :return: A(lambda) with keys identical to ``k_dict``
    """

    log_verbose(log, "starting ...", verbose=verbose)

    k_arr = np.array(list(k_dict.values()))

    A_arr = k_arr * EBV
    A_dict = dict(zip(list(k_dict.keys()), A_arr))

    log_verbose(log, "finished.", verbose=verbose)
    return A_dict


def line_ratio_atten(ratio: Union[float, np.ndarray],
                     EBV: Union[float, np.ndarray],
                     wave_top: str, wave_bottom: str,
                     verbose: bool = False,
                     log: type(log_stdout) = log_stdout()) -> \
        Union[float, np.ndarray]:
    """
    Determine dust-corrected emission-line ratios

    :param ratio: Float or array of observed flux ratios
    :param EBV: E(B-V) value(s)
    :param wave_top: Emission-line name for flux ratio numerator
    :param wave_bottom: Emission-line name for flux ratio denominator
    :param verbose: Write verbose message to stdout.
           Default: file only
    :param log: LogClass or logging object

    :return: Float or array of dust-corrected flux ratios
    """

    log_verbose(log, "starting ...", verbose=verbose)

    k_top = k_dict[wave_top]
    k_bottom = k_dict[wave_bottom]

    ratio_atten = ratio * 10**(0.4*EBV*(k_top - k_bottom))

    log_verbose(log, "finished.", verbose=verbose)
    return ratio_atten


def Hb_SFR(log_LHb: Union[float, np.ndarray],
           EBV: Union[float, np.ndarray],
           verbose: bool = False,
           log: type(log_stdout) = log_stdout()) -> \
        Union[float, np.ndarray]:
    """
    Determine dust-corrected SFR using the H-beta luminosity and a
    measurement for nebular attenuation

    Equation below is based on Eq. 2 in Ly et al. (2015), ApJ, 805, 45
      DOI: https://doi.org/10.1088/0004-637X/805/1/45

    :param log_LHb: Logarithm of H-beta luminosity in units of erg/s
    :param EBV: E(B-V) value(s)
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return: SFRs in logarithmic units of M_sun/yr
    """

    log_verbose(log, "starting ...", verbose=verbose)

    logSFR = np.log10(4.4e-42 * HaHb_CaseB) + 0.4*EBV*k_HBETA + log_LHb

    log_verbose(log, "finished.", verbose=verbose)
    return logSFR
