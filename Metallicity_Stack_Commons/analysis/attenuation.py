import numpy as np

from .. import k_dict, line_name_short
from ..logging import log_stdout

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


def compute_EBV(ratio, source='HgHb', zero_neg=True, log=None):
    """
    Purpose:
      Determines E(B-V) from Hg/Hb or Hd/Hb flux ratios using Case B assumptions

    :param ratio: float or numpy array containing Hg/Hb or Hd/Hb
    :param source: str indicate ratio type.  Either 'HgHb' or 'HdHb'. Default: 'HgHb'
    :param zero_neg: boolean to indicate whether to zero out negative reddening. Default: True
    :param log: LogClass object

    :return EBV: float or numpy array containing E(B-V).
                 Note: Not correcting for negative reddening
    :return EBV_peak: float or numpy array return when it is a 2-D distribution
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")
    if isinstance(ratio, list):
        log.warning("!!! Incorrect type for input [ratio].  Cannot be list !!!")
        raise TypeError

    if 'HgHb' in source:
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
                log.debug("finished ...")
                return EBV
            if len(EBV.shape) == 2:
                EBV_avg = np.average(EBV, axis=0)  # initial guess
                EBV_err, EBV_peak = compute_onesig_pdf(EBV, EBV_avg, usepeak=True)
                neg_idx = np.where(EBV_peak < 0.0)[0]
                if len(neg_idx) > 0:
                    log.info('EBV distribution shifted for peak')
                    EBV[neg_idx, :] -= EBV_peak[neg_idx].reshape((len(neg_idx), 1))
                    EBV_peak[neg_idx] = 0.0
                log.debug("finished ...")
                return EBV, EBV_peak
    else:
        log.debug("finished ...")
        return EBV


def compute_A(EBV):
    """
    Purpose:
      Compute A(Lambda) for all possible emission lines

    :param EBV: float value of E(B-V)
      Has not been configured to handle a large array.  Some array handling would be needed

    :return A_dict: dict containing A(lambda) with keys identical to k_dict
    """

    k_arr  = np.array(list(k_dict.values()))

    A_arr  = k_arr * EBV
    A_dict = dict(zip(list(k_dict.keys()), A_arr))

    return A_dict


def line_ratio_atten(ratio, EBV, wave_top, wave_bottom):

    k_top    = k_dict[wave_top]
    k_bottom = k_dict[wave_bottom]

    ratio_atten = ratio * 10**(0.4*EBV*(k_top - k_bottom))

    return ratio_atten


def Hb_SFR(log_LHb, EBV):
    """
    Purpose:
      Determine dust-corrected SFR using the H-beta luminosity and a
      measurement for nebular attenuation

    Equation below is based on Eq. 2 in Ly et al. (2015), ApJ, 805, 45
      DOI: https://doi.org/10.1088/0004-637X/805/1/45

    :param log_LHb: numpy array or float containing logarithm of H-beta
           luminosity in units of erg/s
    :param EBV: numpy array or float providing E(B-V)

    :return logSFR: numpy array or float containing the SFR in
            logarithmic units of M_sun/yr
    """

    logSFR = np.log10(4.4e-42 * HaHb_CaseB) + 0.4*EBV*k_HBETA + log_LHb

    return logSFR
