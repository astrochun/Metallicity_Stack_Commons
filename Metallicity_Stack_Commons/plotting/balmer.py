"""
    balmer
    ------

    Generates plots illustrating Balmer recombination lines

    This code was created from:
      https://github.com/astrochun/Zcalbase_gal/blob/master/Analysis/DEEP2_R23_O32/balmer_plots.py
"""

from os.path import join
from typing import Union

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages

from ..analysis.fitting import gauss, double_gauss
from ..analysis.attenuation import compute_EBV

from .. import scalefact, wavelength_dict
from ..column_names import filename_dict
from ..logging import log_stdout, log_verbose

n_rows = 3
n_cols = 3


def extract_fit(astropy_table: Table, line_name: str, balmer: bool = False,
                verbose: bool = False, log: type(log_stdout) = log_stdout()) \
        -> Union[None, dict]:
    """
    Extract best fit from table and fluxes, return a list of
    fitting parameters and fluxes

    :param astropy_table: Astropy table containing fitting result
    :param line_name: Name of Line to extract fit results
    :param balmer: Indicate whether line is a Balmer line
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return: Fitting results
    """

    log_verbose(log, "starting ...", verbose=verbose)

    try:
        astropy_table[line_name + '_Center']
    except KeyError:
        log.warning("Line not present in table")
        log.warning("Exiting!!!")
        return

    xbar = astropy_table[line_name + '_Center']
    sp   = astropy_table[line_name + '_Sigma']
    ap   = astropy_table[line_name + '_Norm']
    con  = astropy_table[line_name + '_Median']

    result_dict = dict()

    result_dict['line_fit'] = [xbar, sp, ap, con]

    result_dict['flux_gauss'] = astropy_table[line_name + '_Flux_Gaussian']
    result_dict['flux_spec']  = astropy_table[line_name + '_Flux_Observed']

    if balmer:
        sn = astropy_table[line_name + '_Abs_Sigma']
        an = astropy_table[line_name + '_Abs_Norm']

        result_dict['line_fit'] += [sn, an]

        result_dict['line_fit_neg'] = [xbar, sn, an, con]

    log_verbose(log, "finished.", verbose=verbose)
    return result_dict


def fitting_result(wave: np.ndarray, y_norm: np.ndarray, lambda_cen: float,
                   line_fit: list, line_fit_neg: list, flux_gauss: float,
                   flux_spec: float, use_revised: bool = False) -> dict:
    """
    Returns fitting results based on inputs of best fit

    :param wave: Array of rest-frame wavelengths
    :param y_norm: Normalize 1-D spectra in units of 10^-17 erg/s/cm2/AA
    :param lambda_cen: Central wavelength in Angstroms
    :param line_fit: List containing Balmer emission fits
    :param line_fit_neg: List containing the absorption ("stellar") Balmer fit
    :param flux_gauss: Flux from Gaussian model
    :param flux_spec: Flux from spectrum (above median)
    :param use_revised: Indicate whether fluxes have been revised.
                        Default: False

    :return: Dictionary of fitting results
    """

    dx = wave[2] - wave[1]

    fit_dict = dict()

    idx_sig = np.where(np.abs((wave - lambda_cen)) / line_fit[1] <= 2.5)[0]
    fit_dict['gauss'] = double_gauss(wave, *line_fit)
    fit_dict['negative'] = gauss(wave, *line_fit_neg)

    gauss0_diff = fit_dict['gauss'] - fit_dict['negative']
    y_norm_diff = y_norm[idx_sig] - fit_dict['negative'][idx_sig]

    # Residuals
    idx_sig_2 = np.where(np.abs((wave - lambda_cen)) / line_fit[1] <= 3.0)[0]
    fit_dict['residual'] = y_norm[idx_sig_2] - fit_dict['gauss'][idx_sig_2] + line_fit[3]

    fit_dict['idx_sig'] = idx_sig_2

    # Fluxes
    if not use_revised:
        fit_dict['flux_gauss'] = np.sum(gauss0_diff * dx)
        fit_dict['flux_spec']  = np.sum(y_norm_diff * dx)
    else:
        fit_dict['flux_gauss'] = flux_gauss
        fit_dict['flux_spec']  = flux_spec

    return fit_dict


# noinspection PyUnboundLocalVariable
def HbHgHd_fits(fitspath: str, out_pdf_prefix: str ='HbHgHd_fits',
                use_revised: bool = False, verbose: bool = False,
                log: type(log_stdout) = log_stdout()):
    """
    Generate PDF plots that illustrate H-delta, H-gamma, and H-beta line
    profiles and best fit

    :param fitspath: Full path
    :param out_pdf_prefix: Prefix for output PDF file
    :param use_revised: Indicate whether to use regular or revised tables
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object
    """

    log_verbose(log, "starting ...", verbose=verbose)

    comp_spec_file = join(fitspath, filename_dict['comp_spec'])
    log.info(f"Reading: {comp_spec_file}")
    stack2D, header = fits.getdata(comp_spec_file, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

    if not use_revised:
        astropy_table_file = join(fitspath, filename_dict['bin_fit'])
        out_pdf = join(fitspath, out_pdf_prefix+'.pdf')
    else:
        astropy_table_file = join(fitspath, filename_dict['bin_fit_MC'])
        out_pdf = join(fitspath, out_pdf_prefix+'.MC.pdf')

    log.info(f"Reading: {astropy_table_file}")
    astropy_table = asc.read(astropy_table_file)

    ID = astropy_table['bin_ID'].data

    pdf_pages = PdfPages(out_pdf)

    for ii in range(len(ID)):
        if ii % n_rows == 0:
            fig, ax_arr = plt.subplots(nrows=n_rows, ncols=n_cols,
                                       squeeze=False)

        y0 = stack2D[ii]
        y_norm = y0/scalefact

        Hb_dict = extract_fit(astropy_table[ii], 'HBETA', balmer=True,
                              verbose=False, log=log)
        Hg_dict = extract_fit(astropy_table[ii], 'HGAMMA', balmer=True,
                              verbose=False, log=log)
        Hd_dict = extract_fit(astropy_table[ii], 'HDELTA', balmer=True,
                              verbose=False, log=log)

        wave_beta  = wavelength_dict['HBETA']
        wave_gamma = wavelength_dict['HGAMMA']
        wave_delta = wavelength_dict['HDELTA']

        # Beta
        Hb_fit_dict = fitting_result(wave, y_norm, wave_beta, **Hb_dict,
                                     use_revised=use_revised)

        # Gamma
        Hg_fit_dict = fitting_result(wave, y_norm, wave_gamma, **Hg_dict,
                                     use_revised=use_revised)

        # Delta
        Hd_fit_dict = fitting_result(wave, y_norm, wave_delta, **Hd_dict,
                                     use_revised=use_revised)

        # Calculate E(B-V)
        EBV_HgHb = compute_EBV(Hg_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss'],
                               source='HgHb', verbose=verbose, log=log)
        EBV_HdHb = compute_EBV(Hd_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss'],
                               source='HdHb', verbose=verbose, log=log)

        row = ii % n_rows

        # Label in upper left once
        ax_arr[row][0].annotate(f'ID: {ID[ii]}', [0.05, 0.95], va='top', ha='left',
                                xycoords='axes fraction', fontsize='8')

        # The below code could be refactored or simplified
        txt0 = fr"+$\sigma$: {Hb_dict['line_fit'][1]:.3f}, " + \
               fr"-$\sigma$: {Hb_dict['line_fit_neg'][1]:.3f}" + "\n"
        txt0 += f"F_G: {Hb_fit_dict['flux_gauss']:.3f} " + \
                f"F_S: {Hb_fit_dict['flux_spec']:.3f}" + "\n"

        ax_arr[row][2].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][2].plot(wave, Hb_fit_dict['gauss'], 'm', linewidth=0.25,
                            label='Beta Fit')
        ax_arr[row][2].set_xlim(4810, 4910)

        ax_arr[row][2].annotate(txt0, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='5')
        ax_arr[row][2].plot(wave[Hb_fit_dict['idx_sig']], Hb_fit_dict['residual'],
                            'r', linestyle='dashed', linewidth=0.2,
                            label='Residuals')

        txt1 = fr"+$\sigma$: {Hg_dict['line_fit'][1]:.3f}, " + \
               fr"-$\sigma$: {Hg_dict['line_fit_neg'][1]:.3f}" + "\n"
        txt1 += f"F_G: {Hg_fit_dict['flux_gauss']:.3f} " + \
                f"F_S: {Hg_fit_dict['flux_spec']:.3f}" + "\n"
        HgHb = Hg_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss']
        txt1 += fr"H$\gamma$/H$\beta$: {HgHb:.2f} E(B-V): {EBV_HgHb:.2f}"

        ax_arr[row][1].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][1].plot(wave, Hg_fit_dict['gauss'], 'm', linewidth=0.25,
                            label='Gamma Fit')
        ax_arr[row][1].set_xlim(4290, 4390)

        ax_arr[row][1].annotate(txt1, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='5')
        ax_arr[row][1].plot(wave[Hg_fit_dict['idx_sig']], Hg_fit_dict['residual'],
                            'r', linestyle='dashed', linewidth=0.2,
                            label='Residuals')

        txt2 = fr"+$\sigma$: {Hd_dict['line_fit'][1]:.3f}, " + \
               fr"-$\sigma$: {Hd_dict['line_fit_neg'][1]:.3f}" + "\n"
        txt2 += f"F_G: {Hd_fit_dict['flux_gauss']:.3f} " + \
                f"F_S: {Hd_fit_dict['flux_spec']:.3f}" + "\n"
        HdHb = Hd_fit_dict['flux_gauss']/Hb_fit_dict['flux_gauss']
        txt2 += fr"H$\delta/H$\beta$: {HdHb:.2f} E(B-V): {EBV_HdHb:.2f}"

        ax_arr[row][0].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][0].plot(wave, Hd_fit_dict['gauss'], 'm', linewidth=0.25,
                            label='Delta Fit')
        ax_arr[row][0].set_xlim(4050, 4150)

        ax_arr[row][0].set_ylim(0, 1.6)
        
        ax_arr[row][0].annotate(txt2, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='5')
        ax_arr[row][0].plot(wave[Hd_fit_dict['idx_sig']], Hd_fit_dict['residual'],
                            'r', linestyle='dashed', linewidth=0.2,
                            label='Residuals')
       
        ax_arr[row][0].set_yticklabels([0, 0.5, 1, 1.5])
        ax_arr[row][1].set_yticklabels([])
        ax_arr[row][2].set_yticklabels([])

        if row != n_rows-1 and ii != stack2D.shape[0]-1:
            ax_arr[row][0].set_xticklabels([])
            ax_arr[row][1].set_xticklabels([])
            ax_arr[row][2].set_xticklabels([])
        else:
            ax_arr[row][0].set_xticklabels([4050, 4075, 4100, 4125])

        for col in range(n_cols):
            ax_arr[row][col].tick_params(direction='in')  # ticks on the inside

        if row == 1:
            ax_arr[row][0].set_ylabel(r"Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]",
                                      fontsize=12)

        if row == n_rows-1 or ii == stack2D.shape[0]-1:
            ax_arr[row][1].set_xlabel(r"Wavelength [$\AA$]", fontsize=12)

        if ii % n_rows == n_rows-1:
            plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99,
                                hspace=0.05, wspace=0.025)
            fig.savefig(pdf_pages, format='pdf')

    log.info(f"Writing : {out_pdf}")
    pdf_pages.close()

    log_verbose(log, "finished.", verbose=verbose)
