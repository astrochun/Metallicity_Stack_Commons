"""
    balmer
    ------

    Generates plots illustrating Balmer recombination lines

    This code was created from:
      https://github.com/astrochun/Zcalbase_gal/blob/master/Analysis/DEEP2_R23_O32/balmer_plots.py
"""

from os.path import join
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages

from ..analysis.fitting import gauss, double_gauss
from ..analysis.attenuation import compute_EBV

from .. import scalefact, wavelength_dict
from ..column_names import filename_dict

n_rows = 3
n_cols = 3


def extract_fit(astropy_table, line_name, balmer=False):
    """
    Purpose:
      Extract best fit from table, return a list of fitting parameters

    :param astropy_table: Astropy table containing fitting result
    :param line_name: line to extract fit results
    :param balmer: boolean to indicate whether line is a Balmer line

    :return:
    param_list: list containing fit to pass in
    """

    try:
        astropy_table[line_name + '_Center']
    except KeyError:
        print("Line not present in table")
        print("Exiting!!!")
        return

    xbar = astropy_table[line_name + '_Center']
    sp   = astropy_table[line_name + '_Sigma']
    ap   = astropy_table[line_name + '_Norm']
    con  = astropy_table[line_name + '_Median']

    param_list = [xbar, sp, ap, con]

    if balmer:
        sn = astropy_table[line_name + '_Abs_Sigma']
        an = astropy_table[line_name + '_Abs_Norm']

        param_list += [sn, an]

        param_list_neg = [xbar, sn, an, con]

        return param_list, param_list_neg
    else:
        return param_list


def fitting_result(wave, y_norm, lambda_cen, balmer_fit, balmer_fit_neg):
    """
    Purpose:
      Returns fitting results based on inputs of best fit

    :param wave: numpy array of rest-frame wavelength
    :param y_norm: Normalize 1-D spectra in units of 10^-17 erg/s/cm2/AA
    :param lambda_cen: Central wavelength in Angstroms
    :param balmer_fit: list containing Balmer emission fits
    :param balmer_fit_neg: list containing the absorption ("stellar") Balmer fit
    :return:
    """

    dx = wave[2] - wave[1]

    x_sigsnip   = np.where(np.abs((wave - lambda_cen)) / balmer_fit[1] <= 2.5)[0]
    gauss0      = double_gauss(wave, *balmer_fit)
    neg0        = gauss(wave, *balmer_fit_neg)
    gauss0_diff = gauss0 - neg0
    y_norm_diff = y_norm[x_sigsnip] - neg0[x_sigsnip]

    # Residuals
    x_sigsnip_2 = np.where(np.abs((wave - lambda_cen)) / balmer_fit[1] <= 3.0)[0]
    resid       = y_norm[x_sigsnip_2] - gauss0[x_sigsnip_2] + balmer_fit[3]

    # Fluxes
    flux_g = np.sum(gauss0_diff * dx)
    flux_s = np.sum(y_norm_diff * dx)

    return gauss0, resid, x_sigsnip_2, flux_g, flux_s


# noinspection PyUnboundLocalVariable
def HbHgHd_fits(fitspath, out_pdf_prefix='HbHgHd_fits',
                use_revised=False):
    """
    Purpose:
      Generate PDF plots that illustrate H-delta, H-gamma, and H-beta line
      profiles and best fit

    :param fitspath: full path (str)
    :param out_pdf_prefix: Prefix for outpute PDF file (str)
    :param use_revised: Indicate whether to use regular or revised tables (bool)
    """

    comp_spec_file = filename_dict['comp_spec']
    stack2D, header = fits.getdata(comp_spec_file, header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])

    if not use_revised:
        astropy_table_file = join(fitspath, filename_dict['bin_fit'])
        out_pdf = join(fitspath, out_pdf_prefix+'.pdf')
    else:
        astropy_table_file = join(fitspath, filename_dict['bin_fit_rev'])
        out_pdf = join(fitspath, out_pdf_prefix+'.revised.pdf')

    astropy_table = asc.read(astropy_table_file)

    ID = astropy_table['bin_ID'].data

    pdf_pages = PdfPages(out_pdf)

    for ii in range(len(ID)):
   
        if ii % n_rows == 0:
            fig, ax_arr = plt.subplots(nrows=n_rows, ncols=n_cols, squeeze=False)

        y0 = stack2D[ii]
        y_norm = y0/scalefact

        Hb_fit, Hb_fit_neg = extract_fit(astropy_table[ii], 'HBETA', balmer=True)
        Hg_fit, Hg_fit_neg = extract_fit(astropy_table[ii], 'HGAMMA', balmer=True)
        Hd_fit, Hd_fit_neg = extract_fit(astropy_table[ii], 'HDELTA', balmer=True)

        wave_beta  = wavelength_dict['HBETA']
        wave_gamma = wavelength_dict['HGAMMA']
        wave_delta = wavelength_dict['HDELTA']

        # Beta
        fit_result_in = [wave, y_norm, wave_beta, Hb_fit, Hb_fit_neg]
        Bgauss0, Bresid, Bx_sigsnip_2, Bflux_g, Bflux_s = fitting_result(*fit_result_in)

        # Gamma
        fit_result_in = [wave, y_norm, wave_gamma, Hg_fit, Hg_fit_neg]
        Ggauss0, Gresid, Gx_sigsnip_2, Gflux_g, Gflux_s = fitting_result(*fit_result_in)

        # Delta
        fit_result_in = [wave, y_norm, wave_delta, Hd_fit, Hd_fit_neg]
        Dgauss0, Dresid, Dx_sigsnip_2, Dflux_g, Dflux_s = fitting_result(*fit_result_in)

        # Calculate E(B-V)
        EBV_HgHb = compute_EBV(Gflux_g/Bflux_g, source='HgHb')
        EBV_HdHb = compute_EBV(Dflux_g/Bflux_g, source='HdHb')

        row = ii % n_rows

        # Label in upper left once
        ax_arr[row][0].annotate(f'ID: {ID[ii]}', [0.05, 0.95], va='top', ha='left',
                                xycoords='axes fraction', fontsize='8')

        # The below code could be refactored or simplified
        txt0 = r'+$\sigma$: %.3f, -$\sigma$: %.3f  ' % (Hb_fit[1], Hb_fit_neg[1]) + '\n'
        txt0 += 'F_G: %.3f F_S: %.3f' % (Bflux_g, Bflux_s)

        ax_arr[row][2].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][2].plot(wave, Bgauss0, 'm', linewidth=0.25, label='Beta Fit')
        ax_arr[row][2].set_xlim(wave_beta-50, wave_beta+50)

        ax_arr[row][2].annotate(txt0, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='5')
        ax_arr[row][2].plot(wave[Bx_sigsnip_2], Bresid, 'r', linestyle='dashed',
                            linewidth=0.2, label='Residuals')

        txt1 = r'+$\sigma$: %.3f, -$\sigma$: %.3f  ' % (Hg_fit[1], Hg_fit_neg[1]) + '\n'
        txt1 += 'F_G: %.3f F_S: %.3f' % (Gflux_g, Gflux_s) + '\n'
        txt1 += r'H$\gamma$/H$\beta$: %.2f E(B-V): %.2f' % (Gflux_g/Bflux_g, EBV_HgHb)

        ax_arr[row][1].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][1].plot(wave, Ggauss0, 'm', linewidth=0.25, label='Gamma Fit')
        ax_arr[row][1].set_xlim(wave_gamma-50, wave_gamma+50)

        ax_arr[row][1].annotate(txt1, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='5')
        ax_arr[row][1].plot(wave[Gx_sigsnip_2], Gresid, 'r', linestyle='dashed',
                            linewidth=0.2, label='Residuals')

        txt2 = r'+$\sigma$: %.3f, -$\sigma$: %.3f  ' % (Hd_fit[1], Hd_fit_neg[1]) + '\n'
        txt2 += 'F_G: %.3f F_S: %.3f' % (Dflux_g, Dflux_s) + '\n'
        txt2 += r'H$\delta$/H$\beta$: %.2f E(B-V): %.2f' % (Dflux_g/Bflux_g, EBV_HdHb)

        ax_arr[row][0].plot(wave, y_norm, 'k', linewidth=0.3, label='Emission')
        ax_arr[row][0].plot(wave, Dgauss0, 'm', linewidth=0.25, label='Delta Fit')
        ax_arr[row][0].set_xlim(wave_delta-50, wave_delta+50)

        ax_arr[row][0].set_ylim(0, 1.5)
        
        ax_arr[row][0].annotate(txt2, [0.95, 0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize='5')
        ax_arr[row][0].plot(wave[Dx_sigsnip_2], Dresid, 'r', linestyle='dashed',
                            linewidth=0.2, label='Residuals')
       
        ax_arr[row][0].set_yticklabels([0, 0.5, 1, 1.5])
        ax_arr[row][1].set_yticklabels([])
        ax_arr[row][2].set_yticklabels([])

        if row != n_rows-1 and ii != stack2D.shape[0]-1:
            ax_arr[row][0].set_xticklabels([])
            ax_arr[row][1].set_xticklabels([])
            ax_arr[row][2].set_xticklabels([])

        if ii % n_rows == n_rows-1:
            plt.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95)
            fig.savefig(pdf_pages, format='pdf')

    pdf_pages.close()
