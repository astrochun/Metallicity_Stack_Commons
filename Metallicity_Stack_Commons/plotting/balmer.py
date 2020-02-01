'''
    balmer
    ------

    Generates plots illustrating Balmer recombination lines

    This code was created from:
      https://github.com/astrochun/Zcalbase_gal/blob/master/Analysis/DEEP2_R23_O32/balmer_plots.py
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages

from ..fitting import gauss, double_gauss

from .. import scalefact, wavelength_dict

nrows = 3
ncols = 3

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
        xbar = astropy_table[line_name + '_X_bar'].data
    except KeyError:
        print("Line not present in table")
        print("Exiting!!!")
        return

    xbar = astropy_table[line_name + '_X_bar'].data
    sp   = astropy_table[line_name + '_Pos_Sig'].data
    ap   = astropy_table[line_name + '_Pos_Amp'].data
    con  = astropy_table[line_name + '_Const'].data

    param_list = [xbar, sp, ap, con]

    if balmer:
        sn = astropy_table[line_name + '_Neg_Sig'].data
        an = astropy_table[line_name + '_Neg_Amp'].data

        param_list += [sn, an]

        param_list_neg = [xbar, sn, an, con]

        return param_list, param_list_neg
    else:
        return param_list


def fitting_result(wave, y_norm, lambda_cen, balmer_fit, balmer_fit_neg):
    x_sigsnip   = np.where(np.abs((wave - lambda_cen)) / balmer_fit[1] <= 2.5)[0]
    gauss0      = double_gauss(wave, *balmer_fit)
    neg0        = gauss(wave, *balmer_fit_neg)
    gauss0_diff = gauss0 - neg0
    y_norm_diff = y_norm[x_sigsnip] - Bneg0[x_sigsnip]

    # Residuals
    x_sigsnip_2 = np.where(np.abs((wave - lambda_cen)) / balmer_fit[1] <= 3.0)[0]
    resid       = y_norm[x_sigsnip_2] - gauss0[x_sigsnip_2] + balmer_fit[3]

    # Fluxes
    flux_g = np.sum(gauss0_diff * dx)
    flux_s = np.sum(y_norm_diff * dx)

    return gauss0, resid, x_sigsnip_2, flux_g, flux_s


def HbHgHd_fits(Stack_name, astropy_table_file, out_pdf):

    stack2D, header = fits.getdata(Stack_name,header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    dispersion = header['CDELT1']

    astropy_table = asc.read(astropy_table_file)

    ID = astropy_table['ID'].data

    pdfpages = PdfPages(out_pdf)

    for ii in range(len(ID)):
   
        if ii % nrows ==0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)

        y0 = stack2D[ii]
        y_norm = y0/scalefact
        dx = wave[2]-wave[1]

        Hb_fit, Hb_fit_neg = extract_fit(astropy_table[ii], 'HBETA', balmer=True)
        Hg_fit, Hg_fit_neg = extract_fit(astropy_table[ii], 'Hgamma', balmer=True)
        Hd_fit, Hd_fit_neg = extract_fit(astropy_table[ii], 'HDELTA', balmer=True)

        wave_beta  = wavelength_dict['HBETA']
        wave_gamma = wavelength_dict['HGAMMA']
        wave_delta = wavelength_dict['HDELTA']

        ##Beta
        fit_result_in = [wave, y_norm, wave_beta, Hb_fit, Hb_fit_neg]
        Bgauss0, Bresid, Bx_sigsnip_2, Bflux_g, Bflux_s = fitting_result(fit_result_in)

        ##Gamma
        fit_result_in = [wave, y_norm, wave_gamma, Hg_fit, Hg_fit_neg]
        Ggauss0, Gresid, Gx_sigsnip_2, Gflux_g, Gflux_s = fitting_result(fit_result_in)

        ##Delta
        fit_result_in = [wave, y_norm, wave_delta, Hd_fit, Hd_fit_neg]
        Dgauss0, Dresid, Dx_sigsnip_2, Dflux_g, Dflux_s = fitting_result(fit_result_in)

        row = ii % nrows

        txt0 = r'ID: %i' % (ID[ii]) +'\n'
        txt0 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  '% (Hb_fit[1], Hb_fit_neg[1]) + '\n'
        txt0 += 'F_G: %.3f F_S: %.3f' %(Bflux_g, Bflux_s)

        ax_arr[row][2].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
        ax_arr[row][2].plot(wave,Bgauss0, 'm', linewidth= 0.25, label= 'Beta Fit')
        ax_arr[row][2].set_xlim(wave_beta-50, wave_beta+50)

        ax_arr[row][2].annotate(txt0, [0.95,0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize= '5')
        ax_arr[row][2].plot(wave[Bx_sigsnip_2],Bresid, 'r', linestyle='dashed',
                            linewidth = 0.2, label= 'Residuals')

        txt1 = r'ID: %i' % (ID[ii]) +'\n'
        txt1 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  '% (Hg_fit[1], Hg_fit_neg[1]) + '\n'
        txt1 += 'F_G: %.3f F_S: %.3f' %(Gflux_g, Gflux_s)
    
        ax_arr[row][1].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
        ax_arr[row][1].plot(wave,Ggauss0, 'm', linewidth= 0.25, label= 'Gamma Fit')
        ax_arr[row][1].set_xlim(wave_gamma-50, wave_gamma+50)

        ax_arr[row][1].annotate(txt1, [0.95,0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize= '5')
        ax_arr[row][1].plot(wave[Gx_sigsnip_2],Gresid, 'r', linestyle='dashed',
                            linewidth = 0.2, label= 'Residuals')

        txt2 = r'ID: %i' % (ID[ii]) +'\n'
        txt2 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  '% (Hd_fit[1], Hd_fit_neg[1]) + '\n'
        txt2 += 'F_G: %.3f F_S: %.3f' %(Dflux_g, Dflux_s)

        ax_arr[row][0].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
        ax_arr[row][0].plot(wave,Dgauss0, 'm', linewidth= 0.25, label= 'Detla Fit')
        ax_arr[row][0].set_xlim(wave_delta-50, wave_delta+50)

        ax_arr[row][0].set_ylim(0,1.5)
        
        ax_arr[row][0].annotate(txt0, [0.95,0.95], xycoords='axes fraction',
                                va='top', ha='right', fontsize= '5')
        ax_arr[row][0].plot(wave[Dx_sigsnip_2],Dresid, 'r', linestyle='dashed',
                            linewidth = 0.2, label= 'Residuals')

       
        ax_arr[row][0].set_yticklabels([0,0.5,1,1.5])
        ax_arr[row][1].set_yticklabels([])
        ax_arr[row][2].set_yticklabels([])

        if row != nrows-1 and ii != stack2D.shape[0]-1:
            ax_arr[row][0].set_xticklabels([])
            ax_arr[row][1].set_xticklabels([])
            ax_arr[row][2].set_xticklabels([])

        if ii % nrows == nrows-1: fig.savefig(pdfpages, format='pdf')

    pdfpages.close()
