'''
    balmer
    ------

    Generates plots illustrating Balmer recombination lines

    This code was created from:
      https://github.com/astrochun/Zcalbase_gal/blob/master/Analysis/DEEP2_R23_O32/balmer_plots.py
'''


import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii as asc
from matplotlib.backends.backend_pdf import PdfPages

from ..fitting import gauss, double_gauss

from .. import scalefact, wavelength_dict

nrows = 3
ncols = 3

def extract_fit(astropy_table, line_name, balmer=False):
    '''
    Purpose:
      Extract best fit from table and returns a dictionary

    :param astropy_table: Astropy table
    :param line_name: line to extract fit results
    :param balmer: boolean to indicate whether line is a Balmer line

    :return:
    param_list: list containing fit to pass in
    '''

    try:
        xbar = combine_asc[line_name + '_X_bar']
    except KeyError:
        print("Line not present in table")
        print("Exiting!!!")
        return

    xbar = combine_asc[line_name + '_X_bar']
    sp   = combine_asc[line_name + '_Pos_Sig']
    ap   = combine_asc[line_name + '_Pos_Amp']
    con  = combine_asc[line_name + '_Const']

    param_list = [xbar, sp, ap, con]

    if balmer:
        sn = combine_asc[line_name + '_Neg_Sig']
        an = combine_asc[line_name + '_Neg_Amp']

        param_list += [sn, an]

        param_list_neg = [xbar, sn, an, con]

        return param_list, param_list_neg
    else:
        return param_list


def HbHgHd_fits(fitspath, nrow, ncol,Stack_name,combine_flux_tab, out_pdf):

    stack2D, header = fits.getdata(Stack_name,header=True)
    wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
    dispersion = header['CDELT1']

    combine_asc = asc.read(combine_flux_tab)

    ID = combine_asc['ID']

    pdfpages = PdfPages(out_pdf)

    for ii in range(len(ID)):
   
        if ii % nrows ==0:
            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols, squeeze = False)

        y0 = stack2D[ii]
        y_norm = y0/scalefact
        dx = wave[2]-wave[1]

        Hb_fit, Hb_fit_neg = extract_fit(combine_asc[ii], 'HBETA', balmer=True)
        Hg_fit, Hg_fit_neg = extract_fit(combine_asc[ii], 'Hgamma', balmer=True)
        Hd_fit, Hd_fit_neg = extract_fit(combine_asc[ii], 'HDELTA', balmer=True)

        ##Beta
        wave_beta    = wavelength_dict['HBETA']
        Bx_sigsnip   = np.where(np.abs((wave_beta))/Hb_fit[1]<=2.5 )[0]
        Bgauss0      = double_gauss(wave, *Hb_fit)
        Bneg0        = gauss(wave, *Hb_fit_neg)
        Bgauss0_diff = Bgauss0 - Bneg0
        By_norm_diff = y_norm[Bx_sigsnip]-Bneg0[Bx_sigsnip]

        #Residuals
        Bx_sigsnip_2 = np.where(np.abs((wave-wave_beta))/Hb_fit[1]<=3.0 )[0]
        Bresid       = y_norm[Bx_sigsnip_2]-Bgauss0[Bx_sigsnip_2] + Hb_fit[3]

        #Fluxes
        Bflux_g = np.sum(Bgauss0_diff*dx)
        Bflux_s = np.sum(By_norm_diff*dx)
        
        ##Gamma
        wave_gamma   = wavelength_dict['HGAMMA']
        Gx_sigsnip   = np.where(np.abs((wave-wave_gamma))/Hg_fit[1]<=2.5 )[0]
        Ggauss0      = double_gauss(wave, *Hg_fit)
        Gneg0        = gauss(wave, *Hg_fit_neg)
        Ggauss0_diff = Ggauss0 - Gneg0
        Gy_norm_diff = y_norm[Gx_sigsnip]-Gneg0[Gx_sigsnip]

        #Residuals
        Gx_sigsnip_2 = np.where(np.abs((wave-wave_gamma))/Hg_fit[1]<=3.0)
        Gresid       = y_norm[Gx_sigsnip_2]-Ggauss0[Gx_sigsnip_2] + Hg_fit[3]

        #Fluxes
        Gflux_g = np.sum(Ggauss0_diff*dx)
        Gflux_s = np.sum(Gy_norm_diff*dx)

        ##Delta
        wave_delta   = wavelength_dict['HDELTA']
        Dx_sigsnip   = np.where(np.abs((wave-wave_delta))/Hd_fit[1]<=2.5 )[0]
        Dgauss0      = double_gauss(wave, *Hd_fit)
        Dneg0        = gauss(wave, *Hd_fit_neg)
        Dgauss0_diff = Dgauss0 - Dneg0
        Dy_norm_diff = y_norm[Dx_sigsnip]-Dneg0[Dx_sigsnip]

        #Residuals
        Dx_sigsnip_2 = np.where(np.abs((wave-wave_delta))/Hd_fit[1]<=3.0)
        Dresid = y_norm[Dx_sigsnip_2]-Dgauss0[Dx_sigsnip_2] + Hd_fit[3]

        #Fluxes
        Dflux_g = np.sum(Dgauss0_diff*dx)
        Dflux_s = np.sum(Dy_norm_diff*dx)

        row = ii % nrows

        txt0 = r'ID: %i' % (ID[ii]) +'\n'
        txt0 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  '% (Hb_fit[1], Hb_fit_neg[1]) + '\n'
        txt0 += 'F_G: %.3f F_S: %.3f' %(Bflux_g, Bflux_s)
        #txt0 += 'o1[2]: %.3f o1[4]: %.3f  o1[5]: %.3f'% (o1[2], o1[4], o1[5]) + 
    
        ax_arr[row][2].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
        ax_arr[row][2].plot(wave,Bgauss0, 'm', linewidth= 0.25, label= 'Beta Fit')
        ax_arr[row][2].set_xlim(wave_beta-50, wave_beta+50)
        #ax_arr[row][2].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2,fontsize = 3)
        ax_arr[row][2].annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')
        ax_arr[row][2].plot(wave[Bx_sigsnip_2],Bresid, 'r', linestyle= 'dashed', linewidth = 0.2, label= 'Residuals')

        txt1 = r'ID: %i' % (ID[ii]) +'\n'
        txt1 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  '% (Hg_fit[1], Hg_fit_neg[1]) + '\n'
        txt1 += 'F_G: %.3f F_S: %.3f' %(Gflux_g, Gflux_s)
    
        ax_arr[row][1].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
        ax_arr[row][1].plot(wave,Ggauss0, 'm', linewidth= 0.25, label= 'Gamma Fit')
        ax_arr[row][1].set_xlim(wave_gamma-50, wave_gamma+50)
        #ax_arr[row][1].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2,fontsize = 3)
        ax_arr[row][1].annotate(txt1, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')
        ax_arr[row][1].plot(wave[Gx_sigsnip_2],Gresid, 'r', linestyle= 'dashed', linewidth = 0.2, label= 'Residuals')

        txt2 = r'ID: %i' % (ID[ii]) +'\n'
        txt2 += r'+$\sigma$: %.3f, -$\sigma$: %.3f  '% (Hd_fit[1], Hd_fit_neg[1]) + '\n'
        txt2 += 'F_G: %.3f F_S: %.3f' %(Dflux_g, Dflux_s)

        ax_arr[row][0].plot(wave, y_norm,'k', linewidth=0.3, label= 'Emission')
        ax_arr[row][0].plot(wave,Dgauss0, 'm', linewidth= 0.25, label= 'Detla Fit')
        ax_arr[row][0].set_xlim(wave_delta-50, wave_delta+50)

        ax_arr[row][0].set_ylim(0,1.5)
        
        #ax_arr[row][0].legend(bbox_to_anchor=(0.25,0.1), borderaxespad=0, ncol=2,fontsize = 3)
        ax_arr[row][0].annotate(txt0, [0.95,0.95], xycoords='axes fraction', va='top', ha='right', fontsize= '5')
        ax_arr[row][0].plot(wave[Dx_sigsnip_2],Dresid, 'r', linestyle= 'dashed', linewidth = 0.2, label= 'Residuals')

       
        ax_arr[row][0].set_yticklabels([0,0.5,1,1.5])
        ax_arr[row][1].set_yticklabels([])
        ax_arr[row][2].set_yticklabels([])

        if row != nrows-1 and ii != stack2D.shape[0]-1:
            ax_arr[row][0].set_xticklabels([])
            ax_arr[row][1].set_xticklabels([])
            ax_arr[row][2].set_xticklabels([])
        
        
        if ii% nrows == nrows-1: fig.savefig(pdfpages, format='pdf')

    pdfpages.close()
