from os.path import join, exists

from chun_codes import random_pdf, compute_onesig_pdf
from . import line_name
from astropy.io import ascii as asc
import numpy as np

from .column_names import filename_dict
from .ratios import error_prop_flux_ratios
from .temp_metallicity_calc import temp_calculation, metallicity_calculation

"""
def construct_pdf(values, RMS, seed_i=1, n_iter=1000):
    '''
    Constructs probability distribution function (PDF) based on input
    values and their associated uncertainty

    :param values: list or numpy array of values/parameters
    :param RMS: 1-sigma errors associated with values (same dimension)
    :param seed_i: integer value for initial seed for np.random. Default: 1
    :param n_iter: Number of iterations. Default: 1000

    :return pdf_arr: numpy array of size (size of values, n_iter)
    '''

    pdf_arr = random_pdf(values, RMS, seed_i=seed_i, n_iter=n_iter, silent=False)

    return pdf_arr
"""


def error_prop_chuncodes(path, binned_data=True):

    # Define files to read in for binned data
    if binned_data:
        flux_file = join(path, filename_dict['bin_fit'])
        prop_file = join(path, filename_dict['bin_derived_prop'])
        verify_file = join(path, filename_dict['bin_valid'])

    flux_tab0  = asc.read(flux_file)
    prop_tab0  = asc.read(prop_file)

    if binned_data:
        verify_tab = asc.read(verify_file)
        detect = verify_tab['Detection']

        # For now we are only considering those with reliable detection and
        # excluding those with reliable non-detections (detect = 0.5)
        detection = np.where((detect == 1))[0]

        ID = verify_tab['bin_ID'].data
        ID_detect = ID[detection]
        print(ID_detect)

        flux_tab = flux_tab0[detection]
        prop_tab = prop_tab0[detection]

    flux_cols     = [str0+'_Flux_Gaussian' for str0 in line_name]
    flux_rms_cols = [str0+'_RMS' for str0 in line_name]

    Temp = prop_tab['T_e'].data
    com_O_log = prop_tab['12+log(O/H)'].data
    O_s_ion = prop_tab['O+/H'].data
    O_d_ion = prop_tab['O++/H'].data
    log_O_s = prop_tab['log(O+/H)'].data
    log_O_d = prop_tab['log(O++/H)'].data

    # Error calculation

    # Initialize Dictionary for flux_gpdf
    flux_pdf_dict = dict()
    flux_peak = dict()
    flux_lowhigh = dict()

    # flux_data = [OII_flux, Hbeta_flux, Hdelta_flux, Hgamma_flux, OIII4363_flux, OIII4959_flux, OIII5007_flux]
    # RMS_data = [OII_RMS, Hbeta_RMS, Hdelta_RMS, Hgamma_RMS, OIII4363_RMS, OIII4959_RMS, OIII5007_RMS]

    for aa, flux, rms in zip(range(len(flux_cols)), flux_cols, flux_rms_cols):
        flux_gpdf = random_pdf(flux_tab[flux], flux_tab[rms], seed_i=aa,
                               n_iter=1000)
        err, xpeak = compute_onesig_pdf(flux_gpdf, flux_tab[flux], usepeak=True)

        # Fill In Dictionary
        flux_pdf_dict[line_name[aa]] = flux_gpdf
        flux_peak[line_name[aa] + '_xpeak'] = xpeak
        flux_lowhigh[line_name[aa] + '_lowhigh_error'] = err

        flux_tab0[line_name[aa] + '_Flux_Gaussian'][detection] = xpeak

    # Edit ASCII Table
    new_flux_file = join(path, filename_dict['bin_fit_rev'])

    if exists(new_flux_file):
        print("Overwriting: "+new_flux_file)
    else:
        print("Writing: "+new_flux_file)
    asc.write(flux_tab0, new_flux_file, overwrite=True, format='fixed_width_two_line')

    # Save npz files
    np.savez(path + 'flux_propdist.npz', **flux_pdf_dict)
    np.savez(path + 'flux_errors.npz', **flux_lowhigh)
    np.savez(path + 'flux_peak.npz', **flux_peak)
    # np.savez(path + 'Te_errors.npz', **Te_lowhigh)

    # Obtain distributions of line ratios: logR23, logO32, two_beta, three_beta, R
    flux_ratios_dict = error_prop_flux_ratios(flux_pdf_dict)

    Te_dict = temp_calculation(flux_ratios_dict['R'])
    com_O_log_pdf, metal_dict = \
        metallicity_calculation(Te_dict, flux_ratios_dict['two_beta'],
                                flux_ratios_dict['three_beta'])
