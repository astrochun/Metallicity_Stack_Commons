from os.path import join

from chun_codes import random_pdf, compute_onesig_pdf
from . import line_name
from astropy.io import ascii as asc
import numpy as np

from .column_names import filename_dict


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


def error_prop_chuncodes(path):
    flux_file = join(path, filename_dict['bin_fit'])
    flux_tab0  = asc.read(flux_file)

    prop_file = join(path, filename_dict['bin_derived_prop'])
    prop_tab0  = asc.read(prop_file)

    verify_file = join(path, filename_dict['bin_valid'])
    verify_tab = asc.read(verify_file)

    detect    = verify_tab['Detection']
    detection = np.where((detect == 1))[0]

    ID = verify_tab['ID']
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
    flux_propdist_dict = {}
    flux_xpeak = {}
    flux_lowhigh = {}

    # flux_data = [OII_flux, Hbeta_flux, Hdelta_flux, Hgamma_flux, OIII4363_flux, OIII4959_flux, OIII5007_flux]
    # RMS_data = [OII_RMS, Hbeta_RMS, Hdelta_RMS, Hgamma_RMS, OIII4363_RMS, OIII4959_RMS, OIII5007_RMS]

    for aa, flux, rms in zip(range(len(flux_cols)), flux_cols, flux_rms_cols):
        flux_gpdf = random_pdf(flux_tab[flux], flux_tab[rms], seed_i=aa,
                               n_iter=1000, silent=False)
        err, xpeak = compute_onesig_pdf(flux_gpdf, flux_tab[flux], usepeak=True,
                                        silent=True, verbose=True)

        # Fill In Dictionary
        flux_propdist_dict[line_name[aa]] = flux_gpdf
        flux_xpeak[line_name[aa] + '_xpeak'] = xpeak
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
    np.savez(path + 'flux_propdist.npz', **flux_propdist_dict)
    np.savez(path + 'flux_errors.npz', **flux_lowhigh)
    np.savez(path + 'flux_xpeak.npz', **flux_xpeak)
    # np.savez(path + 'Te_errors.npz', **Te_lowhigh)

