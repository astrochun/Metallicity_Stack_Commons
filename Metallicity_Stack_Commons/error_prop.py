from os.path import join, exists

from chun_codes import random_pdf, compute_onesig_pdf
from . import line_name
from astropy.io import ascii as asc
import numpy as np

from .column_names import filename_dict, temp_metal_names0, npz_filename_dict
from .ratios import error_prop_flux_ratios
from .temp_metallicity_calc import temp_calculation, metallicity_calculation


def write_npz(path, npz_files, dict_list):
    for file, dict_input in zip(npz_files, dict_list):
        npz_outfile = join(path, file)
        if exists(npz_outfile):
            print("Overwriting : "+npz_outfile)
        else:
            print("Writing : "+npz_outfile)
        np.savez(npz_outfile, **dict_input)


def fluxes_derived_prop(path, binned_data=True):
    """
    Purpose:
      Use measurements and their uncertainties to perform a randomization
      approach. The randomization is performed on individual emission lines.
      It carries that information to derived flux ratios and then
      determines electron temperature and metallicity

    :param path: str of full path
    :param binned_data: bool for whether to analysis binned data. Default: True

    :return:
    """

    # Define files to read in for binned data
    if binned_data:
        flux_file = join(path, filename_dict['bin_fit'])
        prop_file = join(path, filename_dict['bin_derived_prop'])
        verify_file = join(path, filename_dict['bin_valid'])

    flux_tab0 = asc.read(flux_file)
    prop_tab0 = asc.read(prop_file)

    if binned_data:
        verify_tab = asc.read(verify_file)
        detection = verify_tab['Detection'].data

        # For now we are only considering those with reliable detection and
        # excluding those with reliable non-detections (detect = 0.5)
        detect_idx = np.where((detection == 1))[0]

        ID = verify_tab['bin_ID'].data
        ID_detect = ID[detect_idx]
        print(ID_detect)

        flux_tab = flux_tab0[detect_idx]
        prop_tab = prop_tab0[detect_idx]

    flux_cols     = [str0+'_Flux_Gaussian' for str0 in line_name]
    flux_rms_cols = [str0+'_RMS' for str0 in line_name]

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
        err, peak = compute_onesig_pdf(flux_gpdf, flux_tab[flux], usepeak=True)

        # Fill In Dictionary
        flux_pdf_dict[line_name[aa]] = flux_gpdf
        flux_peak[line_name[aa] + '_peak'] = peak
        flux_lowhigh[line_name[aa] + '_lowhigh_error'] = err

        flux_tab0[line_name[aa] + '_Flux_Gaussian'][detect_idx] = peak

    # Edit ASCII Table
    new_flux_file = join(path, filename_dict['bin_fit_rev'])

    if exists(new_flux_file):
        print("Overwriting: "+new_flux_file)
    else:
        print("Writing: "+new_flux_file)
    asc.write(flux_tab0, new_flux_file, overwrite=True, format='fixed_width_two_line')

    # Save npz files
    npz_files = [npz_filename_dict['flux_pdf'],
                 npz_filename_dict['flux_errors'],
                 npz_filename_dict['flux_peak']]
    dict_list = [flux_pdf_dict, flux_lowhigh, flux_peak]
    write_npz(path, npz_files, dict_list)
    # np.savez(path + 'Te_errors.npz', **Te_lowhigh)

    # Obtain distributions of line ratios: logR23, logO32, two_beta, three_beta, R
    flux_ratios_dict = error_prop_flux_ratios(flux_pdf_dict)

    Te_dict = temp_calculation(flux_ratios_dict['R'])
    metal_dict = metallicity_calculation(Te_dict, flux_ratios_dict['two_beta'],
                                         flux_ratios_dict['three_beta'])

    # Loop for each derived properties (T_e, metallicity, etc.)
    metal_error = dict()
    metal_peak = dict()
    for names0 in temp_metal_names0:
        arr0 = prop_tab[names0].data

        pdf_arr = Te_dict if names0 == 'T_e' else metal_dict[names0]
        err_prop, peak_prop = compute_onesig_pdf(pdf_arr, arr0, usepeak=True)

        metal_error[names0+'_lowhigh_error'] = err_prop
        metal_peak[names0+'_peak'] = peak_prop

    npz_files = [npz_filename_dict['metal_errors'],
                 npz_filename_dict['metal_peak']]
    dict_list = [metal_error, metal_peak]
    write_npz(path, npz_files, dict_list)
