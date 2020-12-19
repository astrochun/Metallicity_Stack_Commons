from os.path import join, exists

from chun_codes import random_pdf, compute_onesig_pdf
from .. import line_name
from ..logging import log_stdout

from astropy.io import ascii as asc
from astropy.table import Table, hstack
import numpy as np

from ..column_names import filename_dict, temp_metal_names0, npz_filename_dict, \
    bin_names0, bin_ratios0, dust0
from .ratios import flux_ratios
from .temp_metallicity_calc import temp_calculation, metallicity_calculation
from .attenuation import compute_EBV


def write_npz(path, npz_files, dict_list, log=None):
    """
    Purpose:
      Write numpy files with provided dictionaries

    :param path: str - prefix for filename output
    :param npz_files: list - contains npz file names
    :param dict_list: list - contains dictionaries for each corresponding npz file
    :param log: LogClass or logging object

    :return: Write npz files
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    for file, dict_input in zip(npz_files, dict_list):
        npz_outfile = join(path, file)
        if exists(npz_outfile):
            log.info(f"Overwriting : {npz_outfile}")
        else:
            log.info(f"Writing : {npz_outfile}")
        np.savez(npz_outfile, **dict_input)

    log.debug("finished ...")


def fluxes_derived_prop(path, raw=False, binned_data=True, apply_dust=False,
                        revised=True, log=None):
    """
    Purpose:
      Use measurements and their uncertainties to perform a randomization
      approach. The randomization is performed on individual emission lines.
      It carries that information to derived flux ratios and then
      determines electron temperature and metallicity

    :param path: str of full path
    :param raw: bool to do a simple calculation, no randomization. Default: False
    :param binned_data: bool for whether to analysis binned data. Default: True
    :param apply_dust: bool for whether to apply dust attenuation. Default: False
    :param revised: bool to indicate if revised validation table is used. Default: True
    :param log: LogClass or logging object
    """

    if log is None:
        log = log_stdout()

    log.debug("starting ...")

    # Define files to read in for binned data
    if binned_data:
        flux_file = join(path, filename_dict['bin_fit'])
        bin_ratios = bin_ratios0
        dust = dust0
    else:
        bin_ratios = [ratios0.replace('_composite', '')
                      for ratios0 in bin_ratios0]
        dust = [t_dust.replace('_composite', '') for t_dust in dust0]

    # Define verification table
    if revised:
        rev_s = ''
        log.info("Using REVISED validation table")
        verify_file = join(path, filename_dict['bin_valid_rev'])
    else:
        rev_s = '_v1'
        log.info("Using validation table")
        verify_file = join(path, filename_dict['bin_valid'])

    if not raw:
        # Set rev_stuff based on revised=True/False
        prop_file = join(path, filename_dict['bin_derived_prop' + rev_s])

        log.info(f"Reading : {prop_file}")
        prop_tab0 = asc.read(prop_file)

    two_beta_name = bin_ratios[2]
    three_beta_name = bin_ratios[3]
    R_name = bin_ratios[4]

    log.info(f"Reading : {flux_file}")
    flux_tab0 = asc.read(flux_file)

    if binned_data:
        if not raw:
            log.info(f"Reading : {verify_file}")
            verify_tab = asc.read(verify_file)
            detection = verify_tab['Detection'].data

            # For now we are only considering those with reliable detection and
            # excluding those with reliable non-detections (detect = 0.5)
            detect_idx = np.where((detection == 1))[0]

            ID = verify_tab['bin_ID'].data
            ID_detect = ID[detect_idx]
            log.info(ID_detect)

            flux_tab = flux_tab0[detect_idx]
            prop_tab = prop_tab0[detect_idx]
    else:
        log.warning("Individual sources not supported yet")
        return

    flux_cols     = [str0+'_Flux_Gaussian' for str0 in line_name]
    flux_rms_cols = [str0+'_RMS' for str0 in line_name]

    #
    # EMISSION-LINE SECTION
    #

    if raw:
        flux_dict = dict()
        for aa, flux in zip(range(len(flux_cols)), flux_cols):

            # Fill in dictionary
            flux_dict[line_name[aa]] = flux_tab0[flux].data

        flux_ratios_dict = flux_ratios(flux_dict, binned_data=binned_data,
                                       log=log)

        #
        # Get EBV from Balmer decrement if apply_dust set
        #

        if apply_dust:
            EBV = compute_EBV(flux_ratios_dict[dust[0]], source=dust[0],
                              log=log)
            EBV_HdHb = compute_EBV(flux_ratios_dict[dust[1]], source=dust[1],
                                   log=log)
        else:
            EBV = None
            EBV_HdHb = None

        derived_prop_dict = dict()

        # Using Hg/Hb Balmer decrement for dust corrections
        Te = temp_calculation(flux_ratios_dict[R_name], EBV=EBV, log=log)
        derived_prop_dict[temp_metal_names0[0]] = Te
        metal_dict = metallicity_calculation(Te, flux_ratios_dict[two_beta_name],
                                             flux_ratios_dict[three_beta_name],
                                             EBV=EBV, log=log)
        derived_prop_dict.update(metal_dict)

        if apply_dust:
            derived_prop_dict[dust[2]] = EBV
            derived_prop_dict[dust[3]] = EBV_HdHb

        # Write Astropy table
        # Get bin IDs and composite measurements in one dictionary and write to table
        tbl_dict = {bin_names0[0]: np.int_(flux_tab0[bin_names0[0]].data)}
        tbl_dict.update(flux_ratios_dict)
        tbl_dict.update(derived_prop_dict)

        if not apply_dust:
            outfile = join(path, filename_dict['bin_derived_prop' + rev_s])
        else:
            outfile = join(path, filename_dict['bin_derived_prop_dust' + rev_s])
        if not exists(outfile):
            log.info(f"Writing : {outfile}")
        else:
            log.info(f"Overwriting : {outfile}")
        asc.write(tbl_dict, output=outfile, overwrite=True,
                  format='fixed_width_two_line')
    else:
        # Initialize dictionaries

        flux_pdf_dict = dict()
        flux_peak_dict = dict()
        flux_error_dict = dict()

        # Randomization for emission-line fluxes
        for aa, flux, rms in zip(range(len(flux_cols)), flux_cols,
                                 flux_rms_cols):
            flux_pdf = random_pdf(flux_tab[flux], flux_tab[rms],
                                  seed_i=aa+1, n_iter=1000)
            err, peak = compute_onesig_pdf(flux_pdf, flux_tab[flux],
                                           usepeak=True)

            # Fill in dictionary
            flux_pdf_dict[line_name[aa]] = flux_pdf
            flux_peak_dict[line_name[aa] + '_peak'] = peak
            flux_error_dict[line_name[aa] + '_error'] = err

            # Update values
            flux_tab0[line_name[aa] + '_Flux_Gaussian'][detect_idx] = peak

        # Write revised emission-line fit ASCII table
        new_flux_file = join(path, filename_dict['bin_fit_MC'])
        if exists(new_flux_file):
            log.info(f"Overwriting: {new_flux_file}")
        else:
            log.info(f"Writing: {new_flux_file}")
        asc.write(flux_tab0, new_flux_file, overwrite=True,
                  format='fixed_width_two_line')

        # Save flux npz files
        npz_files = [npz_filename_dict['flux_pdf' + rev_s],
                     npz_filename_dict['flux_errors' + rev_s],
                     npz_filename_dict['flux_peak' + rev_s]]
        dict_list = [flux_pdf_dict, flux_error_dict, flux_peak_dict]
        write_npz(path, npz_files, dict_list)

        # Obtain line ratio distributions: logR23, logO32, two_beta, three_beta, R
        # Also include Balmer decrements
        flux_ratios_dict = flux_ratios(flux_pdf_dict, binned_data=binned_data,
                                       log=log)

        flux_ratios_err_dict = dict()

        for name0 in flux_ratios_dict.keys():
            err_prop, peak_prop = compute_onesig_pdf(flux_ratios_dict[name0],
                                                     prop_tab0[name0],
                                                     usepeak=True)

            # Update values
            prop_tab0[name0][detect_idx] = peak_prop

            flux_ratios_err_dict[name0 + '_low_err'] = np.repeat(np.nan, len(prop_tab0))
            flux_ratios_err_dict[name0 + '_high_err'] = np.repeat(np.nan, len(prop_tab0))
            flux_ratios_err_dict[name0 + '_low_err'][detect_idx] = err_prop[:, 0]
            flux_ratios_err_dict[name0 + '_high_err'][detect_idx] = err_prop[:, 1]

        # Add flux ratio errors to table
        flux_ratios_err_table = Table(flux_ratios_err_dict)
        prop_tab0 = hstack([prop_tab0, flux_ratios_err_table])

        #
        # Get EBV from Balmer decrement if apply_dust set
        #

        if apply_dust:
            # Hg/Hb case
            EBV_dict = {dust[2]: np.repeat(np.nan, len(prop_tab0)),
                        dust[2] + '_low_err': np.repeat(np.nan, len(prop_tab0)),
                        dust[2] + '_high_err': np.repeat(np.nan, len(prop_tab0))}

            EBV, EBV_peak = compute_EBV(flux_ratios_dict[dust[0]],
                                        source=dust[0], log=log)

            err_prop, peak_prop = compute_onesig_pdf(EBV, EBV_peak,
                                                     usepeak=True)

            EBV_dict[dust[2]][detect_idx] = peak_prop
            EBV_dict[dust[2] + '_low_err'][detect_idx] = err_prop[:, 0]
            EBV_dict[dust[2] + '_high_err'][detect_idx] = err_prop[:, 1]
            EBV_table = Table(EBV_dict)  # Add in at the end

            # Hd/Hb case
            EBV_HdHb_dict = {dust[3]: np.repeat(np.nan, len(prop_tab0)),
                             dust[3] + '_low_err': np.repeat(np.nan, len(prop_tab0)),
                             dust[3] + '_high_err': np.repeat(np.nan, len(prop_tab0))}

            EBV_HdHb, EBV_HdHb_peak = compute_EBV(flux_ratios_dict[dust[1]],
                                                  source=dust[1], log=log)
            log.info(EBV_HdHb_peak)

            err_prop, peak_prop = compute_onesig_pdf(EBV_HdHb, EBV_HdHb_peak,
                                                     usepeak=True)

            EBV_HdHb_dict[dust[3]][detect_idx] = peak_prop
            EBV_HdHb_dict[dust[3] + '_low_err'][detect_idx] = err_prop[:, 0]
            EBV_HdHb_dict[dust[3] + '_high_err'][detect_idx] = err_prop[:, 1]

            # Add EBV and errors to table
            EBV_HdHb_table = Table(EBV_HdHb_dict)  # Add in at the end
        else:
            EBV = None

        #
        # DERIVED PROPERTIES SECTION
        #

        # Initialize dictionaries
        derived_prop_pdf_dict = dict()
        derived_prop_error_dict = dict()
        derived_prop_peak_dict = dict()

        # Calculate temperature distribution
        Te_pdf = temp_calculation(flux_ratios_dict[R_name], EBV=EBV, log=log)
        derived_prop_pdf_dict[temp_metal_names0[0]] = Te_pdf

        # Calculate metallicity distribution
        metal_dict = metallicity_calculation(Te_pdf,
                                             flux_ratios_dict[two_beta_name],
                                             flux_ratios_dict[three_beta_name],
                                             EBV=EBV, log=log)
        derived_prop_pdf_dict.update(metal_dict)

        prop_err_dict = dict()  # Initialize

        # Loop for each derived properties (T_e, metallicity, etc.)
        for names0 in temp_metal_names0:
            prop_err_dict[names0 + '_low_err'] = np.repeat(np.nan, len(prop_tab0))
            prop_err_dict[names0 + '_high_err'] = np.repeat(np.nan, len(prop_tab0))

            arr0 = prop_tab[names0].data

            err_prop, peak_prop = \
                compute_onesig_pdf(derived_prop_pdf_dict[names0], arr0,
                                   usepeak=True)

            # Fill in dictionary
            derived_prop_error_dict[names0+'_error'] = err_prop
            derived_prop_peak_dict[names0+'_peak'] = peak_prop

            prop_err_dict[names0 + '_low_err'][detect_idx] = err_prop[:, 0]
            prop_err_dict[names0 + '_high_err'][detect_idx] = err_prop[:, 1]

            # Update values
            prop_tab0[names0][detect_idx] = peak_prop

        # Add errors to table
        prop_err_table = Table(prop_err_dict)
        prop_tab0 = hstack([prop_tab0, prop_err_table])

        if apply_dust:
            prop_tab0 = hstack([prop_tab0, EBV_table])
            prop_tab0 = hstack([prop_tab0, EBV_HdHb_table])

        # Write revised properties ASCII table
        if not apply_dust:
            new_prop_file = join(path,
                                 filename_dict['bin_derived_prop_MC' + rev_s])
        else:
            new_prop_file = join(path,
                                 filename_dict['bin_derived_prop_MC_dust' + rev_s])
        if exists(new_prop_file):
            log.info(f"Overwriting: {new_prop_file}")
        else:
            log.info(f"Writing: {new_prop_file}")

        '''
        # Order of table:
        tab_order = [prop_tab0.colnames[0]]
        tab_order += flux_ratios_dict.keys()
        tab_order += temp_metal_names0
        if apply_dust:
            tab_order += EBV_dict.keys()
            tab_order += EBV_HdHb_dict.keys()

        if not raw:
            tab_order += prop_err_dict.keys()

        asc.write(prop_tab0[tab_order], new_prop_file, overwrite=True,
        '''
        asc.write(prop_tab0, new_prop_file, overwrite=True,
                  format='fixed_width_two_line')

        # Save derived properties npz files
        if not apply_dust:
            npz_files = [npz_filename_dict['der_prop_pdf' + rev_s],
                         npz_filename_dict['der_prop_errors' + rev_s],
                         npz_filename_dict['der_prop_peak' + rev_s]]
        else:
            npz_files = [npz_filename_dict['der_prop_dust_pdf' + rev_s],
                         npz_filename_dict['der_prop_dust_errors' + rev_s],
                         npz_filename_dict['der_prop_dust_peak' + rev_s]]

        dict_list = [derived_prop_pdf_dict, derived_prop_error_dict,
                     derived_prop_peak_dict]
        write_npz(path, npz_files, dict_list)

    log.debug("finished ...")
