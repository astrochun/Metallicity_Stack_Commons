from os.path import join
from os.path import exists

import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table

from .temp_metallicity_calc import metallicity_calculation
from ..column_names import bin_names0, indv_names0, temp_metal_names0
from ..column_names import filename_dict
from .ratios import flux_ratios
from .. import line_name
from ..logging import log_stdout, log_verbose

ID_name = indv_names0[0]
bin_ID_name = bin_names0[0]

logR23_name = indv_names0[1]
logO32_name = indv_names0[2]
two_beta_name = indv_names0[5]
three_beta_name = indv_names0[6]


def main(fitspath, dataset, revised=False, det3=True, verbose=False, log=None):
    """
    Purpose:
      Reads in composite table(s) containing bin information to
      determine temperature-based metallicity from composite average
      T_e and individual line ratios ([OII]/H-beta, [OIII]/H-beta)

    :param fitspath: str containing folder path
    :param dataset: str containing sub-folder (specific to stacking approach)
    :param revised: Bool indicates whether to use revised bin properties
                    (e.g., *.revised.tbl files)
    :param det3: Bool indicates whether individual galaxy files is limited to
                 those satisfying emission-line det3 requirement
                 Default: True
    :param verbose: bool to write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    Files identified by default
    composite_file: str containing filename of composite data
      e.g., '[dataset]/bin_derived_properties.tbl' or
            '[dataset]/bin_derived_properties.revised.tbl'
    indv_em_line_file: str containing filename that contains
                     emission-line information for each galaxy
      e.g., 'individual_properties.tbl'
    indv_bin_file: str containing filename tha contains bin information
                   for each galaxy
      e.g., '[dataset]/individual_bin_info.tbl'
    outfile: str containing filename of output file
      e.g., '[dataset]/individual_derived_properties.tbl'
    """

    if log is None:
        log = log_stdout()

    log_verbose(log, "starting ...", verbose=verbose)

    # Define [composite_file]
    t_comp = filename_dict['bin_derived_prop'] if not revised else \
        filename_dict['bin_derived_prop']
    composite_file = join(fitspath, dataset, t_comp)
    if not exists(composite_file):
        log.warning(f"ERROR: File not found! {composite_file}")
        return

    # Read in composite table
    log.info(f"Reading: {composite_file}")
    composite_table = asc.read(composite_file)

    bin_id = composite_table['bin_ID'].data
    bin_temp = composite_table['T_e'].data

    # Read in validation table
    valid_file = join(fitspath, dataset, filename_dict['bin_valid'])
    if not exists(valid_file):
        log.warning(f"ERROR: File not found! {valid_file}")
        return
    log.info(f"Reading: {valid_file}")
    valid_table = asc.read(valid_file)

    # Define [indv_em_line_file]
    indv_em_line_file = join(fitspath, dataset, filename_dict['indv_prop'])
    if not exists(indv_em_line_file):
        log.warning(f"ERROR: File not found! {indv_em_line_file}")
        return

    # Read in tables containing line ratios, etc.
    log.info(f"Reading: {indv_em_line_file}")
    indv_em_line_table = asc.read(indv_em_line_file)

    # Define [indv_bin_file]
    indv_bin_file = join(fitspath, dataset, filename_dict['indv_bin_info'])
    if not exists(indv_bin_file):
        log.warning(f"ERROR: File not found! {indv_bin_file}")
        return

    # Read in tables containing bin info for individual
    log.info(f"Reading: {indv_bin_file}")
    indv_bin_info_table = asc.read(indv_bin_file)

    # Populate composite temperature for individual galaxies
    adopted_temp = np.zeros(len(indv_em_line_table))
    bin_id_indv = np.zeros(len(indv_em_line_table))
    detect_indv = np.zeros(len(indv_em_line_table))
    for comp_bin, comp_temp, detect in zip(bin_id, bin_temp,
                                           valid_table['Detection']):
        bin_idx = np.where(indv_bin_info_table['bin_ID'].data == comp_bin)[0]
        adopted_temp[bin_idx] = comp_temp
        bin_id_indv[bin_idx]  = comp_bin
        detect_indv[bin_idx]  = detect

    flux_dict = dict()
    for line in line_name:
        flux_dict[line] = indv_em_line_table[line+'_Flux_Gaussian'].data

    flux_ratios_dict = flux_ratios(flux_dict, get_R=False)

    two_beta = flux_ratios_dict[two_beta_name]
    three_beta = flux_ratios_dict[three_beta_name]
    logR23 = flux_ratios_dict[logR23_name]
    logO32 = flux_ratios_dict[logO32_name]

    if not det3:
        metal_dict = \
            metallicity_calculation(adopted_temp, two_beta, three_beta)
    else:
        det3_idx = np.where((detect_indv == 1.0) | (detect_indv == 0.5))[0]
        metal_dict = \
            metallicity_calculation(adopted_temp, two_beta, three_beta,
                                    det3=det3_idx)

    # Define [indv_derived_prop_table] to include ID, bin_ID, composite T_e,
    # and 12+log(O/H)
    arr0 = [indv_em_line_table[ID_name], bin_id_indv, logR23, logO32,
            two_beta, three_beta, adopted_temp]
    names0 = [ID_name, bin_ID_name, logR23_name, logO32_name, two_beta_name,
              three_beta_name, temp_metal_names0[0]]

    # Include other metallicities
    arr0 += list(metal_dict.values())
    names0 += metal_dict.keys()
    indv_derived_prop_table = Table(arr0, names=names0)

    outfile = join(fitspath, dataset, filename_dict['indv_derived_prop'])

    # Write Astropy ASCII table containing composite T_e and derived metallicity
    if exists(outfile):
        log.info(f"File exists! Overwriting : {outfile}")
    else:
        log.info(f"Writing : {outfile}")
    indv_derived_prop_table.write(outfile, overwrite=True,
                                  format='ascii.fixed_width_two_line')

    log_verbose(log, "finished.", verbose=verbose)
