from . import line_name, line_type
import sys

# Get python version
py_vers = sys.version_info.major


# Need to define here
def line_fit_suffix_add(line_name0, line_type0):
    """
    Purpose:
      Simple list comprehension combining emission line fit suffixes with
      the emission line.  This works for individual lines

    :param line_name0: str containing the line name
    :param line_type0: str containing the emisison-line type (e.g., 'Balmer')
    :return: gauss_lines_names: list with str formatted as [LINE]_[SUFFIX]
    """

    gauss_lines_names = ['{}_{}'.format(line_name0, suffix) for suffix in gauss_names0]
    if line_type0 == 'Balmer':
        gauss_lines_names += ['{}_{}'.format(line_name0, suffix) for suffix in balmer_names0]

    return gauss_lines_names


# These are common/general column names

# Column names for bin information
bin_names0 = ['bin_ID', 'N_stack', 'Detection']

# Column names for individual galaxies/spectra
indv_names0 = ['ID', 'logR23', 'logO32', 'logM', 'logLHb', 'two_beta', 'three_beta', 'R']

# Dust attenuation
dust0 = ['HgHb', 'HdHb', 'EBV_HgHb', 'EBV_HdHb']

# Column names for bin information in stellar mass and H-beta luminosity
bin_mzevolve_names0 = ['logM_min', 'logM_max', 'logM_avg', 'logM_median',
                       'logLHb_min', 'logLHb_max', 'logLHb_avg', 'logLHb_median']

# Column names for bin information in R23 and O32 line ratios
bin_zcalbase_names0 = ['logR23_min', 'logR23_max', 'logR23_avg', 'logR23_median',
                       'logO32_min', 'logO32_max', 'logO32_avg', 'logO32_median']

# Column names for composite line ratios
bin_ratios0 = ['logR23_composite', 'logO32_composite',
               'two_beta_composite', 'three_beta_composite', 'R_composite']

# Column names for Gaussian fitting
# This is just the suffix
gauss_names0 = ['Flux_Gaussian', 'Flux_Observed', 'S/N', 'Center', 'Norm',
                'Median', 'Sigma', 'RMS']
balmer_names0 = ['Abs_Norm', 'Abs_Sigma']

# Emission-line fit column names with [LINE] prefix
gauss_lines_names0 = []
for line0, type0 in zip(line_name, line_type):
    gauss_lines_names0 += line_fit_suffix_add(line0, type0)


# Temperature and metallicity properties
temp_metal_names0 = ['T_e', '12+log(O/H)', 'log(O+/H)', 'log(O++/H)', 'O+/H', 'O++/H']

# Validation Table
valid_table_names0 = ['bin_ID','N_stack','Detection', 'OIII_4363_Flux_Observed', 'OIII_4363_S/N']

# Dictionary containing filenames
filename_dict = dict()

# Bin-related files
filename_dict['comp_spec'] = 'composite_spectra.fits'
filename_dict['bin_info'] = 'bin_info.tbl'
filename_dict['bin_valid'] = 'bin_validation.tbl'
filename_dict['bin_valid_rev'] = filename_dict['bin_valid'].replace('.tbl', '.revised.tbl')
filename_dict['bin_fit'] = 'bin_emission_line_fit.tbl'
filename_dict['bin_fit_rev'] = filename_dict['bin_fit'].replace('.tbl', '.MC.tbl')
filename_dict['bin_derived_prop'] = 'bin_derived_properties.tbl'
filename_dict['bin_derived_prop_rev'] = filename_dict['bin_derived_prop'].replace('.tbl', '.MC.tbl')
filename_dict['bin_derived_prop_dust'] = filename_dict['bin_derived_prop'].replace('.tbl', '.dustcorr.tbl')
filename_dict['bin_derived_prop_rev_dust'] = filename_dict['bin_derived_prop_rev'].replace('.tbl', '.dustcorr.tbl')

# Individual galaxy/spectra-related files
filename_dict['indv_prop'] = 'individual_properties.tbl'
filename_dict['indv_bin_info'] = 'individual_bin_info.tbl'
filename_dict['indv_derived_prop'] = 'individual_derived_properties.tbl'

# numpy files
npz_filename_dict = dict()
npz_filename_dict['flux_pdf'] = 'flux_pdf.npz'
npz_filename_dict['flux_errors'] = 'flux_errors.npz'
npz_filename_dict['flux_peak'] = 'flux_peak.npz'
npz_filename_dict['der_prop_pdf'] = 'derived_properties_pdf.npz'
npz_filename_dict['der_prop_errors'] = 'derived_properties_errors.npz'
npz_filename_dict['der_prop_peak'] = 'derived_properties_peak.npz'


def merge_column_names(*args):
    """
    Purpose:
      Merges multiple lists containing column names.

    Usage:
      column_names = merge_column_names(bin_names0, indv_names0)

    :param args: An undefined number of lists
    :return merge_list:
    """

    merge_list = list()

    arg_count = len(args)
    if arg_count > 0:
        for elem in args:
            merge_list += elem

    return merge_list


def remove_from_list(list0, remove_entries):
    """
    Purpose:
      Remove entries from list

    :param list0: list of column names
    :param remove_entries: list of column names to remove
    """

    if py_vers == 3:
        dup_list0 = list0.copy()
    if py_vers == 2:
        dup_list0 = list(list0)

    for entry in remove_entries:
        dup_list0.remove(entry)

    return dup_list0


def indv_R23_O32():
    """
    Purpose:
      Use remove_from_list() to provide simplified list that contains ID, logR23 and logO32

    :return: list containing just ID, logR23, logO32
    """

    return remove_from_list(indv_names0, ['logM', 'logLHb'])


def indv_M_LHb():
    """
    Purpose:
      Use remove_from_list() to provide simplified list that contains ID, logM and logLHb

    :return: list containing just ID, logM, logLHb
    """

    return remove_from_list(indv_names0, ['logR23', 'logO32'])
