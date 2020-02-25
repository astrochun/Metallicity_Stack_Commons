from . import line_name, line_type

# These are common/general column names

# Column names for bin information
bin_names0 = ['bin_ID', 'N_stack', 'Detection']

# Column names for individual galaxies/spectra
indv_names0 = ['ID', 'logR23', 'logO32', 'logM', 'logLHb']

# Column names for bin information in stellar mass and H-beta luminosity
bin_mzevolve_names0 = ['logM_min', 'logM_max', 'logM_avg', 'logM_median',
                       'logLHb_min', 'logLHb_max', 'logLHb_avg', 'logLHb_median']

# Column names for bin information in R23 and O32 line ratios
bin_zcalbase_names0 = ['logR23_min', 'logR23_max', 'logR23_avg', 'logR23_median',
                       'logO32_min', 'logO32_max', 'logO32_avg', 'logO32_median']

# Column names for Gaussian fitting
# This is just the suffix
gauss_names0 = ['Flux_Gaussian', 'Flux_Observed', 'S/N', 'Center', 'Norm',
                'Median', 'Sigma']
balmer_names0 = ['Abs_Norm', 'Abs_Sigma']
gauss_lines_names0 = []
for line0, type0 in zip(line_name, line_type):
    gauss_lines_names0 += ['{}_{}'.format(line0, suffix) for suffix in gauss_names0]
    if type0 == 'Balmer':
        gauss_lines_names0 += ['{}_{}'.format(line0, suffix) for suffix in balmer_names0]

# Temperature and metallicity properties
temp_metal_names0 = ['T_e', '12+log(O/H)', 'log(O+/H)', 'log(O++/H)', 'O+/H', 'O++/H']

# Dictionary containing filenames
filename_dict = dict()

# Bin-related files
filename_dict['bin_info'] = 'bin_info.tbl'
filename_dict['bin_valid'] = 'bin_validation.tbl'
filename_dict['bin_fit'] = 'bin_emission_line_fit.tbl'
filename_dict['bin_fit_rev'] = filename_dict['bin_fit'].replace('.tbl', '.revised.tbl')
filename_dict['bin_derived_prop'] = 'bin_derived_properties.tbl'
filename_dict['bin_derived_prop_rev'] = filename_dict['bin_derived_prop'].replace('.tbl', '.revised.tbl')

# Individual galaxy/spectra-related files
filename_dict['indv_prop'] = 'individual_properties.tbl'
filename_dict['indv_bin_info'] = 'individual_bin_info.tbl'
filename_dict['indv_derived_prop'] = 'individual_derived_properties.tbl'


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

    dup_list0 = list0.copy()

    for entry in remove_entries:
        dup_list0.remove(entry)

    return dup_list0
