from . import line_name

# These are common/general column names

# Column names for bin information
bin_names0 = ['bin_ID', 'N_stack', 'Detection']

# Column names for individual galaxies/spectra
indv_names0 = ['ID', 'R23', 'O32']

# Column names for bin information in stellar mass and H-beta luminosity
bin_mzevolve_names0 = ['logM_min', 'logM_max', 'logM_avg', 'logLHb_min', 'logLHb_max', 'logLHb_avg']

# Column names for bin information in R23 and O32 line ratios
bin_zcalbase_names0 = ['logR23_min', 'logR23_max', 'logR23_avg', 'logO32_min', 'logO32_max', 'logO32_avg']

# Column names for Gaussian fitting
gauss_names0 = ['Flux_Gaussian', 'Flux_Observed', 'S/N', 'Sigma', 'Norm']  # This is just the suffix
gauss_lines_names0 = []
for line0 in line_name:
    gauss_lines_names0 += ['{}_{}'.format(line0, suffix) for suffix in gauss_names0]

# Temperature and metallicity properties
temp_metal_names0 = ['T_e', '12+log(O/H)', 'log(O+/H)', 'log(O++/H)', 'O+/H', 'O++/H']

# Dictionary containing filenames
filename_dict = {}
filename_dict['bin_info'] = 'bin_info.tbl'
filename_dict['bin_fit'] = 'bin_emission_line_fit.tbl'
filename_dict['bin_fit_revised'] = filename_dict['bin_fit'].replace('.tbl', '.revised.tbl')


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
