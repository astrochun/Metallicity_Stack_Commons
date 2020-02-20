# These are common/general column names

bin_names0 = ['bin_ID', 'N_stack', 'Detection']

indv_names0 = ['ID', 'R23', 'O32']

bin_mzevolve_names0 = ['logM_min', 'logM_max', 'logM_avg', 'logLHb_min', 'logLHb_max', 'logLHb_avg']

bin_zcalbase_names0 = ['logR23_min', 'logR23_max', 'logR23_avg', 'logO32_min', 'logO32_max', 'logO32_avg']


def merge_column_names(*args):
    """
    Purpose:
      Merges multiple lists containing column names.

    Usage:
      column_names = merge_column_names(bin_names0, indv_names0)

    :param args: An undefined number of lists
    :return merge_list:
    """

    merge_list = []

    arg_count = len(args)
    if arg_count > 0:
        for elem in args:
            merge_list += elem

    return merge_list
