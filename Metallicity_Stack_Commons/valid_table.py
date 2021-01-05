# from os.path import exists
import numpy as np
from os.path import join, exists

from astropy.io import ascii as asc
from astropy.table import Table, Column

from .column_names import filename_dict, valid_table_names0  # , bin_names0, remove_from_list


def make_validation_table(fitspath):
    """
    Purpose:
        This function creates a validation table for a given binning set. The validation table
        contains a OIII4363 detection column where 1.0 means detection, 0.5 means non-detection with
        reliable OIII5007, and 0.0 means unreliable non-detection.
        This function will be run every time the analysis is completed and will create a validation
        table for every analysis.

    Usage:
        valid_table.make_validation_table(fitspath, bin_type_str)

    Params:
        fitspath --> a string of the file path where the input file is and where the output file
            will be placed.
        not in use:
        bin_type_str --> a string describing the binning type. (e.g. 'massLHbetabin' or 'massbin')
                     --> This is dataset for Zcalbase_gal analysis

    Returns:
        None

    Outputs:
        fitspath + 'bin_validation.tbl' --> a validation table containing bin IDs;
                                        number of galaxies in each bin; and
                                        column indicating  OIII4363 detection/non-detection
                                        OIII4363_Flux_Observed
                                        OIII4363_S/N
    """

    bin_table = asc.read(join(fitspath, filename_dict['bin_info']))
    em_table = asc.read(join(fitspath, filename_dict['bin_fit']))

    bin_ID = em_table['bin_ID'].data
    raw_OIII4363 = em_table['OIII_4363_Flux_Observed'].data
    O_4363_SN = em_table['OIII_4363_S/N'].data
    O_4363_sigma = em_table['OIII_4363_Sigma'].data
    O_5007_SN = em_table['OIII_5007_S/N'].data
    
    N_stack = bin_table['N_stack'].data
    Hgamma_SN = em_table['HGAMMA_S/N'].data
    Hgamma = em_table['HGAMMA_Flux_Observed'].data

    detection  = np.zeros(len(bin_ID))
    OIII4363 = np.zeros(len(bin_ID))
    up_limit = (Hgamma/Hgamma_SN) * 3

    valid_stacks_idx = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma < 1.6))[0] 
    reliable_5007_stacks = np.where((O_4363_SN < 3) & (O_5007_SN > 100))[0]
    wide_lines_valid = np.where((O_4363_SN >= 3) & (O_5007_SN > 100) & (O_4363_sigma >= 1.6))[0]
    detection[valid_stacks_idx] = 1
    detection[reliable_5007_stacks] = 0.5
    detection[wide_lines_valid] = 0.5
    print(detection)

    for ii in range(len(OIII4363)):
        if detection[ii] == 1:
            OIII4363[ii] = raw_OIII4363[ii]
        if detection[ii] == 0.5:
            OIII4363[ii] = up_limit[ii]
        if detection[ii] == 0:
            OIII4363[ii] = up_limit[ii]

    ver_tab = fitspath + filename_dict['bin_valid']
    tab1 = Table([bin_ID, N_stack, detection, OIII4363, O_4363_SN], names=valid_table_names0)
    asc.write(tab1, ver_tab, format='fixed_width_two_line')
    '''
    # Write revised file for human editing
    ver_tab_revised = fitspath + filename_dict['bin_valid_rev']
    if not exists(ver_tab_revised):
        asc.write(tab1, ver_tab_revised, format='fixed_width_two_line')
        print("   ")
        print("URGENT!!! HUMAN EDITING OF FILE NEEDED : "+ver_tab_revised)
        print("   ")
    else:
        print("   ")
        print("ERROR!!! FILE EXISTS!!!  WILL NOT OVERWRITE !!!")
        print("ERROR!!! PLEASE RENAME/DELETE FILE TO REGENERATE !!!")
        print("   ")
    '''


def compare_to_by_eye(fitspath, dataset):
    """
    Purpose -> This function takes the automated validation table and checks it against
               inputted measurement that are determined by eye. These inputted measurements
               are in the np.where statements. It outputs a revised validation table based
               on the inputted measurements.

    Usage -> valid_table.make_validation_table(fitspath, dataset)

    Params:
        fitspath --> a string of the file path where the input file is and where the output file
            will be placed
        dataset  --> a string that is used to determine which eye measurements to use

    Returns:
        None

    Outputs:
        fitspath + 'bin_validation_revised.tbl' and '.csv' -->
                                        a validation table containing bin IDs;
                                        number of galaxies in each bin; and
                                        column indicating  OIII4363 detection/non-detection
                                        OIII4363_Flux_Observed
                                        OIII4363_S/N
                                        Notes

    """

    valid_rev_file = join(fitspath, filename_dict['bin_valid_rev'])
    if exists(valid_rev_file):
        print("!!! Revised validation table exists. Not overwriting! : ", valid_rev_file)
        raise FileExistsError
    else:
        valid_file = join(fitspath, filename_dict['bin_valid'])
        valid_tab = asc.read(valid_file)
        indicate = valid_tab['Detection']
        ID = valid_tab['bin_ID']

        # Will need to provide a config file that contains the manual inspections

        # Will need to modify table Detection Column

        asc.write(valid_tab, valid_rev_file, format='fixed_width_two_line')
