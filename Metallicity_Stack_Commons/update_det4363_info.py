from chun_codes import match_nosort
from astropy.io import fits


def get_index(det4363_table, input_table, column_name):
    """
    Uses either OBJNO or AP/SLIT info to get index for an existing table

    Inputs:
    det4363_table: astropy table containing DEEP2 [OIII]4363-detected sample

    input_table: astropy table containing the entire sample to be updated

    Return
    Numpy index arrays containing for det4363_table and input_table
    """

    if column_name != 'OBJNO' or column_name != 'SLIT' or column_name != 'AP':
        print("column_name not understood")
        print("Exiting!!!")
        return

    det4363_id = det4363_table[column_name]
    input_id   = input_table[column_name]

    det4363_idx, input_idx = match_nosort(det4363_id, input_id)

    return det4363_idx, input_idx
