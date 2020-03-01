from os.path import join
from os.path import exists

import numpy as np
from astropy.io import ascii as asc
from astropy.table import Column

from ..temp_metallicity_calc import metallicity_calculation
from .. import OIII_r


def main(fitspath, dataset, composite_file, indv_em_line_file, indv_bin_file, outfile, det3=True):
    """
    Purpose:
      Reads in composite table(s) containing bin information to
      determine temperature-based metallicity from composite average
      T_e and individual line ratios ([OII]/H-beta, [OIII]/H-beta)

    :param fitspath: str containing folder path
    :param dataset: str containing sub-folder (specific to stacking approach)
    :param composite_file: str containing filename of composite data
             e.g., '[dataset]/bin_derived_properties.tbl' or
                   '[dataset]/bin_derived_properties.revised.tbl'
    :param indv_em_line_file: str containing filename that contains
                              emission-line information for each galaxy
             e.g., 'individual_properties.tbl'
    :param indv_bin_file: str containing filename tha contains bin information
                          for each galaxy
             e.g., '[dataset]/individual_bin_info.tbl'
    :param outfile: str containing filename of output file
             e.g., '[dataset]/individual_derived_properties.tbl'
    :param det3: Bool indicates whether individual galaxy files is limited to
                 those satisfying emission-line det3 requirement
                 Default: True
    """

    # Read in composite table
    composite_table = asc.read(composite_file)

    bin_id = composite_table['bin_ID'].data
    bin_temp = composite_table['T_e'].data

    # Read in tables containing line ratios, bins, etc.
    indv_em_line_table = asc.read(join(fitspath, indv_em_line_file))
    indv_bin_info_table = asc.read(join(fitspath, dataset+indv_bin_file))
    # Not used for now
    # average_table = asc.read(join(fitspath, dataset+'_Average_R23_O32_Values.tbl'))

    # Populate composite temperature for individual galaxies
    adopted_temp = np.zeros(len(indv_em_line_table))
    for comp_bin, comp_temp in zip(bin_id, bin_temp):
        bin_idx = np.where(indv_bin_info_table['bin_ID'].data == comp_bin)[0]
        adopted_temp[bin_idx] = comp_temp

    O2 = indv_em_line_table['OII_3727_Flux_Observed'].data            # [OII]3726,3728 fluxes
    O3 = indv_em_line_table['OIII_5007_Flux_Observed'].data * OIII_r  # [OIII]4959,5007 fluxes (Assume 3.1:1 ratio)
    Hb = indv_em_line_table['HBETA_Flux_Observed'].data               # H-beta fluxes

    if det3:
        com_O_log, metal_dict = metallicity_calculation(adopted_temp, O2/Hb, O3/Hb)
    else:
        det3 = np.where(indv_bin_info_table['Individual_Detections'])[0]
        temp_com_O_log, temp_metal_dict = \
            metallicity_calculation(adopted_temp[det3], O2[det3]/Hb[det3],
                                    O3[det3]/Hb[det3])
        com_O_log = np.zeros(len(indv_em_line_table))
        com_O_log[det3] = temp_com_O_log

    # Update [indv_em_line_table] to include two new columns
    col_temp = Column(adopted_temp, name='T_e')
    col_metal = Column(com_O_log, name='12+log(O/H)')
    indv_em_line_table.add_columns([col_temp, col_metal])  # Add at the end (default)

    # Write Astropy ASCII table containing composite T_e and derived metallicity

    if exists(outfile):
        print("File exists! Overwriting : ", outfile)
    else:
        print("Writing : ", outfile)
    indv_em_line_table.write(outfile, overwrite=True, format='ascii.fixed_width_two_line')
