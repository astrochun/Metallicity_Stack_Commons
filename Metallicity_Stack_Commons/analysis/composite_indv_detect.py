from os.path import join
from os.path import exists

import numpy as np
from astropy.io import ascii as asc
from astropy.table import Table

from ..temp_metallicity_calc import metallicity_calculation
from .. import OIII_r
from ..column_names import bin_names0, indv_names0, temp_metal_names0
from ..column_names import filename_dict

ID_name = indv_names0[0]
bin_ID_name = bin_names0[0]


def main(fitspath, dataset, outfile,
         revised=False, det3=True):
    """
    Purpose:
      Reads in composite table(s) containing bin information to
      determine temperature-based metallicity from composite average
      T_e and individual line ratios ([OII]/H-beta, [OIII]/H-beta)

    :param fitspath: str containing folder path
    :param dataset: str containing sub-folder (specific to stacking approach)

    Files identified by default
    composite_file: str containing filename of composite data
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

    # Define [composite_file]
    t_comp = filename_dict['bin_derived_prop'] if not revised else \
        filename_dict['bin_derived_prop']
    composite_file = join(fitspath, dataset, t_comp)
    if not exists(composite_file):
        print("ERROR: File not found! "+composite_file)
        return

    # Read in composite table
    composite_table = asc.read(composite_file)

    bin_id = composite_table['bin_ID'].data
    bin_temp = composite_table['T_e'].data

    # Define [indv_em_line_file]
    indv_em_line_file = join(fitspath, dataset, filename_dict['indv_prop'])
    if not exists(indv_em_line_file):
        print("ERROR: File not found! "+indv_em_line_file)
        return

    # Read in tables containing line ratios, etc.
    indv_em_line_table = asc.read(indv_em_line_file)

    # Define [indv_bin_file]
    indv_bin_file = join(fitspath, dataset, filename_dict['indv_bin_info'])
    if not exists(indv_bin_file):
        print("ERROR: File not found! "+indv_bin_file)
        return

    # Read in tables containing bin info for individual
    indv_bin_info_table = asc.read(indv_bin_file)

    # Populate composite temperature for individual galaxies
    adopted_temp = np.zeros(len(indv_em_line_table))
    bin_id_indv = np.zeros(len(indv_em_line_table))
    for comp_bin, comp_temp in zip(bin_id, bin_temp):
        bin_idx = np.where(indv_bin_info_table['bin_ID'].data == comp_bin)[0]
        adopted_temp[bin_idx] = comp_temp
        bin_id_indv[bin_idx]  = comp_bin

    O2 = indv_em_line_table['OII_3727_Flux_Gaussian'].data            # [OII]3726,3728 fluxes
    O3 = indv_em_line_table['OIII_5007_Flux_Gaussian'].data * OIII_r  # [OIII]4959,5007 fluxes (Assume 3.1:1 ratio)
    Hb = indv_em_line_table['HBETA_Flux_Gaussian'].data               # H-beta fluxes

    if det3:
        com_O_log, metal_dict = metallicity_calculation(adopted_temp, O2/Hb, O3/Hb)
    else:
        det3 = np.where((indv_bin_info_table['Detection'] == 1.0) | (indv_bin_info_table['Detection'] == 0.5))[0]
        temp_com_O_log, temp_metal_dict = \
            metallicity_calculation(adopted_temp[det3], O2[det3]/Hb[det3],
                                    O3[det3]/Hb[det3])
        com_O_log = np.zeros(len(indv_em_line_table))
        com_O_log[det3] = temp_com_O_log

        metal_dict = dict()
        for key0 in temp_metal_dict.keys():
            metal_dict[key0] = np.zeros(len(indv_em_line_table))
            metal_dict[key0][det3] = temp_metal_dict[key0]

    # Define [indv_derived_prop_table] to include ID, bin_ID, composite T_e,
    # and 12+log(O/H)
    arr0 = [indv_em_line_table[ID_name], bin_id_indv, adopted_temp, com_O_log]
    names0 = [ID_name, bin_ID_name] + temp_metal_names0[:2]

    # Include other metallicities
    arr0 += list(metal_dict.values())
    names0 += metal_dict.keys()
    indv_derived_prop_table = Table(arr0, names=names0)

    # Write Astropy ASCII table containing composite T_e and derived metallicity
    if exists(outfile):
        print("File exists! Overwriting : ", outfile)
    else:
        print("Writing : ", outfile)
    indv_derived_prop_table.write(outfile, overwrite=True, format='ascii.fixed_width_two_line')
