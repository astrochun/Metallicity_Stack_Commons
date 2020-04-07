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


def main(fitspath, dataset, revised=False, det3=True):
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

    # Read in validation table
    valid_file = join(fitspath, dataset, filename_dict['bin_valid'])
    if not exists(valid_file):
        print("ERROR: File not found! "+valid_file)
        return
    valid_table = asc.read(valid_file)

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
    detect_indv = np.zeros(len(indv_em_line_table))
    for comp_bin, comp_temp, detect in zip(bin_id, bin_temp, valid_table['Detection']):
        bin_idx = np.where(indv_bin_info_table['bin_ID'].data == comp_bin)[0]
        adopted_temp[bin_idx] = comp_temp
        bin_id_indv[bin_idx]  = comp_bin
        detect_indv[bin_idx]  = detect

    O2 = indv_em_line_table['OII_3727_Flux_Gaussian'].data   # [OII]3726,3728 fluxes
    O3 = indv_em_line_table['OIII_5007_Flux_Gaussian'].data  # [OIII]5007 fluxes
    O3 = O3 * (1+1/OIII_r)  # Scale to include OIII4959; Assume 3.1:1 ratio
    Hb = indv_em_line_table['HBETA_Flux_Gaussian'].data      # H-beta fluxes

    if not det3:
        metal_dict = metallicity_calculation(adopted_temp, O2/Hb, O3/Hb)
    else:
        det3_idx = np.where((detect_indv == 1.0) | (detect_indv == 0.5))[0]
        metal_dict = \
            metallicity_calculation(adopted_temp, O2/Hb, O3/Hb, det3=det3_idx)

    # Define [indv_derived_prop_table] to include ID, bin_ID, composite T_e,
    # and 12+log(O/H)
    arr0 = [indv_em_line_table[ID_name], bin_id_indv, adopted_temp]
    names0 = [ID_name, bin_ID_name, temp_metal_names0[0]]

    # Include other metallicities
    arr0 += list(metal_dict.values())
    names0 += metal_dict.keys()
    indv_derived_prop_table = Table(arr0, names=names0)

    outfile = join(fitspath, dataset, filename_dict['indv_derived_prop'])

    # Write Astropy ASCII table containing composite T_e and derived metallicity
    if exists(outfile):
        print("File exists! Overwriting : ", outfile)
    else:
        print("Writing : ", outfile)
    indv_derived_prop_table.write(outfile, overwrite=True, format='ascii.fixed_width_two_line')
