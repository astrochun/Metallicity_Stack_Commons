from astropy.io import ascii as asc
from astropy.table import vstack
from astropy.table import Table
import glob

from ..temp_metallicity_calc import metallicity_calculation


def main(fitspath, dataset, composite_file):
    """
    Purpose:
      Reads in composite table(s) containing bin information to
      determine temperature-based metallicity from composite average
      T_e and individual line ratios ([OII]/H-beta, [OIII]/H-beta)

    :param fitspath: str containing folder path
    :param dataset: str containing sub-folder (specific to stacking approach)
    :param composite_file: str containing filename of composite data

    :return: TBD
    """

    # Read in composite table
    composite_table = asc.read(composite_file)
    ID = composite_table['ID'].data

    det3_table    = asc.read(fitspath + 'get_det3_table2.tbl')
    bin_table     = asc.read(fitspath + dataset + '_2d_binning_datadet3.tbl')
    average_table = asc.read(fitspath + dataset + '_Average_R23_O32_Values.tbl')
    stack_table   = asc.read(fitspath + dataset + '_temperatures_metalicity.tbl')


def run_ind_detection(fitspath, dataset, average_value_ascii):

    N_gal_tab = asc.read(average_value_ascii)
    ID = N_gal_tab['ID']
    for aa in range(len(ID)):
        ind_detection(fitspath, dataset, ID[aa])
    new_name = fitspath + 'Individual_ratio_temperature.tbl'
    vertical_stacking(fitspath, dataset, new_name)
    print('run complete')


def ind_detection(fitspath, dataset, bin_id):
    get_det3_tab = asc.read(fitspath + 'get_det3_table2.tbl')
    bin_tab = asc.read(fitspath + dataset + '_2d_binning_datadet3.tbl')
    N_gal_tab = asc.read(fitspath + dataset + '_Average_R23_O32_Values.tbl')
    stackmeas_tab = asc.read(fitspath + dataset + '_temperatures_metalicity.tbl')

    # From tables
    Source_id = get_det3_tab['Individual_IDs']
    O4959 = get_det3_tab['O4959']
    O5007 = get_det3_tab['O5007']
    Bin_number = bin_tab['Bin_number']
    O2 = get_det3_tab['O2']
    O3 = get_det3_tab['O3']
    Hb = get_det3_tab['Hb']
    N_Galaxies = N_gal_tab['N_Galaxies']
    temp_bin = stackmeas_tab['Temperature']

    R23 = get_det3_tab['R23']
    O32 = get_det3_tab['O32']

    # Initializing Arrays

    Source_IDs = []
    Bin_ID = []
    two_beta = []
    three_beta = []
    OIII4959 = []
    OIII5007 = []
    HBeta = []
    average_temp = []
    R23_ind = []
    O32_ind = []

    for ii in range(len(O2)):
        if Bin_number[ii] == bin_id:
            # print 'Bin_number:', Bin_number[ii], 'O2:', O2[ii], 'O3:', O3[ii], 'Hb:', Hb[ii]
            Source_ID.append(Source_id[ii])
            Bin_ID.append(bin_id)
            R23_ind.append(R23[ii])
            O32_ind.append(O32[ii])
            two_beta.append(O2[ii] / Hb[ii])
            three_beta.append(O3[ii] / Hb[ii])
            OIII4959.append(O4959[ii])
            OIII5007.append(O5007[ii])
            HBeta.append(Hb[ii])
            average_temp.append(temp_bin[bin_id])

    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/' + str(
        bin_id) + '_individual_ratios_temp.tbl'
    n = ('Source_ID', 'Bin_ID', 'Individual_R23', 'Individual_O32', 'two_beta', 'three_beta', 'OIII4959', 'OIII5007',
         'HBeta', 'Temperature')  # 'ID', 'R23_Average', 'O32_Average'
    ind_tab = Table(
        [Source_ID, Bin_ID, R23_ind, O32_ind, two_beta, three_beta, OIII4959, OIII5007, HBeta, average_temp],
        names=n)  # ID, R23, O32,
    asc.write(ind_tab, individual_ascii, format='fixed_width_two_line')


def individual_galaxy_table_stacking(fitspath, dataset, new_name):
    individual_ascii = '/Users/reagenleimbach/Desktop/Zcalbase_gal/individual_detection/*_individual_ratios_temp.tbl'
    table_files = glob.glob(individual_ascii)
    table_files.sort()

    for ii in range(len(table_files)):
        asc_tab = asc.read(table_files[ii])
        if ii == 0:
            vstacking = asc_tab
        else:
            vstacking = vstack([vstacking, asc_tab])
    asc.write(vstacking, new_name, format='fixed_width_two_line', overwrite=True)
