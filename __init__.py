from chun_codes.cardelli import cardelli
import astropy.units as u
from datetime import date
import os

lambda0   = [3726.18, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']
line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958', 'OIII_5007']

fitting_lines_dict = {"lambda0":lambda0, "line_type":line_type, "line_name":line_name}

fitspath_reagen = '/Users/reagenleimbach/Desktop/Zcalbase_gal/'

fitspath_caroline = 'C:\\Users\\carol\\Google Drive\\'

scalefact = 1e-17

# Define k values for dust attenuation
k_values = cardelli(lambda0 * u.Angstrom)
k_dict   = dict(zip(line_name,k_values))


def exclude_outliers(objno):
    '''
    Exclude spectra that are identified as outliers.

    Generally this is because the spectra have very high S/N on the continuum.

    :param objno: list or numpy array of eight-digit identifier

    :return:
     flag: numpy array of zeros and ones
    '''

    flag = np.zeros(len(objno), dtype=int)
    bad_data = np.array(['32007727', '32101412', '42006031', '32035286', '14023705'])
    for ii in range(len(bad_data)):
        idx = [xx for xx in range(len(objno)) if bad_data[ii] == str(objno[xx])][0]
        flag[idx] = 1

    return flag


def dir_date(org_name, path_init='', year=False):
    '''
    Purpose:
        This function finds and returns the path to a directory named after the
        current date (MMDDYYYY). If the directory doesn't exist yet, it creates
        a new directory named after the current date in the provided org_name
        directory.

        From https://github.com/rafia37/Evolution-of-Galaxies/blob/master/general.py
    Usage:
        fitspath = general.get_time(org_name, path_init='', year=True)

    Params:
        org_name --> a string of the directory that the date subdirectory will be in.

    Returns:
        fitspath --> the path to the date directory.

    Outputs:
        "Path already exists" --> prints this message if the current date directory already exists.
        fitspath --> prints the path to the directory.

    '''

    today = date.today()

    list_path = [path_init, org_name, "%02i%02i" % (today.month, today.day), '']
    if year:
        list_path[-2] += "i%02i%02i" % today.year

    fitspath = os.path.join(*list_path)
    try:
        os.mkdir(fitspath)
    except FileExistsError:
        print("Path already exists : ", fitspath)

    return fitspath