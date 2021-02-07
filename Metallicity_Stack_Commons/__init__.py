from typing import Union

from chun_codes.cardelli import cardelli
import astropy.units as u
from datetime import date
import os
import getpass
import numpy as np

from .logging import log_stdout, log_verbose

version = "1.3.1"

lambda0   = [3726.18, 4101.73, 4340.46, 4363.21, 4861.32, 4958.91, 5006.84]
line_type = ['Oxy2', 'Balmer', 'Balmer', 'Single', 'Balmer', 'Single', 'Single']
line_name = ['OII_3727', 'HDELTA', 'HGAMMA', 'OIII_4363', 'HBETA', 'OIII_4958',
             'OIII_5007']

line_name_short = {
    "OII": line_name[0],
    "4363": line_name[3],
    "HB": line_name[4],
    "OIII": line_name[-1],
    "HG": line_name[2],
    "HD": line_name[1]
}

fitting_lines_dict = {
    "lambda0": lambda0,
    "line_type": line_type,
    "line_name": line_name
}

all_lambda0   = [lambda0[0]] + [3728.91] + lambda0[1:]
all_line_name = ['OII_3726', 'OII_3729'] + line_name[1:]
wavelength_dict = dict(zip(all_line_name, all_lambda0))

fitspath_dict = {
    'reagenleimbach': '/Users/reagenleimbach/GoogleDrive/Research/',
    'carol': 'C:/Users/carol/Google Drive/',
    'cly': '/Users/cly/GoogleDrive/Research/',
    'travis': '/home/travis/',
    'runner': '/home/runner/'
}

scalefact = 1e-17

# Flux ratio of [OIII]5007 to [OIII]4959
OIII_r = 3.1

# Define k values for dust attenuation
k_values = cardelli(lambda0 * u.Angstrom)
k_dict   = dict(zip(line_name, k_values))


def exclude_outliers(objno: Union[list, np.ndarray], verbose: bool = False,
                     log: type(log_stdout) = log_stdout()) -> np.ndarray:
    """
    Exclude spectra that are identified as outliers.

    Generally this is because the spectra have very high S/N on the continuum.

    :param objno: Array of eight-digit identifier
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return flag: Array of zeros (not flagged) and ones (flagged
    """

    log_verbose(log, "starting ...", verbose=verbose)

    flag = np.zeros(len(objno), dtype=int)
    bad_data = np.array(['32007727', '32101412', '42006031',
                         '32035286', '14023705'])
    for ii in range(len(bad_data)):
        idx = [xx for xx in range(len(objno)) if
               bad_data[ii] in str(objno[xx])]
        flag[idx] = 1

    log_verbose(log, "finished.", verbose=verbose)
    return flag


def dir_date(folder_name: str, path_init: str = '', year: bool = False,
             verbose: bool = False, log: type(log_stdout) = log_stdout()) \
        -> str:
    """
    This function finds and returns the path to a directory named after the
    current date (MMDDYYYY). If the directory doesn't exist yet, it creates
    a new directory named after the current date in the provided
    ``folder_name`` directory.

    Originally from https://github.com/rafia37/Evolution-of-Galaxies/blob/master/general.py

    Usage:
        fitspath = dir_date(folder_name, year=True)

    :param folder_name: Directory for date subdirectory will be in
    :param path_init: root path. Default: empty string
    :param year: Indicate whether to include year in date folder. Default: False
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return fitspath: Full path to the date directory
    """

    log_verbose(log, "starting ...", verbose=verbose)

    today = date.today()

    list_path = [path_init, folder_name,
                 f"{today.month:02d}{today.day:02d}", '']
    if year:
        list_path[-2] = f"{today.year:d}" + list_path[-2]

    fitspath = os.path.join(*list_path)
    try:
        os.makedirs(fitspath)
    except OSError:
        log.warning(f"Path already exists : {fitspath}")

    log_verbose(log, "finished.", verbose=verbose)
    return fitspath


def get_user(username: Union[None, str] = None,
             verbose: bool = False, log=None) -> str:
    """
    Get the corresponding path for a given ``username``

    :param username: Optional input for username
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return fitspath: Full path to the date directory
    """

    log_verbose(log, "starting ...", verbose=verbose)

    if username is None:
        username = getpass.getuser()

    if username in fitspath_dict.keys():
        fitspath = fitspath_dict[username]
    else:
        log.warning("Incorrect username input")
        raise ValueError("Incorrect username input")

    log_verbose(log, "finished.", verbose=verbose)
    return fitspath
