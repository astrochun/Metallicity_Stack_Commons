from os.path import exists, join
from os import symlink
import numpy as np

from astropy.io import fits
from astropy.table import Table

from ..logging import log_stdout
from .. import line_name
from ..column_names import filename_dict


def main(infile: str, log=None):
    """
    Import previous DEEP2 Ly et al. (2015) dataset and export
    an astropy Table called bin_fit.tbl

    :param infile: Input filename
    :param log: LogClass or logging object
    """
    if log is None:
        log = log_stdout()

    log.info(f"Reading : {infile}")
    orig_tab = fits.getdata(infile)

    flux_cols     = [str0 + '_Flux_Gaussian' for str0 in line_name]
    flux_rms_cols = [str0 + '_RMS' for str0 in line_name]

    # Column name in orig_table corresponding to flux_cols
    prefixes = ['OII', 'HD', 'HG', 'OIIIA', 'HB', 'OIIIB', 'OIIIR']
    orig_flux_cols = [str0 + '_FLUX_MOD' for str0 in prefixes]
    orig_snr_cols = [str0 + '_SNR' for str0 in prefixes]

    reformatted_tab = Table()
    reformatted_tab['bin_ID'] = 1 + np.arange(len(orig_tab))

    for ii in range(len(flux_cols)):
        flux = orig_tab[orig_flux_cols[ii]]
        rms  = flux / orig_tab[orig_snr_cols[ii]]
        reformatted_tab[flux_cols[ii]] = flux
        reformatted_tab[flux_rms_cols[ii]] = rms

    out_dir = "tests_data/DEEP2_Ly2015"
    output_table = "DEEP2_Ly2015_fluxes.tbl"
    out_file = join(out_dir, output_table)
    log.info(f"Writing: {out_file}")
    reformatted_tab.write(out_file, format='ascii.fixed_width_two_line')

    symlink_file = f"tests_data/DEEP2_Ly2015/{filename_dict['bin_fit']}"
    if not exists(symlink_file):
        symlink(output_table, symlink_file)
        log.info(f"Symlink created to : {output_table}")
    else:
        log.info(f"Symlink exists: {symlink_file}")
