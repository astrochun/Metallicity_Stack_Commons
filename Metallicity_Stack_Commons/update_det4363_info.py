from typing import Tuple
import numpy as np

from chun_codes import match_nosort

from astropy.table import Table

from .logging import log_stdout, log_verbose


def get_index(det4363_table: Table, input_table: Table, column_name: str,
              verbose: bool = False, log: type(log_stdout) = log_stdout()) -> \
        Tuple[np.ndarray, np.ndarray]:
    """
    Uses either OBJNO or AP/SLIT info to get index for an existing table

    :param det4363_table: Astropy table containing DEEP2 [OIII]4363-detected sample
    :param input_table: Astropy table containing the entire sample to be updated
    :param column_name: Column name for cross-matching
    :param verbose: Write verbose message to stdout. Default: file only
    :param log: LogClass or logging object

    :return: Index arrays for ``det4363_table``, ``input_table``
    """

    log_verbose(log, "starting ...", verbose=verbose)

    if column_name != 'OBJNO' or column_name != 'SLIT' or column_name != 'AP':
        log.warning("column_name not understood")
        log.warning("Exiting!!!")
        raise KeyError("column_name not understood")

    det4363_id = det4363_table[column_name]
    input_id   = input_table[column_name]

    det4363_idx, input_idx = match_nosort(det4363_id, input_id)

    log_verbose(log, "finished.", verbose=verbose)
    return det4363_idx, input_idx
