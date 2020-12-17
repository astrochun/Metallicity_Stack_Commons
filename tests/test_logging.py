from os import remove
from os.path import join

from Metallicity_Stack_Commons import logging


def test_get_user_hostname():

    sys_info = logging.get_user_hostname()
    assert isinstance(sys_info, dict)
    for k in ['user', 'hostname', 'ip', 'os']:
        assert isinstance(sys_info[k], str)


def test_log_stdout():

    log = None
    if log is None:
        log = logging.log_stdout()

        log.info("INFO MESSAGE")
        log.warning("WARNING MESSAGE")
        log.debug("DEBUG MESSAGE")


def test_LogClass():

    log_dir = ''
    logfile = 'testlog.log'

    log = logging.LogClass(log_dir, logfile).get_logger()

    log.info("INFO MESSAGE")
    log.warning("WARNING MESSAGE")
    log.debug("DEBUG MESSAGE")

    remove(join(log_dir, logfile))
