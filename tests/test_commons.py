from Metallicity_Stack_Commons import column_names, dir_date
from Metallicity_Stack_Commons import get_user, fitspath_reagen, fitspath_caroline

from os.path import exists
from os import rmdir


def test_dir_date():
    mmdd = dir_date('', '')
    mmddyyyy = dir_date('', '', year=True)

    assert exists(mmdd)
    assert exists(mmddyyyy)

    # Check path existence
    mmdd = dir_date('', '')

    if exists(mmdd):
        rmdir(mmdd)

    if exists(mmddyyyy):
        rmdir(mmddyyyy)

    assert len(mmdd) == 5
    assert len(mmddyyyy) == 9


def test_get_user():

    assert get_user('reagenleimbach') == fitspath_reagen
    assert get_user('carol') == fitspath_caroline


def test_column_names():
    # Check python version is integer
    assert isinstance(column_names.py_vers, int)

    # Check type and size for variables in column_names
    assert isinstance(column_names.bin_names0, list)
    assert len(column_names.bin_names0) == 3

    assert isinstance(column_names.indv_names0, list)
    assert len(column_names.indv_names0) == 8

    assert isinstance(column_names.dust0, list)
    assert len(column_names.dust0) == 4

    assert isinstance(column_names.bin_mzevolve_names0, list)
    assert len(column_names.bin_mzevolve_names0) == 8

    assert isinstance(column_names.bin_zcalbase_names0, list)
    assert len(column_names.bin_zcalbase_names0) == 8

    assert isinstance(column_names.bin_ratios0, list)
    assert len(column_names.bin_ratios0) == 5

    assert isinstance(column_names.gauss_names0, list)
    assert len(column_names.gauss_names0) == 8

    assert isinstance(column_names.balmer_names0, list)
    assert len(column_names.balmer_names0) == 2

    assert isinstance(column_names.temp_metal_names0, list)
    assert len(column_names.temp_metal_names0) == 6

    assert isinstance(column_names.valid_table_names0, list)
    assert len(column_names.valid_table_names0) == 5

    assert isinstance(column_names.filename_dict, dict)
    assert len(column_names.filename_dict.keys()) == 13

    assert isinstance(column_names.npz_filename_dict, dict)
    assert len(column_names.npz_filename_dict.keys()) == 6

    column_merge = column_names.merge_column_names(['AAA'], ['BBB'])
    assert isinstance(column_merge, list)
    assert len(column_merge) == 2
