import numpy as np
from os.path import join

from Metallicity_Stack_Commons.analysis import attenuation
from Metallicity_Stack_Commons import line_name_short
from chun_codes import random_pdf

import pytest


def test_compute_EBV():

    # Test Error handling when list is provided
    with pytest.raises(TypeError):
        attenuation.compute_EBV([0.45, 0.38], verbose=True)

    # Test Error handling when source is incorrect
    with pytest.raises(KeyError):
        attenuation.compute_EBV(0.45, source='test', verbose=True)

    dx = 0.05  # This is for randomization

    offsets = [0, -0.01, 0.01]

    Balmer_list = [attenuation.HgHb_CaseB, attenuation.HdHb_CaseB]
    source_list = ['HgHb', 'HdHb']
    zip_data = zip(Balmer_list, source_list)

    for value, source in zip_data:
        for zero in [False, True]:
            for offset in offsets:
                # Test float input
                Balmer = value + offset
                EBV = attenuation.compute_EBV(Balmer, source=source,
                                              zero_neg=zero, verbose=True)

                assert isinstance(EBV, float)
                if offset == 0:
                    assert EBV == 0
                if offset < 0:
                    assert EBV > 0

                if offset > 0:
                    if not zero:
                        assert EBV < 0
                    else:
                        assert EBV == 0

                # Test numpy array with single record
                EBV = attenuation.compute_EBV(np.array([Balmer]),
                                              source=source, zero_neg=zero,
                                              verbose=True)

                assert isinstance(EBV, (np.ndarray, np.generic))
                if offset == 0:
                    assert EBV == 0
                if offset < 0:
                    assert EBV > 0

                if offset > 0:
                    if not zero:
                        assert EBV < 0
                    else:
                        assert EBV == 0

            # Test EBV distribution case
            values = [value, value - dx, value + dx]
            Balmer_dist = random_pdf(values, [dx] * len(values), seed_i=1,
                                     n_iter=5000)
            if not zero:
                EBV_dist = attenuation.compute_EBV(Balmer_dist, source=source,
                                                   zero_neg=zero, verbose=True)
            else:
                EBV_dist, EBV_peak = attenuation.compute_EBV(Balmer_dist,
                                                             source=source,
                                                             zero_neg=zero,
                                                             verbose=True)

            '''
            # For writing initial file
            npz_outfile = join('tests_data', f'EBV_dist_{source}_{zero}.npz')
            print(f"Writing : {npz_outfile}")
            np.savez(npz_outfile, EBV_dist=EBV_dist)
            '''
            assert isinstance(EBV_dist, (np.ndarray, np.generic))

            # Read in reference data
            npz_infile = join('tests_data', f'EBV_dist_{source}_{zero}.npz')
            print(f"Reading : {npz_infile}")
            npz_reference = np.load(npz_infile)
            # assert np.array_equal(EBV_dist, npz_reference['EBV_dist'])


def test_compute_A():

    for EBV in [0.0, 0.25]:
        A_dict = attenuation.compute_A(EBV, verbose=True)
        assert isinstance(A_dict, dict)
        for key in A_dict.keys():
            if EBV == 0:
                assert A_dict[key] == 0.0
            else:
                assert A_dict[key] > 0


def test_line_ratio_atten():

    ratio = 2.0

    for EBV in [0.0, 0.25]:
        # [OII]/H-beta
        ratio_atten = attenuation.line_ratio_atten(ratio, EBV,
                                                   line_name_short['OII'],
                                                   line_name_short['HB'],
                                                   verbose=True)
        assert isinstance(ratio_atten, float)
        if EBV == 0:
            assert ratio_atten == ratio
        else:
            assert ratio_atten > ratio

        # [OIII]/[OII]
        ratio_atten = attenuation.line_ratio_atten(ratio, EBV,
                                                   line_name_short['OIII'],
                                                   line_name_short['OII'],
                                                   verbose=True)
        assert isinstance(ratio_atten, float)
        if EBV == 0:
            assert ratio_atten == ratio
        else:
            assert ratio_atten < ratio


def test_Hb_SFR():

    logSFR = attenuation.Hb_SFR(41.0, 0.25, verbose=True)
    assert isinstance(logSFR, float)
