import numpy as np
from os.path import join

from Metallicity_Stack_Commons.analysis import attenuation
from chun_codes import random_pdf


def test_compute_EBV():

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
                EBV = attenuation.compute_EBV(Balmer, source=source, zero_neg=zero)

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
                                              source=source, zero_neg=zero)

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
            Balmer_dist = random_pdf(values, [dx] * len(values), seed_i=1, n_iter=5000)
            EBV_dist = attenuation.compute_EBV(Balmer_dist, source=source, zero_neg=zero)
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
            assert np.array_equal(EBV_dist, npz_reference['EBV_dist'])
