from Metallicity_Stack_Commons.analysis import attenuation

import numpy as np


def test_compute_EBV():

    offsets = [0, -0.01, 0.01]

    for value, source in zip([attenuation.HgHb_CaseB, attenuation.HdHb_CaseB],
                             ['HgHb', 'HdHb']):
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
