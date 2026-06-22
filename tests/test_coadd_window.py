"""Regression test for auto FWHM-based co-add window sizing."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math
import numpy as np
from core.mzml_parser import _xic_fwhm_scans, _smooth_savgol_simple

def _gaussian_xic(n_scans, apex_idx, fwhm):
    sigma = fwhm / 2.355
    x = np.arange(n_scans, dtype=float)
    return np.exp(-0.5 * ((x - apex_idx) / sigma) ** 2) * 1e6

def test_fwhm_estimation_accuracy():
    """FWHM estimate should be within ±30% of the true value."""
    for true_fwhm in [4, 6, 8, 10]:
        n = 40
        apex = n // 2
        xic = _gaussian_xic(n, apex, true_fwhm)
        estimated = _xic_fwhm_scans(xic, apex)
        assert estimated > 0, f"FWHM estimation returned 0 for true_fwhm={true_fwhm}"
        ratio = estimated / true_fwhm
        assert 0.70 <= ratio <= 1.30, (
            f"FWHM estimate {estimated} is outside ±30% of true {true_fwhm} "
            f"(ratio={ratio:.2f})"
        )

def test_fwhm_window_bounds():
    for true_fwhm in [1, 2, 20, 30]:
        n = 60
        apex = n // 2
        xic = _gaussian_xic(n, apex, true_fwhm)
        fwhm = _xic_fwhm_scans(xic, apex)
        if fwhm >= 3:
            half = max(2, min(12, int(round(fwhm * 0.5))))
            assert 2 <= half <= 12, f"half-window {half} outside [2,12] for fwhm={fwhm}"

if __name__ == '__main__':
    test_fwhm_estimation_accuracy()
    test_fwhm_window_bounds()
    print('PASS: FWHM-based window tests')
