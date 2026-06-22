"""Regression: a pure-parent spectrum must not produce a deamidation hit."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from core.peak_detector import PeakList, Peak
from core.peak_matcher import match_peaks, theoretical_isotope_envelope

SEQUENCE = 'REQLKHKLEQLRNSCA'
PROTON = 1.00728
C13 = 1.003355

def _make_synthetic_peaks(parent_mass, charges=(1, 2, 3)):
    theo = theoretical_isotope_envelope(parent_mass, n_iso=5)
    peaks = []
    base_int = 1_000_000.0
    for z in charges:
        m0 = (parent_mass + z * PROTON) / z
        for k, ratio in enumerate(theo[:6]):
            mz = m0 + k * C13 / z
            intensity = base_int * ratio / z
            if intensity > 1000:
                peaks.append(Peak(mz=round(mz, 5), intensity=intensity,
                                  rt_min=5.0, scan_number=258))
    return peaks

def test_no_deamidation_on_pure_parent():
    from core.mass_calculator import calculate_monoisotopic_mass
    parent_mass = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)
    peaks = _make_synthetic_peaks(parent_mass)
    pl = PeakList(scan_number=258, rt_min=5.0, peaks=peaks)
    mr = match_peaks(
        pl, SEQUENCE,
        tolerance_ppm=15.0,
        c_terminal_amide=True,
        coupling_reagent='HATU',
    )
    deam_matches = [m for m in mr.matches if m.impurity_type == 'deamidation']
    assert len(deam_matches) == 0, (
        f"Expected 0 deamidation matches on pure-parent spectrum, got {len(deam_matches)}: "
        + str([(m.charge_state, m.observed_mz, m.envelope_score) for m in deam_matches])
    )

if __name__ == '__main__':
    test_no_deamidation_on_pure_parent()
    print('PASS: no deamidation on pure-parent spectrum')
