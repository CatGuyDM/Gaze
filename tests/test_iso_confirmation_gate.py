"""HIGH-risk impurity at >=5% without isotope confirmation must downgrade to MEDIUM."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from core.peak_detector import PeakList, Peak
from core.peak_matcher import match_peaks

PROTON = 1.00728
SEQUENCE = 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'

def test_high_intensity_unconfirmed_downgraded():
    from core.mass_calculator import calculate_monoisotopic_mass
    parent_mass = calculate_monoisotopic_mass(SEQUENCE, c_terminal_amide=True)

    peaks = []
    base_int = 1_000_000.0
    m0_z4 = (parent_mass + 4 * PROTON) / 4
    C13 = 1.003355
    for k in range(6):
        ratio = [1.0, 0.9, 0.65, 0.38, 0.19, 0.08][k]
        peaks.append(Peak(mz=round(m0_z4 + k * C13 / 4, 5),
                          intensity=base_int * ratio,
                          rt_min=5.0, scan_number=100))

    del_mass = parent_mass - 71.03711
    del_m0_z4 = (del_mass + 4 * PROTON) / 4
    peaks.append(Peak(mz=round(del_m0_z4, 5), intensity=base_int * 0.10,
                      rt_min=5.0, scan_number=100))

    pl = PeakList(scan_number=100, rt_min=5.0, peaks=peaks)
    mr = match_peaks(
        pl, SEQUENCE,
        tolerance_ppm=15.0,
        c_terminal_amide=True,
        coupling_reagent='HATU',
        min_envelope_score=0.0,
        min_confidence=0.0,
    )

    high_unconfirmed = [
        m for m in mr.matches
        if m.risk_level == 'HIGH' and not m.isotope_confirmed
        and 100.0 * m.observed_intensity / base_int >= 5.0
    ]
    assert len(high_unconfirmed) == 0, (
        f"Expected 0 HIGH-risk unconfirmed matches at ≥5%, got {len(high_unconfirmed)}: "
        + str([(m.impurity_name, m.risk_level, m.iso_warning) for m in high_unconfirmed])
    )

    del_matches = [m for m in mr.matches
                   if 'del' in m.impurity_type.lower()
                   and 100.0 * m.observed_intensity / base_int >= 5.0]
    for dm in del_matches:
        assert dm.risk_level in ('MEDIUM', 'LOW'), (
            f"Expected MEDIUM/LOW for unconfirmed high-intensity deletion, got {dm.risk_level}")
        assert dm.iso_warning, "Expected iso_warning=True for unconfirmed high-intensity match"

if __name__ == '__main__':
    test_high_intensity_unconfirmed_downgraded()
    print('PASS: high-intensity unconfirmed match downgraded to MEDIUM')
