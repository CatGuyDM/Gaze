"""Regression test: per-charge XIC apex picker must select the parent's apex."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from core.mzml_parser import _build_xic_from_scan_slices, _xic_apex_index

PROTON = 1.00728
C13 = 1.003355

def _make_slices(n_scans, peaks_by_scan):
    slices = []
    for scan_peaks in peaks_by_scan:
        if scan_peaks:
            mzs = np.array([p[0] for p in scan_peaks])
            ints = np.array([p[1] for p in scan_peaks])
            idx = np.argsort(mzs)
            slices.append((mzs[idx], ints[idx]))
        else:
            slices.append((np.array([]), np.array([])))
    return slices

def test_per_charge_apex_selects_parent():
    import math
    parent_mass = 2500.3
    ox_mass = parent_mass + 15.995

    n_scans = 12
    parent_z3_m0 = (parent_mass + 3 * PROTON) / 3
    ox_z2_m0 = (ox_mass + 2 * PROTON) / 2

    def gauss(center, fwhm, i):
        sigma = fwhm / 2.355
        return math.exp(-0.5 * ((i - center) / sigma) ** 2)

    peaks_by_scan = []
    for i in range(n_scans):
        scan_peaks = []
        parent_int = gauss(7, 3.0, i) * 1e6
        if parent_int > 1000:
            scan_peaks.append((parent_z3_m0, parent_int))
            scan_peaks.append((parent_z3_m0 + C13/3, parent_int * 0.9))
            scan_peaks.append((parent_z3_m0 + 2*C13/3, parent_int * 0.65))
        ox_int = gauss(3, 3.0, i) * 1e6
        if ox_int > 1000:
            scan_peaks.append((ox_z2_m0, ox_int))
            scan_peaks.append((ox_z2_m0 + C13/2, ox_int * 0.85))
        peaks_by_scan.append(scan_peaks)

    slices = _make_slices(n_scans, peaks_by_scan)

    target_mzs = [(2, ox_z2_m0), (3, parent_z3_m0)]
    per_charge = {}
    for z, m0_theo in target_mzs:
        xic_z, hits_z = _build_xic_from_scan_slices(slices, [(z, m0_theo)], tolerance_ppm=25.0)
        apex = _xic_apex_index(xic_z, min_fwhm_scans=2, snr_threshold=2.0)
        if apex is not None:
            apex_int = float(xic_z[apex])
            n_iso = max((h[3] for h in hits_z[apex]), default=0) if apex < len(hits_z) and hits_z[apex] else 0
            per_charge[z] = {'apex': apex, 'score': apex_int * max(n_iso, 1)}

    assert 3 in per_charge, "z=3 parent XIC should have a detectable apex"
    best_z = max(per_charge, key=lambda z: per_charge[z]['score'])
    best_apex = per_charge[best_z]['apex']

    assert best_z == 3, f"Expected best charge z=3, got z={best_z}"
    assert 6 <= best_apex <= 8, f"Expected apex near scan 7, got {best_apex}"

if __name__ == '__main__':
    test_per_charge_apex_selects_parent()
    print('PASS: per-charge apex picker selects parent apex')
