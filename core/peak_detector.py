import numpy as np
from bisect import bisect_left
from dataclasses import dataclass, field
from typing import Optional, List
from core.mzml_parser import Spectrum

ISOTOPE_SPACING = 1.003355


@dataclass
class Peak:
    mz:           float
    intensity:    float
    rt_min:       Optional[float]
    scan_number:  int
    charge_state: Optional[int] = None
    is_monoisotopic: bool = False

    def neutral_mass(self, proton=1.00728):
        if self.charge_state is None or self.charge_state < 1:
            return None
        return self.mz * self.charge_state - self.charge_state * proton


@dataclass
class PeakList:
    scan_number:   int
    rt_min:        Optional[float]
    peaks:         List[Peak]        = field(default_factory=list)
    n_peaks:       int               = 0
    threshold:     float             = 0.0
    base_peak_mz:  Optional[float]  = None
    base_peak_int: float             = 0.0

    def __post_init__(self):
        self.n_peaks = len(self.peaks)
        if self.peaks:
            bp = max(self.peaks, key=lambda p: p.intensity)
            self.base_peak_mz  = bp.mz
            self.base_peak_int = bp.intensity


def _find_local_maxima(int_arr):
    if len(int_arr) < 3:
        return np.array([], dtype=int)
    left  = int_arr[1:-1] > int_arr[:-2]
    right = int_arr[1:-1] > int_arr[2:]
    return np.where(left & right)[0] + 1


def _centroid_peak(mz_arr, int_arr, apex_idx, window=2):
    lo = max(0, apex_idx - window)
    hi = min(len(mz_arr), apex_idx + window + 1)
    mz_w, int_w = mz_arr[lo:hi], int_arr[lo:hi]
    total = int_w.sum()
    if total == 0:
        return float(mz_arr[apex_idx])
    return float(np.dot(mz_w, int_w) / total)


def _is_centroided(mz_arr, sample=200):
    if len(mz_arr) < 10:
        return True
    n = min(len(mz_arr), sample)
    median_gap = float(np.median(np.diff(mz_arr[:n])))
    return median_gap > 0.08


def _expected_m1_over_m0(neutral_mass):
    return max(neutral_mass / 1800.0, 0.05)


def _infer_charge_states(peaks, max_charge=7, n_iso_check=6):
    if len(peaks) < 2:
        return peaks

    sort_idx = sorted(range(len(peaks)), key=lambda i: peaks[i].mz)
    mz_arr   = [peaks[i].mz for i in sort_idx]
    int_arr  = [peaks[i].intensity for i in sort_idx]
    sorted_peaks = [peaks[i] for i in sort_idx]

    spacings = {z: ISOTOPE_SPACING / z for z in range(1, max_charge + 1)}
    n_peaks  = len(mz_arr)
    bisect_left_local = bisect_left

    for i, peak in enumerate(sorted_peaks):
        if peak.charge_state is not None:
            continue

        m0_int = peak.intensity
        candidates = []

        for z, spacing in spacings.items():
            score = 0
            iso_idxs = []
            m1_intensity = None
            for k in range(1, n_iso_check + 1):
                target = peak.mz + k * spacing
                tol_da = 0.5 * spacing
                ins = bisect_left_local(mz_arr, target)
                hit = ins
                if hit >= n_peaks:
                    hit = n_peaks - 1
                if hit > 0 and (hit == n_peaks or
                                 abs(mz_arr[hit - 1] - target) < abs(mz_arr[hit] - target)):
                    hit = hit - 1
                if abs(mz_arr[hit] - target) < tol_da and hit > i:
                    iso_ratio = int_arr[hit] / max(m0_int, 1.0)
                    if 0.05 < iso_ratio < 4.0:
                        score += 1
                        iso_idxs.append(hit)
                        if k == 1:
                            m1_intensity = int_arr[hit]
                        continue
                break

            if score >= 2:
                neutral_mass = peak.mz * z - z * 1.00728
                expected_r1  = _expected_m1_over_m0(neutral_mass)
                obs_r1       = (m1_intensity / max(m0_int, 1.0)) if m1_intensity else 0.0
                fit_err = abs(obs_r1 - expected_r1) / max(expected_r1, 0.05)
                candidates.append((score, fit_err, z, iso_idxs))

        if not candidates:
            continue

        candidates.sort(key=lambda c: (-c[0], c[1], -c[2]))
        score, fit_err, best_z, best_iso_idxs = candidates[0]

        peak.charge_state = best_z
        peak.is_monoisotopic = True
        for j in best_iso_idxs:
            if sorted_peaks[j].charge_state is None:
                sorted_peaks[j].charge_state = best_z
                sorted_peaks[j].is_monoisotopic = False

    return peaks


def detect_peaks(spectrum, threshold_fraction=0.01, centroid_window=2,
                 infer_charges=True, min_peaks=0, triage_mode=False):
    if triage_mode:
        threshold_fraction = 0.0005

    if spectrum.is_empty or len(spectrum.mz_array) == 0:
        return PeakList(scan_number=spectrum.scan_number, rt_min=spectrum.rt_min, peaks=[])

    mz_arr, int_arr = spectrum.mz_array, spectrum.intensity_array
    base_int  = float(int_arr.max())
    threshold = base_int * threshold_fraction

    is_cent = getattr(spectrum, 'centroided', None)
    if is_cent is None:
        is_cent = _is_centroided(mz_arr)

    if is_cent:
        keep  = int_arr >= threshold
        peaks = [
            Peak(mz=round(float(mz_arr[i]), 5), intensity=float(int_arr[i]),
                 rt_min=spectrum.rt_min, scan_number=spectrum.scan_number)
            for i in np.where(keep)[0]
        ]
    else:
        maxima = _find_local_maxima(int_arr)
        maxima = maxima[int_arr[maxima] >= threshold]
        if len(maxima) == 0:
            return PeakList(scan_number=spectrum.scan_number, rt_min=spectrum.rt_min,
                            peaks=[], threshold=threshold)
        peaks = []
        for idx in maxima:
            centroid_mz = _centroid_peak(mz_arr, int_arr, idx, centroid_window)
            peaks.append(Peak(mz=round(centroid_mz, 5), intensity=float(int_arr[idx]),
                              rt_min=spectrum.rt_min, scan_number=spectrum.scan_number))

    peaks.sort(key=lambda p: p.mz)

    if infer_charges and len(peaks) >= 2:
        peaks = _infer_charge_states(peaks)

    if len(peaks) < min_peaks:
        return PeakList(scan_number=spectrum.scan_number, rt_min=spectrum.rt_min,
                        peaks=[], threshold=threshold)

    return PeakList(scan_number=spectrum.scan_number, rt_min=spectrum.rt_min,
                    peaks=peaks, threshold=threshold)
