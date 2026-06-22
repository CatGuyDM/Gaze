import math
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List

from core.mass_calculator import (
    calculate_monoisotopic_mass, calculate_mz, calculate_ppm_error,
    NON_STANDARD_RESIDUES, NON_STANDARD_RESIDUES_AVG,
    peptide_isotope_envelope,
)
from core.peak_detector import PeakList, Peak

PROTON = 1.00728
_C13_ABUND = 0.01103
_N15_ABUND = 0.00368
_O18_ABUND = 0.00205
_S34_ABUND = 0.04253


def theoretical_isotope_ratios(neutral_mass):
    n = neutral_mass / 111.1
    r1 = n * (4.9384 * _C13_ABUND + 1.3577 * _N15_ABUND + 1.4773 * _O18_ABUND + 0.0417 * _S34_ABUND)
    return float(r1), float(r1 ** 2 / 2.0)


def theoretical_isotope_envelope(neutral_mass, n_iso=5):
    r1, r2 = theoretical_isotope_ratios(neutral_mass)
    n_c = neutral_mass / 111.1 * 4.9384
    lam = n_c * _C13_ABUND

    out = []
    for k in range(0, n_iso + 1):
        out.append(math.exp(-lam) * (lam ** k) / math.factorial(k))
    if out[0] > 0:
        out = [v / out[0] for v in out]
    return out


def envelope_fit_score(peaks, m0_mz, charge, neutral_mass,
                       tolerance_ppm=15.0, n_iso=4):
    if neutral_mass < 800.0:
        return 0.5, {'reason': 'mass below envelope reliability threshold'}

    spacing = 1.00335 / charge
    m0_peak = None
    for p in peaks:
        if calculate_ppm_error(m0_mz, p.mz) < tolerance_ppm:
            if m0_peak is None or p.intensity > m0_peak.intensity:
                m0_peak = p
    if m0_peak is None or m0_peak.intensity <= 0:
        return 0.0, {'reason': 'M0 not found'}

    theo = theoretical_isotope_envelope(neutral_mass, n_iso=n_iso)

    obs_intensities = [m0_peak.intensity]
    for k in range(1, n_iso + 1):
        target = m0_mz + k * spacing
        peak_k = None
        for p in peaks:
            if calculate_ppm_error(target, p.mz) < tolerance_ppm:
                if peak_k is None or p.intensity > peak_k.intensity:
                    peak_k = p
        obs_intensities.append(peak_k.intensity if peak_k else 0.0)

    obs_ratios = [v / m0_peak.intensity for v in obs_intensities]

    score_terms  = []
    weight_sum   = 0.0
    weighted_sum = 0.0
    for k in range(1, n_iso + 1):
        theo_k = theo[k]
        obs_k  = obs_ratios[k]
        if theo_k < 0.02:
            continue
        w = theo_k
        if obs_k == 0.0:
            rel_err = 1.0
        else:
            rel_err = min(1.0, abs(obs_k - theo_k) / theo_k)
        contribution = (1.0 - rel_err) * w
        score_terms.append({
            'k': k, 'theo_ratio': theo_k, 'obs_ratio': obs_k, 'rel_err': rel_err,
        })
        weight_sum   += w
        weighted_sum += contribution

    if weight_sum < 1e-6:
        return 0.5, {'reason': 'no informative theoretical isotopes',
                     'theoretical': theo}
    score = weighted_sum / weight_sum

    if charge >= 2 and obs_ratios[1] == 0.0 and theo[1] > 0.1:
        score *= 0.3

    return float(round(score, 3)), {
        'theoretical': [round(t, 3) for t in theo],
        'observed': [round(v, 3) for v in obs_ratios],
        'per_iso': score_terms,
    }


def check_isotope_pattern(peaks, apex_mz, charge, neutral_mass,
                          tolerance_ppm=10.0, ratio_tolerance=0.30):
    if neutral_mass < 800.0:
        return False
    score, details = envelope_fit_score(
        peaks, apex_mz, charge, neutral_mass,
        tolerance_ppm=tolerance_ppm, n_iso=3,
    )
    if score < 0.5:
        return False
    if 'observed' in details and 'theoretical' in details:
        obs1 = details['observed'][1] if len(details['observed']) > 1 else 0
        the1 = details['theoretical'][1] if len(details['theoretical']) > 1 else 0
        if the1 > 0 and abs(obs1 - the1) / the1 > ratio_tolerance:
            return False
    return True


@dataclass
class PeakMatch:
    observed_mz:        float
    observed_intensity: float
    theoretical_mz:     float
    ppm_error:          float
    charge_state:       int
    impurity_name:      str
    impurity_type:      str
    risk_level:         str
    confidence:         float
    isotope_confirmed:  bool
    delta_mono:         float
    predicted_mass:     float
    notes:              str = ''
    envelope_score:     float = 0.0
    position:           str   = ''
    n_degenerate:       int   = 1
    scan_persistence:   float = 1.0
    specificity:        float = 1.0
    triage_tier:        str   = ''
    iso_warning:        bool  = False


@dataclass
class MatchResult:
    sequence:               str
    parent_mass:            float
    scan_number:            int
    n_peaks_input:          int
    n_peaks_matched:        int
    n_impurities_predicted: int
    matches:                List[PeakMatch] = field(default_factory=list)
    unmatched_peaks:        List[float]     = field(default_factory=list)
    unmatched_impurities:   List[str]       = field(default_factory=list)
    bear_case_applied:      bool            = False
    tolerance_ppm:          float           = 10.0
    notes:                  str             = ''

    total_intensity:                   float = 0.0
    parent_envelope_intensity:         float = 0.0
    matched_impurity_intensity:        float = 0.0
    unexplained_intensity:             float = 0.0
    parent_envelope_fraction:          float = 0.0
    matched_impurity_fraction:         float = 0.0
    unexplained_fraction:              float = 0.0
    n_unexplained_peaks_above_1pct:    int   = 0
    top_unexplained_peaks:             list  = field(default_factory=list)

    native_parent_intensity:           float = 0.0
    deamidated_parent_intensity:       float = 0.0
    native_parent_fraction:            float = 0.0
    deamidated_parent_fraction:        float = 0.0
    iemm_applied:                      bool  = False

    @property
    def unique_impurity_types_matched(self):
        return len({m.impurity_type for m in self.matches
                    if m.confidence >= 4.0 and m.impurity_type not in ('racemization',)})

    @property
    def sensitivity(self):
        if self.n_impurities_predicted == 0:
            return 1.0
        _EXCLUDE = {'racemization'}
        matched_types   = {m.impurity_type for m in self.matches
                           if m.impurity_type not in _EXCLUDE}
        unmatched_types = {t for t in self.unmatched_impurities
                           if t not in _EXCLUDE and t not in matched_types}
        total_types = len(matched_types) + len(unmatched_types)
        if total_types == 0:
            return 1.0
        return min(1.0, len(matched_types) / total_types)

    @property
    def precision(self):
        if not self.matches:
            return 1.0
        confirmed = sum(1 for m in self.matches if m.confidence >= 4.0)
        return confirmed / max(len(self.matches), 1)

    def to_dataframe(self):
        if not self.matches:
            return pd.DataFrame()
        return pd.DataFrame([{
            'sequence': self.sequence, 'impurity_name': m.impurity_name,
            'impurity_type': m.impurity_type, 'risk_level': m.risk_level,
            'observed_mz': m.observed_mz, 'theoretical_mz': m.theoretical_mz,
            'ppm_error': m.ppm_error, 'charge_state': m.charge_state,
            'confidence': m.confidence, 'isotope_confirmed': m.isotope_confirmed,
            'delta_mono': m.delta_mono, 'predicted_mass': m.predicted_mass,
        } for m in self.matches])


def _score(ppm_error, risk_level, intensity, base_peak_intensity, isotope_confirmed):
    base        = {'HIGH': 5.0, 'MEDIUM': 4.0, 'LOW': 3.0}.get(risk_level, 3.0)
    iso_bonus   = 1.5 if isotope_confirmed else 0.0

    rel_int = intensity / max(base_peak_intensity, 1.0)
    int_bonus = 1.5 if rel_int > 0.05 else (1.0 if rel_int > 0.01 else (0.5 if rel_int > 0.003 else 0.0))
    ppm_bonus   = 1.0 if ppm_error < 2.0 else (0.5 if ppm_error < 4.0 else 0.0)
    ppm_malus   = -1.0 if ppm_error > 8.0 else 0.0
    return round(max(0.0, min(10.0, base + iso_bonus + int_bonus + ppm_bonus + ppm_malus)), 2)


def _calibrate_from_anchors(anchors, parent_mass):
    if not anchors:
        return lambda mz: mz

    by_charge = {}
    for hit in anchors:
        if len(hit) < 3 or hit[0] is None:
            continue
        z, obs_m0, intensity = hit[0], hit[1], hit[2]
        n_iso = hit[3] if len(hit) > 3 else 2
        by_charge.setdefault(z, []).append((obs_m0, intensity, n_iso))

    collapsed = []
    for z, hits in by_charge.items():
        obs_arr  = np.array([h[0] for h in hits])
        wgt_arr  = np.array([h[1] * h[2] for h in hits])
        if wgt_arr.sum() <= 0:
            continue
        obs_med = float(np.average(obs_arr, weights=wgt_arr))
        collapsed.append((z, obs_med, float(wgt_arr.sum())))

    if len(collapsed) < 2:
        if len(collapsed) == 1:
            z, obs, _ = collapsed[0]
            theo = (parent_mass + z * 1.00728) / z
            scale = obs / theo
            return lambda mz, _s=scale: mz / _s
        return lambda mz: mz

    obs_arr  = np.array([o for _, o, _ in collapsed])
    theo_arr = np.array([(parent_mass + z * 1.00728) / z for z, _, _ in collapsed])
    wgt_arr  = np.array([w for _, _, w in collapsed])

    err_ppm = (obs_arr - theo_arr) / theo_arr * 1e6

    mx = np.average(obs_arr, weights=wgt_arr)
    my = np.average(err_ppm, weights=wgt_arr)
    cov = np.average((obs_arr - mx) * (err_ppm - my), weights=wgt_arr)
    var = np.average((obs_arr - mx) ** 2, weights=wgt_arr)
    if var < 1e-9:
        offset_ppm = my
        return lambda mz, _o=offset_ppm: mz - _o * mz * 1e-6
    slope     = cov / var
    intercept = my - slope * mx

    if len(collapsed) >= 3:
        residuals = err_ppm - (slope * obs_arr + intercept)
        spread = float(np.sqrt(np.average(residuals ** 2, weights=wgt_arr)))
        if spread > 0:
            keep = np.abs(residuals) <= max(3 * spread, 5.0)
            if keep.sum() >= 2 and keep.sum() < len(collapsed):
                obs_k = obs_arr[keep]; err_k = err_ppm[keep]; w_k = wgt_arr[keep]
                mx2 = np.average(obs_k, weights=w_k)
                my2 = np.average(err_k, weights=w_k)
                cov2 = np.average((obs_k - mx2) * (err_k - my2), weights=w_k)
                var2 = np.average((obs_k - mx2) ** 2, weights=w_k)
                if var2 > 1e-9:
                    slope     = cov2 / var2
                    intercept = my2 - slope * mx2

    def correct(mz, _a=slope, _b=intercept):
        err_ppm_at_mz = _a * mz + _b
        return mz - err_ppm_at_mz * mz * 1e-6

    return correct


def _iemm_deconvolve_charge(observed_intensities, theoretical_native,
                             theoretical_deamidated):
    O = np.asarray(observed_intensities, dtype=float)
    Tn = np.asarray(theoretical_native, dtype=float)
    Td = np.asarray(theoretical_deamidated, dtype=float)
    if len(O) == 0 or len(Tn) == 0 or len(Td) == 0:
        return 0.0, 0.0
    n = min(len(O), len(Tn), len(Td))
    O, Tn, Td = O[:n], Tn[:n], Td[:n]

    A11 = float(np.sum(Tn * Tn))
    A12 = float(np.sum(Tn * Td))
    A22 = float(np.sum(Td * Td))
    B1  = float(np.sum(O * Tn))
    B2  = float(np.sum(O * Td))
    det = A11 * A22 - A12 * A12

    if det <= 0 or A11 <= 0:
        if A11 > 0:
            return max(0.0, B1 / A11), 0.0
        return 0.0, 0.0

    a = (A22 * B1 - A12 * B2) / det
    b = (A11 * B2 - A12 * B1) / det

    if a < 0 and b < 0:
        return 0.0, 0.0
    if a < 0:
        b = max(0.0, B2 / A22) if A22 > 0 else 0.0
        return 0.0, b
    if b < 0:
        a = max(0.0, B1 / A11) if A11 > 0 else 0.0
        return a, 0.0
    return a, b


def _self_calibrate(peak_list, sequence, search_ppm=300,
                    c_terminal_amide=False, fixed_modification_da=0.0,
                    extra_residues=None, disulfide_bonds=0):
    from core.mass_calculator import calculate_mz as _calc_mz, calculate_monoisotopic_mass

    try:
        parent_mass = calculate_monoisotopic_mass(
            sequence, c_terminal_amide=c_terminal_amide,
            disulfide_bonds=disulfide_bonds,
            fixed_modification_da=fixed_modification_da,
            extra_residues=extra_residues,
        )
    except Exception:
        return lambda mz: mz

    C13    = 1.003355
    peak_mz  = np.array([p.mz       for p in peak_list.peaks])
    peak_int = np.array([p.intensity for p in peak_list.peaks])

    anchors = []
    for z in range(2, 8):
        theo_m0 = _calc_mz(parent_mass, z)
        spacing = C13 / z

        z_anchors = []
        for k in range(0, 10):
            theo_mk = theo_m0 + k * spacing
            window_mk = theo_mk * search_ppm * 1e-6
            cands = np.where(np.abs(peak_mz - theo_mk) < window_mk)[0]
            if len(cands) == 0:
                continue
            best  = cands[np.argmin(np.abs(peak_mz[cands] - theo_mk))]
            obs_mk = float(peak_mz[best])
            err = (obs_mk - theo_mk) / theo_mk * 1e6
            if abs(err) > search_ppm:
                continue

            if z_anchors:
                prev_obs_mk = z_anchors[-1][1]
                prev_k      = z_anchors[-1][3]
                obs_gap     = obs_mk - prev_obs_mk
                k_gap       = k - prev_k
                expected    = k_gap * spacing
                if abs(obs_gap - expected) > expected * 0.10:
                    continue

            z_anchors.append((theo_mk, obs_mk, float(peak_int[best]), k))

        if len(z_anchors) < 2:
            continue

        for theo_mk, obs_mk, w, k in z_anchors:
            anchors.append((theo_mk, obs_mk, w))

    if len(anchors) == 0:
        return lambda mz: mz

    if len(anchors) == 1:
        theo, obs, _ = anchors[0]
        scale = obs / theo
        return lambda mz, _s=scale: mz / _s

    def _wlfit(obs_arr, err_arr, weights):
        mx  = np.average(obs_arr, weights=weights)
        my  = np.average(err_arr, weights=weights)
        cov = np.average((obs_arr - mx) * (err_arr - my), weights=weights)
        var = np.average((obs_arr - mx) ** 2,             weights=weights)
        if var < 1e-9:
            return 0.0, my
        slope     = cov / var
        intercept = my - slope * mx
        return slope, intercept

    obs_arr = np.array([a[1] for a in anchors])
    err_ppm = np.array([(a[1] - a[0]) / a[0] * 1e6 for a in anchors])
    weights = np.array([a[2] for a in anchors])

    slope, intercept = _wlfit(obs_arr, err_ppm, weights)

    if len(anchors) >= 3:
        residuals = err_ppm - (slope * obs_arr + intercept)
        spread = np.sqrt(np.average(residuals ** 2, weights=weights))
        if spread > 0:
            keep = np.abs(residuals) <= max(3 * spread, 5.0)
            if keep.sum() >= 2 and keep.sum() < len(anchors):
                slope, intercept = _wlfit(obs_arr[keep], err_ppm[keep], weights[keep])

    def correct(mz, _a=slope, _b=intercept):
        err_ppm = _a * mz + _b
        return mz - err_ppm * mz * 1e-6

    return correct


def _two_stage_calibrate(peak_list, sequence,
                         c_terminal_amide=False, fixed_modification_da=0.0,
                         extra_residues=None, disulfide_bonds=0,
                         stage1_search_ppm=1500, stage2_search_ppm=200):
    stage1 = _self_calibrate(
        peak_list, sequence,
        search_ppm=stage1_search_ppm,
        c_terminal_amide=c_terminal_amide,
        fixed_modification_da=fixed_modification_da,
        extra_residues=extra_residues,
        disulfide_bonds=disulfide_bonds,
    )

    test_mz = 1000.0
    if abs(stage1(test_mz) - test_mz) < 1e-6:
        return stage1

    corrected_peaks = [
        Peak(mz=stage1(p.mz), intensity=p.intensity,
             rt_min=p.rt_min, scan_number=p.scan_number,
             charge_state=getattr(p, 'charge_state', None),
             is_monoisotopic=getattr(p, 'is_monoisotopic', False))
        for p in peak_list.peaks
    ]
    corrected_pl = PeakList(
        scan_number=peak_list.scan_number,
        rt_min=peak_list.rt_min,
        peaks=corrected_peaks,
    )

    stage2 = _self_calibrate(
        corrected_pl, sequence,
        search_ppm=stage2_search_ppm,
        c_terminal_amide=c_terminal_amide,
        fixed_modification_da=fixed_modification_da,
        extra_residues=extra_residues,
        disulfide_bonds=disulfide_bonds,
    )

    def combined_correct(mz, _s1=stage1, _s2=stage2):
        return _s2(_s1(mz))

    return combined_correct


def compute_spectrum_coverage(
    result, peak_list, parent_mass, max_charge, tolerance_ppm,
    cal_correct_fn=None,
    parent_isotope_window=10,
    above_pct_threshold=1.0,
    n_top_unexplained=10,
    c_terminal_amide=False,
):
    if peak_list is None or peak_list.n_peaks == 0:
        return

    C13 = 1.003355
    PROTON = 1.00728
    correct = cal_correct_fn or (lambda mz: mz)

    raw_mz_arr = np.fromiter((p.mz        for p in peak_list.peaks), dtype=float,
                             count=peak_list.n_peaks)
    int_arr    = np.fromiter((p.intensity for p in peak_list.peaks), dtype=float,
                             count=peak_list.n_peaks)
    if cal_correct_fn is None:
        corr_mz_arr = raw_mz_arr
    else:
        corr_mz_arr = np.array([correct(mz) for mz in raw_mz_arr])

    sort_idx       = np.argsort(corr_mz_arr)
    corr_mz_sorted = corr_mz_arr[sort_idx]
    int_sorted     = int_arr[sort_idx]

    def _peaks_in_window_indices(target_mz, half_width):
        lo = np.searchsorted(corr_mz_sorted, target_mz - half_width, side='left')
        hi = np.searchsorted(corr_mz_sorted, target_mz + half_width, side='right')
        return lo, hi

    parent_iso_mz = []
    for z in range(1, max_charge + 1):
        m0 = (parent_mass + z * PROTON) / z
        for k in range(parent_isotope_window + 1):
            parent_iso_mz.append(m0 + k * C13 / z)
    parent_iso_arr = np.sort(np.array(parent_iso_mz))

    is_parent_arr = np.zeros(len(corr_mz_sorted), dtype=bool)
    if len(parent_iso_arr) > 0 and len(corr_mz_sorted) > 0:
        ins = np.searchsorted(parent_iso_arr, corr_mz_sorted)
        ins_lo = np.clip(ins - 1, 0, len(parent_iso_arr) - 1)
        ins_hi = np.clip(ins,     0, len(parent_iso_arr) - 1)
        d_lo = np.abs(corr_mz_sorted - parent_iso_arr[ins_lo])
        d_hi = np.abs(corr_mz_sorted - parent_iso_arr[ins_hi])
        d_min = np.minimum(d_lo, d_hi)
        safe_mz = np.where(corr_mz_sorted > 0, corr_mz_sorted, 1.0)
        ppm_min = d_min / safe_mz * 1e6
        is_parent_arr = ppm_min < tolerance_ppm

    total_int  = float(int_sorted.sum())
    parent_int = float(int_sorted[is_parent_arr].sum())

    matched_obs_mz = {round(m.observed_mz, 4) for m in result.matches}
    raw_mz_sorted  = raw_mz_arr[sort_idx]

    matched_int = sum(m.observed_intensity for m in result.matches)

    if matched_obs_mz:
        raw_rounded = np.round(raw_mz_sorted, 4)
        is_matched_arr = np.array(
            [bool(mz in matched_obs_mz) for mz in raw_rounded.tolist()],
            dtype=bool,
        )
    else:
        is_matched_arr = np.zeros(len(corr_mz_sorted), dtype=bool)

    unaccounted_mask = ~is_parent_arr & ~is_matched_arr
    unaccounted_mz   = raw_mz_sorted[unaccounted_mask]
    unaccounted_int  = int_sorted[unaccounted_mask]

    unexplained = max(0.0, total_int - parent_int - matched_int)

    result.total_intensity            = total_int
    result.parent_envelope_intensity  = parent_int
    result.matched_impurity_intensity = matched_int
    result.unexplained_intensity      = unexplained
    if total_int > 0:
        result.parent_envelope_fraction  = parent_int / total_int
        result.matched_impurity_fraction = matched_int / total_int
        result.unexplained_fraction      = unexplained / total_int

    base_int = peak_list.base_peak_int or 1.0

    n_iso_iemm  = 6
    iemm_native = 0.0
    iemm_deam   = 0.0
    iemm_charges_used = 0
    Tn = None
    if getattr(result, 'sequence', None):
        try:
            Tn = peptide_isotope_envelope(
                result.sequence, c_terminal_amide=c_terminal_amide,
                n_iso=n_iso_iemm,
            )
        except Exception:
            Tn = None
    if Tn is None:
        Tn = theoretical_isotope_envelope(parent_mass, n_iso=n_iso_iemm)
    Tn_arr = np.asarray(Tn, dtype=float)
    Td_arr = np.concatenate(([0.0], Tn_arr[:-1]))
    sum_Tn = float(Tn_arr.sum())
    sum_Td = float(Td_arr.sum())

    for z in range(1, max_charge + 1):
        m0_theo = (parent_mass + z * PROTON) / z
        spacing = C13 / z
        obs_env = np.zeros(n_iso_iemm + 1, dtype=float)
        for k in range(n_iso_iemm + 1):
            target_mz = m0_theo + k * spacing
            half_w    = target_mz * tolerance_ppm * 1e-6
            lo, hi    = _peaks_in_window_indices(target_mz, half_w)
            if lo < hi:
                obs_env[k] = int_sorted[lo:hi].max()
        envelope_total = float(obs_env.sum())
        if envelope_total < base_int * 0.001:
            continue
        a, b = _iemm_deconvolve_charge(obs_env, Tn_arr, Td_arr)
        native_z = a * sum_Tn
        deam_z   = b * sum_Td
        total_explained = native_z + deam_z
        if total_explained > envelope_total * 1.1 and total_explained > 0:
            scale = envelope_total / total_explained
            native_z *= scale
            deam_z   *= scale
        iemm_native += native_z
        iemm_deam   += deam_z
        iemm_charges_used += 1

    if iemm_charges_used > 0:
        iemm_total = iemm_native + iemm_deam
        if iemm_total > 0 and parent_int > 0:
            scale = min(1.0, parent_int / iemm_total)
            iemm_native *= scale
            iemm_deam   *= scale
        result.native_parent_intensity     = iemm_native
        result.deamidated_parent_intensity = iemm_deam
        if total_int > 0:
            result.native_parent_fraction      = iemm_native / total_int
            result.deamidated_parent_fraction  = iemm_deam   / total_int
        result.iemm_applied = True

    if len(unaccounted_int) > 0:
        n_take   = min(n_top_unexplained, len(unaccounted_int))
        top_idx  = np.argpartition(-unaccounted_int, n_take - 1)[:n_take]
        top_idx  = top_idx[np.argsort(-unaccounted_int[top_idx])]
        result.top_unexplained_peaks = [
            (round(float(unaccounted_mz[i]), 4),
             round(float(unaccounted_int[i]), 1),
             round(100.0 * float(unaccounted_int[i]) / base_int, 2))
            for i in top_idx
        ]
        result.n_unexplained_peaks_above_1pct = int(
            (unaccounted_int > base_int * above_pct_threshold / 100).sum()
        )
    else:
        result.top_unexplained_peaks = []
        result.n_unexplained_peaks_above_1pct = 0


def match_peaks(
    peak_list, sequence, tolerance_ppm=10.0, max_charge=3,
    intensity_threshold_fraction=0.001, c_terminal_amide=False,
    disulfide_bonds=0, fixed_modification_da=0.0, extra_residues=None,
    extra_residues_avg=None, coupling_reagent='HATU', apply_bear_case=True,
    synthesis_conditions=None, tfa_min_intensity_fraction=0.005,
    min_envelope_score=0.30,
    min_confidence=3.0,
    source_scans=None,
    min_persistence_fraction=0.30,
    cal_anchors=None,
    triage_mode=False,


):
    if triage_mode:
        min_envelope_score          = 0.0
        min_confidence              = 0.0
        min_persistence_fraction    = 0.0
        tfa_min_intensity_fraction  = 0.0
        intensity_threshold_fraction = min(intensity_threshold_fraction, 0.0001)

    if extra_residues is None:
        extra_residues = NON_STANDARD_RESIDUES

    seq         = sequence.upper().strip()
    parent_mass = calculate_monoisotopic_mass(
        seq, c_terminal_amide=c_terminal_amide, disulfide_bonds=disulfide_bonds,
        fixed_modification_da=fixed_modification_da, extra_residues=extra_residues,
    )

    effective_max_charge = max_charge
    bear_case_applied    = False
    if apply_bear_case:
        if parent_mass >= 6000.0:
            effective_max_charge = max(max_charge, 7)
            bear_case_applied = True
        elif parent_mass >= 3500.0:
            effective_max_charge = max(max_charge, 6)
            bear_case_applied = True
        elif parent_mass >= 2000.0:
            effective_max_charge = max(max_charge, 5)
            bear_case_applied = True
        elif parent_mass >= 1000.0:
            effective_max_charge = max(max_charge, 4)
            bear_case_applied = True

    from core.impurity_engine import enumerate_impurities
    imp_df = enumerate_impurities(
        seq, c_terminal_amide=c_terminal_amide, disulfide_bonds=disulfide_bonds,
        fixed_modification_da=fixed_modification_da, extra_residues=extra_residues,
        extra_residues_avg=extra_residues_avg if extra_residues_avg else extra_residues,
        synthesis_conditions=synthesis_conditions,
        tfa_min_intensity_fraction=tfa_min_intensity_fraction,
    )

    from core.delta_lookup import get_coupling_reagent_impurities
    from core.mass_calculator import calculate_average_mass
    from core.combination_impurities import enumerate_combination_impurities
    avg_extra  = extra_residues_avg if extra_residues_avg else {
        k: v * 1.0006 for k, v in extra_residues.items()
    }
    parent_avg = calculate_average_mass(seq, extra_residues=avg_extra)
    cr_df = get_coupling_reagent_impurities(parent_mass, parent_avg, coupling_reagent)
    if len(cr_df) > 0:
        imp_df = pd.concat([imp_df, cr_df], ignore_index=True)


    combo_df = enumerate_combination_impurities(
        seq, c_terminal_amide=c_terminal_amide,
        disulfide_bonds=disulfide_bonds,
        fixed_modification_da=fixed_modification_da,
        extra_residues=extra_residues,
    )
    if len(combo_df) > 0:
        imp_df = pd.concat([imp_df, combo_df], ignore_index=True).reset_index(drop=True)


    from core.impurity_engine import collapse_mass_degenerate
    imp_df = collapse_mass_degenerate(imp_df, mass_tolerance_da=0.001)

    if peak_list.n_peaks == 0:
        return MatchResult(
            sequence=seq, parent_mass=parent_mass, scan_number=peak_list.scan_number,
            n_peaks_input=0, n_peaks_matched=0, n_impurities_predicted=len(imp_df),
            bear_case_applied=bear_case_applied, tolerance_ppm=tolerance_ppm,
        )


    if cal_anchors:

        charges_in_anchors = set(
            h[0] for h in cal_anchors if len(h) >= 1 and h[0] is not None
        )
        _correct = _calibrate_from_anchors(cal_anchors, parent_mass)
        _test_mz_check = 1000.0
        if abs(_correct(_test_mz_check) - _test_mz_check) < 1e-6:

            _correct = _two_stage_calibrate(
                peak_list, seq,
                c_terminal_amide=c_terminal_amide,
                fixed_modification_da=fixed_modification_da,
                extra_residues=extra_residues,
                disulfide_bonds=disulfide_bonds,
            )
        elif len(charges_in_anchors) < 2:


            anchor_corrected_peaks = [
                Peak(mz=_correct(p.mz), intensity=p.intensity,
                     rt_min=p.rt_min, scan_number=p.scan_number,
                     charge_state=getattr(p, 'charge_state', None),
                     is_monoisotopic=getattr(p, 'is_monoisotopic', False))
                for p in peak_list.peaks
            ]
            anchor_corrected_pl = PeakList(
                scan_number=peak_list.scan_number,
                rt_min=peak_list.rt_min,
                peaks=anchor_corrected_peaks,
            )
            refine = _two_stage_calibrate(
                anchor_corrected_pl, seq,
                c_terminal_amide=c_terminal_amide,
                fixed_modification_da=fixed_modification_da,
                extra_residues=extra_residues,
                disulfide_bonds=disulfide_bonds,
            )
            _anchor_correct = _correct
            _correct = lambda mz, _a=_anchor_correct, _r=refine: _r(_a(mz))
    else:
        _correct = _two_stage_calibrate(
            peak_list, seq,
            c_terminal_amide=c_terminal_amide,
            fixed_modification_da=fixed_modification_da,
            extra_residues=extra_residues,
            disulfide_bonds=disulfide_bonds,
        )

    _test_mz = peak_list.peaks[0].mz if peak_list.peaks else 1000.0
    _delta = abs(_correct(_test_mz) - _test_mz)
    _peaks_already_corrected = False
    if _delta > _test_mz * 5e-6:
        _corr_peaks = [
            Peak(mz=round(_correct(p.mz), 5), intensity=p.intensity,
                 rt_min=p.rt_min, scan_number=p.scan_number,
                 charge_state=p.charge_state if hasattr(p, 'charge_state') else None,
                 is_monoisotopic=p.is_monoisotopic if hasattr(p, 'is_monoisotopic') else False)
            for p in peak_list.peaks
        ]
        peak_list = PeakList(
            scan_number=peak_list.scan_number,
            rt_min=peak_list.rt_min,
            peaks=_corr_peaks,
        )
        _peaks_already_corrected = True

    base_int      = peak_list.base_peak_int
    threshold     = base_int * intensity_threshold_fraction
    filtered      = [p for p in peak_list.peaks if p.intensity >= threshold]

    if not filtered:
        return MatchResult(
            sequence=seq, parent_mass=parent_mass, scan_number=peak_list.scan_number,
            n_peaks_input=peak_list.n_peaks, n_peaks_matched=0,
            n_impurities_predicted=len(imp_df),
            bear_case_applied=bear_case_applied, tolerance_ppm=tolerance_ppm,
        )


    _C13_ISO = 1.003355
    _parent_mz_list = []
    try:
        _parent_mono = calculate_monoisotopic_mass(seq)
        for _z in range(1, 8):
            _m0 = calculate_mz(_parent_mono, _z)
            for _iso in range(7):
                _parent_mz_list.append(_m0 + _iso * _C13_ISO / _z)
    except Exception:
        pass
    _parent_mz_np = np.sort(np.array(_parent_mz_list)) if _parent_mz_list else np.array([0.0])

    def _on_parent_cluster(obs_mz, tol_ppm=10.0):


        idx = np.searchsorted(_parent_mz_np, obs_mz)
        candidates = []
        if idx > 0:
            candidates.append(_parent_mz_np[idx - 1])
        if idx < len(_parent_mz_np):
            candidates.append(_parent_mz_np[idx])
        if not candidates:
            return False
        d_min = min(abs(c - obs_mz) for c in candidates)
        return bool(d_min / max(obs_mz, 1.0) * 1e6 < tol_ppm)

    candidates = []
    for _, row in imp_df.iterrows():


        if abs(float(row['delta_mono'])) < 1e-4:
            continue
        for z in range(1, effective_max_charge + 1):
            theo_mz = calculate_mz(row['mono_mass'], z)
            if 50 < theo_mz < 3000:
                candidates.append((row, z, theo_mz))


    risk_order = {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2}
    candidates.sort(key=lambda c: (risk_order.get(c[0]['risk_level'], 3), c[1]))

    peak_mz_arr  = np.array([p.mz for p in filtered])
    peak_int_arr = np.array([p.intensity for p in filtered])
    matched_peak_mask = np.zeros(len(filtered), dtype=bool)

    matches, matched_imp_charges = [], set()

    _C13_SPACING = 1.003355
    for row, z, theo_mz in candidates:
        imp_id = row['impurity_id']
        if (imp_id, z) in matched_imp_charges:
            continue

        _eff_tol = (min(tolerance_ppm, 5.0)
                    if row.get('impurity_type') == 'deamidation' and z == 1
                    else tolerance_ppm)

        iso_spacing = _C13_SPACING / z
        search_mzs = [theo_mz,
                      theo_mz + iso_spacing,
                      theo_mz + 2 * iso_spacing]

        best_local = None
        best_ppm_err = tolerance_ppm + 1.0
        found_at_iso = -1


        m0_abs_ppm = None
        for iso_idx, search_mz in enumerate(search_mzs):
            diffs_s = np.abs(peak_mz_arr - search_mz)
            ppm_errs_s = diffs_s / search_mz * 1e6
            in_tol_s = np.where((ppm_errs_s <= _eff_tol) & ~matched_peak_mask)[0]
            if len(in_tol_s) == 0:
                continue
            local_best = in_tol_s[np.argmin(diffs_s[in_tol_s])]
            local_ppm = float(ppm_errs_s[local_best])
            if iso_idx > 0:


                if m0_abs_ppm is None:
                    m0_abs_ppm = float(np.abs(peak_mz_arr - theo_mz).min() / theo_mz * 1e6)
                if m0_abs_ppm <= _eff_tol:


                    continue
            if local_ppm < best_ppm_err:
                best_ppm_err = local_ppm
                best_local = local_best
                found_at_iso = iso_idx
            if iso_idx == 0 and best_local is not None:
                break

        if best_local is None:
            continue

        ppm_err = best_ppm_err
        peak    = filtered[best_local]

        if row['impurity_type'] == 'tfa_adduct':
            if peak.intensity < base_int * tfa_min_intensity_fraction:
                continue

        iso   = check_isotope_pattern(filtered, peak.mz, z, row['mono_mass'], tolerance_ppm)


        delta_mono = float(row['delta_mono'])
        delta_iso_units = delta_mono / _C13_SPACING
        delta_int_alignment = abs(delta_iso_units - round(delta_iso_units))
        is_isotope_coincident = (
            abs(delta_mono) < 16.0 and
            delta_int_alignment < 0.05
        )
        on_parent = _on_parent_cluster(peak.mz, tol_ppm=tolerance_ppm)
        if is_isotope_coincident and on_parent and not iso:
            continue

        score = _score(ppm_err, row['risk_level'], peak.intensity, base_int, iso)

        if _on_parent_cluster(peak.mz):
            score = max(0.0, score - 2.0)

        env_score, _env_details = envelope_fit_score(
            filtered, peak.mz, z, row['mono_mass'],
            tolerance_ppm=tolerance_ppm, n_iso=4,
        )

        score = max(0.0, score - (1.0 - env_score) * 2.0)

        if abs(abs(float(row['delta_mono'])) - _C13_SPACING) < 0.05:
            _multi_z = any(
                m.impurity_type == row['impurity_type'] and m.charge_state != z
                for m in matches
            )
            if not _multi_z and env_score < 0.6:
                _r1, _ = theoretical_isotope_ratios(parent_mass)
                _pm0_mz = (parent_mass + z * PROTON) / z
                _pm0_diffs = np.abs(peak_mz_arr - _pm0_mz)
                _pm0_in = np.where(_pm0_diffs < _pm0_mz * tolerance_ppm * 1e-6)[0]
                if len(_pm0_in) > 0:
                    _pm0_int = float(peak_int_arr[_pm0_in[np.argmin(_pm0_diffs[_pm0_in])]])
                    if peak.intensity < 3.0 * _pm0_int * _r1:
                        continue
                else:
                    continue

        if env_score < min_envelope_score:
            continue
        if score < min_confidence:
            continue

        matches.append(PeakMatch(
            observed_mz=round(peak.mz, 5), observed_intensity=peak.intensity,
            theoretical_mz=round(theo_mz, 5), ppm_error=round(ppm_err, 4),
            charge_state=z, impurity_name=row['impurity_name'],
            impurity_type=row['impurity_type'], risk_level=row['risk_level'],
            confidence=score, isotope_confirmed=iso,
            delta_mono=row['delta_mono'], predicted_mass=row['mono_mass'],
            envelope_score=env_score,
            position=str(row.get('position', '')),
            n_degenerate=int(row.get('n_degenerate', 1)),
            notes=str(row.get('notes', '')),
        ))
        matched_peak_mask[best_local] = True
        matched_imp_charges.add((imp_id, z))


    low_risk_candidates = [(row, z, theo_mz) for row, z, theo_mz in candidates
                           if row['risk_level'] == 'LOW'
                           and (row['impurity_id'], z) not in matched_imp_charges]
    for row, z, theo_mz in low_risk_candidates:
        imp_id = row['impurity_id']
        if (imp_id, z) in matched_imp_charges:
            continue
        _eff_tol_lr = (min(tolerance_ppm, 5.0)
                       if row.get('impurity_type') == 'deamidation' and z == 1
                       else tolerance_ppm)
        _C13_sp = _C13_SPACING / z
        for search_mz in [theo_mz, theo_mz + _C13_sp, theo_mz + 2*_C13_sp]:
            diffs_s  = np.abs(peak_mz_arr - search_mz)
            ppm_s    = diffs_s / search_mz * 1e6
            in_tol_s = np.where((ppm_s <= _eff_tol_lr) & ~matched_peak_mask)[0]
            if len(in_tol_s) == 0:
                continue
            best_l = in_tol_s[np.argmin(diffs_s[in_tol_s])]
            ppm_e  = float(ppm_s[best_l])
            peak   = filtered[best_l]
            iso    = check_isotope_pattern(filtered, peak.mz, z, row['mono_mass'], tolerance_ppm)


            delta_mono = float(row['delta_mono'])
            delta_iso_units = delta_mono / _C13_SPACING
            delta_int_alignment = abs(delta_iso_units - round(delta_iso_units))
            is_isotope_coincident = (
                abs(delta_mono) < 16.0 and delta_int_alignment < 0.05
            )
            if is_isotope_coincident and _on_parent_cluster(peak.mz, tol_ppm=tolerance_ppm) and not iso:
                continue

            score  = _score(ppm_e, row['risk_level'], peak.intensity, base_int, iso)
            if _on_parent_cluster(peak.mz):
                score = max(0.0, score - 2.0)

            env_score, _ = envelope_fit_score(
                filtered, peak.mz, z, row['mono_mass'],
                tolerance_ppm=tolerance_ppm, n_iso=4,
            )
            score = max(0.0, score - (1.0 - env_score) * 2.0)

            if abs(abs(float(row['delta_mono'])) - _C13_SPACING) < 0.05:
                _multi_z_lr = any(
                    m.impurity_type == row['impurity_type'] and m.charge_state != z
                    for m in matches
                )
                if not _multi_z_lr and env_score < 0.6:
                    _r1, _ = theoretical_isotope_ratios(parent_mass)
                    _pm0_mz = (parent_mass + z * PROTON) / z
                    _pm0_diffs = np.abs(peak_mz_arr - _pm0_mz)
                    _pm0_in = np.where(_pm0_diffs < _pm0_mz * tolerance_ppm * 1e-6)[0]
                    if len(_pm0_in) > 0:
                        _pm0_int = float(peak_int_arr[_pm0_in[np.argmin(_pm0_diffs[_pm0_in])]])
                        if peak.intensity < 3.0 * _pm0_int * _r1:
                            continue
                    else:
                        continue

            if env_score < min_envelope_score:
                continue
            if score < max(2.0, min_confidence):
                continue
            matches.append(PeakMatch(
                observed_mz=round(peak.mz,5), observed_intensity=peak.intensity,
                theoretical_mz=round(theo_mz,5), ppm_error=round(ppm_e,4),
                charge_state=z, impurity_name=row['impurity_name'],
                impurity_type=row['impurity_type'], risk_level=row['risk_level'],
                confidence=score, isotope_confirmed=iso,
                delta_mono=row['delta_mono'], predicted_mass=row['mono_mass'],
                envelope_score=env_score,
                position=str(row.get('position', '')),
                n_degenerate=int(row.get('n_degenerate', 1)),
                notes=str(row.get('notes', '')),
            ))
            matched_peak_mask[best_l] = True
            matched_imp_charges.add((imp_id, z))
            break

    matched_obs = {round(filtered[i].mz, 3) for i in np.where(matched_peak_mask)[0]}

    unmatched_peaks = [round(p.mz, 5) for p in filtered if round(p.mz, 3) not in matched_obs]

    _matched_imp_id_set = {imp_id for imp_id, _ in matched_imp_charges}
    matched_types_set = {row['impurity_type'] for _, row in imp_df.iterrows()
                         if row['impurity_id'] in _matched_imp_id_set}
    unmatched_imps = list({row['impurity_type'] for _, row in imp_df.iterrows()
                           if row['impurity_type'] not in matched_types_set})


    if source_scans and len(source_scans) > 1:
        kept = []
        for m in matches:
            target_mz = m.observed_mz
            n_present = 0
            for spec in source_scans:
                if spec.is_empty or len(spec.mz_array) == 0:
                    continue

                diffs_ppm = np.abs(spec.mz_array - target_mz) / target_mz * 1e6
                if diffs_ppm.min() < tolerance_ppm:
                    n_present += 1
            persistence = n_present / len(source_scans)
            m.scan_persistence = round(persistence, 3)
            if persistence >= min_persistence_fraction:
                kept.append(m)
        matches = kept


    if triage_mode:
        base_int = peak_list.base_peak_int or 1.0
        for m in matches:
            persistence_ok = (m.scan_persistence is None or m.scan_persistence >= 0.5
                              or not source_scans or len(source_scans) <= 1)
            if (m.envelope_score >= 0.65 and m.confidence >= 5.0
                    and m.isotope_confirmed and persistence_ok):
                m.triage_tier = 'A'
            elif m.envelope_score >= 0.30 and m.confidence >= 3.0:
                m.triage_tier = 'B'
            else:
                m.triage_tier = 'C'
        tier_order = {'A': 0, 'B': 1, 'C': 2}
        matches.sort(key=lambda m: (tier_order.get(m.triage_tier, 3),
                                     -(m.observed_intensity / base_int)))

    _bint = peak_list.base_peak_int or 1.0
    for _m in matches:
        if not _m.isotope_confirmed and (100.0 * _m.observed_intensity / _bint) >= 5.0:
            _m.iso_warning = True
            if _m.risk_level == 'HIGH':
                _m.risk_level = 'MEDIUM'

    result = MatchResult(
        sequence=seq, parent_mass=parent_mass, scan_number=peak_list.scan_number,
        n_peaks_input=len(filtered), n_peaks_matched=len(matches),
        n_impurities_predicted=len(imp_df), matches=matches,
        unmatched_peaks=unmatched_peaks, unmatched_impurities=unmatched_imps,
        bear_case_applied=bear_case_applied, tolerance_ppm=tolerance_ppm,
    )
    compute_spectrum_coverage(
        result, peak_list, parent_mass, effective_max_charge, tolerance_ppm,
        cal_correct_fn=(lambda mz: mz) if _peaks_already_corrected else _correct,
        c_terminal_amide=c_terminal_amide,
    )
    return result


